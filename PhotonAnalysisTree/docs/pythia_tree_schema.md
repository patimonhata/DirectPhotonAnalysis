# Pythia tree schema

Schema version: 4. This is produced by `PythiaPhotonAnalysisTree`; the
single-particle `PhotonAnalysisTree` remains schema version 3.

## Event and vertex

- `event_uid = (source_file_id << 32) | event_in_file`.
- `hepmc_n_subevent` is `PHHepMCGenEventMap::size()` and
  `hepmc_embedding_id` stores every map key in ascending order.
- The reconstruction vertex is `PHHepMCGenEventMap::get(signal_embedding_id)`
  followed by `PHHepMCGenEvent::get_collision_vertex()`. Its x/y/z unit is cm.
- `vertex_embedding_id` records which subevent supplied that vertex. The
  default signal embedding ID is 1; a missing, non-simulated, or invalid signal
  event is rejected instead of silently using `(0,0,0)`.

## Reconstructed quantities shared with single-particle trees

All ordinary `split_*` and `nosplit_*` cluster, pair, shower-shape, and tower
branches have the same definitions and ordering as schema version 3. They are
filled by the shared `PhotonTreeCommon` implementation. Cluster direction and
pair mass use the HepMC collision vertex above.

## Cluster truth matching

For `<c> = split` or `nosplit`, every ordinary cluster has aligned summary
vectors named `<c>_cluster_truth_*`:

- `valid`, `total_edep`, `n_contributor`
- `dominant_g4_track_id`, `dominant_g4_pdg_id`, `dominant_embedding_id`
- `dominant_hepmc_barcode`, `dominant_edep`, `dominant_fraction`
- `dominant_hepmc_valid`, `dominant_hepmc_pdg_id`
- `dominant_photon_category`, `dominant_photon_source`
- `dominant_immediate_parent_pdg`, `dominant_classification_parent_pdg`

All contributors are retained in a flat table. `contributor_offset` has
`ncluster+1` entries; contributor indices for cluster `i` occupy
`[offset[i], offset[i+1])`. The remaining `contributor_*` vectors are mutually
aligned and include `contributor_cluster_index` explicitly.

Contributors are sorted by deposited energy descending and then G4 track ID
ascending. `dominant_*` is the first contributor. Fractions are normalized to
the sum of matched G4 shower edep in that cluster, not to reconstructed cluster
energy.

### Matching chain

1. Iterate every RawCluster tower entry.
2. Read every `(G4 shower ID, edep)` entry in
   `TowerInfo::get_showerEdepMap()`. If that map was stripped in an older DST
   layout, use the corresponding legacy `RawTower::get_g4showers()` map.
   A nonempty cluster with neither provenance map is invalid rather than a
   The event is retained with that cluster's truth `valid=0`, and the metadata
   `n_events_invalid_truth` counter is incremented.
   valid zero-contributor match. This is not a first-hit-only match.
3. Resolve each shower with `PHG4TruthInfoContainer::GetShower()`.
4. Use the official `CaloTruthEval::get_primary_particle(PHG4Shower*)` to
   follow that shower to its G4 primary. Contributions from every descendant
   shower/hit represented in the tower map are accumulated by primary track ID.
5. Use `CaloTruthEval::get_embed(primary)` to select the matching
   `PHHepMCGenEventMap` value, then resolve `primary->get_barcode()` with
   `HepMC::GenEvent::barcode_to_particle()`.

For NO_SPLIT clusters the full tower shower edep is counted. A tower may be
shared by SPLIT clusters, so SPLIT contributions are multiplied by
`clamp(RawCluster tower value / calibrated TowerInfo energy, 0, 1)`. This
proportional allocation is an analysis choice and is recorded in metadata; it
is not supplied by the sPHENIX evaluator.

## Photon classification

The HepMC classifier is versioned separately in metadata. Categories reproduce
the rules in the referenced `PhotonAna::photon_type` implementation:

- `-1`: non-photon
- `0`: photon, unclassified
- `1`: direct (2-to-2 rule)
- `2`: fragmentation (1-to-2 rule retaining the incoming parton)
- `3`: decay (single incoming parent with `abs(PDG) > 37`)

It first follows unique photon-to-photon bookkeeping copies backwards.
`photon_source` is `-1` for non-photon, `1` for a photon parent, `2` for pi0,
`3` for eta, and `0` for another/no parent. These physics categories are local
analysis policy, not an official sPHENIX classification.

## Official framework code versus local policy

sPHENIX provides the DST classes and provenance primitives:
`PHHepMCGenEventMap`, `PHHepMCGenEvent`, `TowerInfo`/`RawTower` shower edep maps,
`PHG4TruthInfoContainer`, and `CaloTruthEval` shower-to-primary/embed lookup.

This package implements the cluster tower traversal, SPLIT allocation,
per-primary accumulation and sorting, HepMC photon category rules, branch
layout, validity policy, and metadata. The distinction is encoded by the
`truth_scheme`, `truth_matcher_version`, `photon_classifier_version`, and
allocation fields in the metadata tree.

Use `macro/check_pythia_tree.C` to validate vector alignment, offsets, edep
sums, contributor fractions, and schema metadata.
