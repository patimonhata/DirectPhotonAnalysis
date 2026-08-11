# Pythia truth spectrum tree schema

Schema version: 1. This compact tree is produced by
`PythiaTruthSpectrumTree` from the `G4Hits` DST only. It is independent of
reconstructed clusters and is intended for generator/G4 truth particle-yield
studies.

## Event identity and weights

- `event_uid = (source_file_id << 32) | event_in_file`.
- `hepmc_event_number` is taken from embedding ID 1 in
  `PHHepMCGenEventMap`.
- `event_weights` retains the complete HepMC weight vector.
- `event_weight` is its first element, or 1 when the vector is empty.
- `event_weight_valid` requires every stored HepMC weight to be finite.

Raw particle counts should use unit weight. Set `use_event_weight=true` in the
plot macro only when the production weight convention has been validated.

## Final HepMC photons

`truth_photon_*` vectors are mutually aligned and contain only particles with
PDG 22, status 1, and no end vertex. The collection therefore counts one final
record rather than every generator bookkeeping copy.

Kinematics are stored in `e`, `pt`, `eta`, and `phi`, with an explicit
`kinematics_valid` flag. `barcode` and `status` retain record identity.

`HepMCPhotonClassifier` runs on each final photon. Its aligned output is:

- `classification_valid`
- `category`: 0 unclassified, 1 direct 2-to-2, 2 fragmentation, 3 decay
- `immediate_parent_count`, `immediate_parent_pdg`

The dumper separately follows every unique photon parent backwards. It stores
`copy_chain_valid`, `copy_depth`, and the parent count at the first
non-photon-copy production vertex. When that vertex has exactly one parent,
`origin_parent_pdg` and `origin_parent_barcode` identify it. This prevents an
intermediate photon-to-photon copy from hiding the physical origin.

The direct-photon spectrum uses final HepMC photon `pt` with valid classifier
category 1. It does not infer directness from G4 PDG alone.

## HepMC pi0

`truth_pi0_*` vectors contain PDG 111 records that have no PDG 111 daughter.
This retains only the last record in a possible pi0-to-pi0 bookkeeping chain.
The collection stores `barcode`, `status`, `kinematics_valid`, `e`, `pt`,
`eta`, and `phi`.

`hepmc_direct_photon_count` is a production diagnostic counting direct photon
daughters attached to that HepMC pi0 record. In the checked Detroit
minimum-bias production the pi0 records are stable (`status=1`) in HepMC and
their decay is delegated to Geant4, so this value is normally zero.

## G4 pi0 decay photons

The decay-photon spectrum cannot be obtained from HepMC daughters in this
production. `truth_pi0_decay_photon_*` therefore selects a G4 secondary only
when all of the following hold:

1. the secondary has PDG 22;
2. its immediate G4 parent has PDG 111;
3. that parent is a G4 primary in embedding ID 1;
4. the parent barcode resolves to a PDG 111 particle in the signal HepMC event.

The collection stores the photon G4 track ID, parent G4 track ID, parent HepMC
barcode, validity, and photon `e`, `pt`, `eta`, and `phi`. Requiring the
immediate parent excludes the many later shower photons that still have a pi0
ancestor. The plotted decay-photon x coordinate is the G4 daughter photon's
production `pt`.

## Metadata diagnostics

The one-entry `metadata` tree records the node names, selection strings,
classifier version, event counts, and aggregate counts for HepMC records,
final photons, terminal pi0s, and G4 pi0 decay photons. It also records photon
copy edges, nonterminal pi0 copies, invalid ancestry/kinematics, and unmatched
G4 candidates.

Use `macro/check_pythia_truth_spectrum_tree.C` to validate the schema, vector
alignment, unique particle IDs, aggregate counts, and G4-parent/HepMC-pi0
barcode consistency.
