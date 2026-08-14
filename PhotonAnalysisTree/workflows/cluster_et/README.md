# Pythia minimum-bias cluster-ET map-reduce workflow

This directory is the complete cluster-ET analysis workflow:

- `Fun4All_PythiaClusterEtSpectra.C`: map synchronized DST inputs to partial histograms;
- `check_pythia_cluster_et_partial.C`: validate each partial;
- `FinalizePythiaClusterEtSpectra.C`: reduce partials and write ROOT/PDF outputs;
- `run_partial.sh` and `submit.job`: execute and submit the map jobs.

`PythiaClusterEtSpectrum` reads the synchronized `DST_CALO_CLUSTER`,
`DST_MBD_EPD`, `DST_TRUTH_JET`, and `G4Hits` streams directly. It writes
eight raw SPLIT-cluster histograms (two inclusive, three geometrical topology,
and three energy-deposit topology) plus one metadata entry; it does not write a
large per-cluster tree.

## Cluster selections

The default selection requires reconstructed cluster energy above 0.2 GeV,
`|eta_cluster| < 0.7`, valid truth provenance, signal embedding ID 1, and a
dominant truth-contribution fraction of at least 0.5.

- Prompt clusters have a dominant G4-primary photon whose HepMC classifier
  category is direct (1) or fragmentation (2), with `|eta_truth| < 0.7`.
- Pi0-origin clusters include two production pathways. A HepMC pi0 transported
  into Geant4 is identified through its matching G4 primary. A pi0 decayed by
  the generator is identified through the G4-primary photon and its HepMC
  photon-copy chain. The eta cut is applied to the parent pi0. G4-secondary
  pi0s produced during detector transport are excluded from both the inclusive
  pi0-origin histogram and all topology classifications.

Pathway-specific counters are retained so that the Detroit-production
assumption can be checked rather than imposed.

## Geometry and topology

For every selected G4-primary or generator-decayed signal pi0 with exactly two
photon daughters, each G4 photon is
projected from its own `PHG4VtxPoint` along its momentum to the CEMC cylinder.
The projected surface `(eta,phi)` is compared with the RawCluster centroid
expressed in the same detector-origin coordinate. Only clusters containing an
ancestry-compatible contribution from that selected pi0 above the configured
fraction are candidates.

The exclusive priority is:

1. `separated`: both photons match different clusters within `DeltaR < 0.03`;
2. `merged`: one cluster lies within `DeltaR < 0.06` of both projections and
   has `0.5 <= ET_cluster/pT_pi0 <= 1.5`;
3. `missing`: exactly one photon matches and its cluster response to that
   daughter is in `[0.5,1.5]`;
4. `none` or `ambiguous`: diagnostic counters only.

## Energy-deposit topology

The `G4Hits` production stream contains `G4HIT_CEMC`, but not the intermediate
`G4CELL_CEMC` or legacy raw truth towers needed by the direct-daughter matcher.
The map macro reconstructs those two transient nodes with the official
`PHG4FullProjSpacalCellReco` and `RawTowerBuilder` chain before running the
analysis. They are not written to the partial output.

Candidate clusters require the summed ancestry-compatible pi0 contributor
fraction to be at least 0.5. For each candidate cluster, CEMC hit deposits are
traced through the reconstructed cell and tower maps to the two direct pi0
daughter track IDs. Each daughter selects the cluster with the largest positive
deposit; the optional minimum fraction of the cluster truth deposit defaults to
zero.

The exclusive energy-deposit priority is:

1. `separated`: the two daughters select two different clusters;
2. `merged`: both select the same cluster and
   `0.5 <= ET_cluster/pT_pi0 <= 1.5`;
3. `missing`: otherwise, at least one selected cluster has
   `0.5 <= ET_cluster/pT_gamma <= 1.5`;
4. `none`: diagnostic counter only.

Separated pi0s fill two cluster entries. Merged and missing pi0s fill one.
Pi0s with no matched cluster cannot be represented on a cluster-ET axis.

## Run and reduce

The Condor defaults process all 200,000 manifest files in chunks of 50:

```bash
condor_submit -maxjobs 4000 workflows/cluster_et/submit.job
```

Each wrapper writes transactionally and runs
`check_pythia_cluster_et_partial.C` before publishing a partial. Reduce the
complete contiguous range with:

```bash
root -l -b -q \
  'workflows/cluster_et/FinalizePythiaClusterEtSpectra.C("output/cluster_et_partial/prompt_primary_generator_pi0_eta07_energy_contribution_0p3/partial_*.root","output/plots/cluster_et_prompt_primary_generator_pi0_eta07",0,200000)'
```


The final ROOT file contains all eight raw-count histograms and their
bin-width-normalized density copies. Density is computed only after summing all
partials, so its unit is `clusters/GeV`. The finalizer writes two raw-count
plots:

- `<output_base>.pdf`: geometrical matching;
- `<output_base>_energy_contribution.pdf`: energy-deposit matching.

Counts are unweighted, matching the policy used by
`minbias_truth_pt_prompt_eta07.pdf`.
