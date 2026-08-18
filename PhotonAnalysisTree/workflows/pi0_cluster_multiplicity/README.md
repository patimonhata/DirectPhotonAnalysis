# Pythia primary/generator pi0-compatible cluster multiplicity workflow

This workflow measures, for each selected truth pi0, how many SPLIT CEMC
clusters satisfy

    f_pi0^cluster = sum_{contributors compatible with pi0} f_contributor.

It is independent of `workflows/cluster_et` and writes to a separate partial
output directory. The four thresholds are evaluated in one DST pass:

- `fmin=0.0` means strictly `f_pi0 > 0`, not `f_pi0 >= 0`;
- `fmin=0.1`, `0.3`, and `0.5` use `f_pi0 >= fmin`.

## Selection

The selected truth pi0 population contains only transported G4-primary pi0
decays and generator-level pi0 decays represented by their two G4-primary
photons. Pi0s created as G4 secondaries during detector transport are excluded.

The defaults require `|eta_truth_pi0| < 0.7`. No cluster-eta acceptance cut is
applied: every cluster in the SPLIT CEMC container can be considered. Clusters
must have finite energy and position and valid truth provenance. There is
deliberately no analysis-level cluster-energy threshold. The RawCluster
reconstruction may still have its own intrinsic thresholds.

`PythiaClusterTruthMatcher` supplies primary-shower energy fractions from the
TowerInfo shower provenance. This workflow does not use the direct-daughter
G4-cell/hit matcher and does not reconstruct transient CEMC cells or truth
towers. The same four synchronized DST streams are registered to preserve the
established minimum-bias input synchronization.

## Outputs

Each partial contains:

- the pi0-by-pi0 cluster multiplicity for all four thresholds;
- multiplicity versus truth pi0 pT;
- multiplicity split into G4-primary and generator pathways;
- the maximum and second-largest compatible fraction per pi0;
- compatible fraction versus reconstructed cluster energy;
- one metadata entry with selection definitions and counters.

The finalizer writes raw and probability histograms, a `threshold_summary` tree
containing the mean, `P(N>=2)`, `P(N>=3)`, and overflow probability, plus:

- `<output_base>.pdf`: four overlaid multiplicity probability distributions;
- `<output_base>_vs_truth_pt.pdf`: the four multiplicity-versus-pT maps;
- `<output_base>_fraction_vs_cluster_energy.pdf`: low-energy-fragment
  diagnostic.

## Build and pilot

Build and install `PhotonAnalysisTree` after the currently running production
no longer depends on changing the installed library:

```bash
PhotonAnalysisTree/src/build.sh
```

`submit.job` intentionally defaults to a 100-file pilot in ten jobs. Review its
paths and then submit manually when ready:

```bash
condor_submit -maxjobs 10 workflows/pi0_cluster_multiplicity/submit.job
```

No submission is performed by repository scripts automatically. Each map job
writes transactionally and validates its partial before publishing it.

Reduce the complete pilot range with:

```bash
root -l -b -q \
  'workflows/pi0_cluster_multiplicity/FinalizePythiaPi0ClusterMultiplicity.C("output/pi0_cluster_multiplicity_partial/pilot_primary_generator_truth_eta07_no_cluster_eta_or_energy_cut/partial_*.root","output/plots/pi0_cluster_multiplicity_primary_generator_truth_eta07_no_cluster_eta_or_energy_cut",0,100)'
```

For full production, change `total_files`, `files_per_job`, `n_chunks`, and the
partial output directory together. Keep
`n_chunks = ceil(total_files/files_per_job) - chunk_offset`.
