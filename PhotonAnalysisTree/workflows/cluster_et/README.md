# Pythia minimum-bias cluster-ET map-reduce workflow

This directory is the complete cluster-ET analysis workflow:

- `Fun4All_PythiaClusterEtSpectra.C`: map synchronized DST inputs to partial histograms;
- `check_pythia_cluster_et_partial.C`: validate each partial;
- `FinalizePythiaClusterEtSpectra.C`: reduce partials and write ROOT/PDF outputs;
- `run_partial.sh` and `submit.job`: execute and submit the map jobs.

`PythiaClusterEtSpectrum` reads the synchronized `DST_CALO_CLUSTER`,
`DST_MBD_EPD`, `DST_TRUTH_JET`, and `G4Hits` streams directly. It writes five
raw SPLIT-cluster histograms plus one metadata entry; it does not write a large
per-cluster tree.

## Cluster selections

The default selection requires reconstructed cluster energy above 0.2 GeV,
`|eta_cluster| < 0.7`, valid truth provenance, signal embedding ID 1, and a
dominant truth-contribution fraction of at least 0.5.

- Prompt clusters have a dominant G4-primary photon whose HepMC classifier
  category is direct (1) or fragmentation (2), with `|eta_truth| < 0.7`.
- Pi0-origin clusters include three production pathways. A HepMC pi0 transported
  into Geant4 is identified through its matching G4 primary. A pi0 decayed by
  the generator is identified through the G4-primary photon and its HepMC
  photon-copy chain. A secondary pi0 produced after a transported long-lived
  particle travels is identified by tracing that pi0 back to the contributing
  signal G4 primary and then using its own two G4 photon daughters. The eta cut
  is applied to the parent pi0. The reduced truth record has no creator-process
  label, so the secondary category includes both displaced decays and secondary
  pi0 production in material; this limitation is retained in the metadata
  selection string and should be treated as a systematic choice.

Pathway-specific counters are retained so that the Detroit-production
assumption can be checked rather than imposed.

## Geometry and topology

For every signal pi0 with exactly two photon daughters, including secondary
pi0s from displaced parent decays, each G4 photon is
projected from its own `PHG4VtxPoint` along its momentum to the CEMC cylinder.
This preserves displaced decays. The projected surface `(eta,phi)` is compared
with the RawCluster centroid expressed in the same detector-origin coordinate.
Only clusters containing an ancestry-compatible pi0 contribution above the
configured fraction are candidates.

The exclusive priority is:

1. `separated`: both photons match different clusters within `DeltaR < 0.03`;
2. `merged`: one cluster lies within `DeltaR < 0.06` of both projections and
   has `0.5 <= ET_cluster/pT_pi0 <= 1.5`;
3. `missing`: exactly one photon matches and its cluster response to that
   daughter is in `[0.5,1.5]`;
4. `none` or `ambiguous`: diagnostic counters only.

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
  'workflows/cluster_et/FinalizePythiaClusterEtSpectra.C("output/cluster_et_partial/prompt_pi0_eta07/partial_*.root","output/plots/minbias_cluster_et_prompt_pi0_eta07",0,200000)'
```


The final ROOT file contains both raw counts and density histograms. Density is
computed only after summing all partials by dividing each bin by its width, so
its unit is `clusters/GeV`. The PDF plots the summed raw histograms with a
`counts/bin` y-axis. Counts are unweighted, matching the policy used by
`minbias_truth_pt_prompt_eta07.pdf`.
