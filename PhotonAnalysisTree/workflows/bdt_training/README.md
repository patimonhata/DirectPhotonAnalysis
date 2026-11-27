# BDT training clusters

`ExtractBdtTrainingDataset.C` reads the schema-11 candidate maps and writes a `clusters` TTree with one row per split cluster satisfying all three strict cuts, `abs(mbd_vertex_z) < 60` cm, `abs(split_cluster_eta) < 0.7`, and `3 < split_cluster_et < 35` GeV. This includes clusters with unavailable shower shapes and clusters excluded from every training version.

## Run

From the repository root, run the following in the analysis ROOT environment to write one file covering all seven Jet and four PhotonJet subsamples.

```sh
root -l -b -q 'PhotonAnalysisTree/workflows/bdt_training/ExtractBdtTrainingDataset.C("PhotonAnalysisTree/output/photon_candidate_selection/productions/20261125_ana565/jet_cluster0p12_partner0p5/maps","PhotonAnalysisTree/output/photon_candidate_selection/productions/20261126_ana565/photonjet_cluster0p12_partner0p5/maps","PhotonAnalysisTree/output/bdt_training/training.root")'
```

To produce individual samples or parallel shards, use the sample entry point.

```sh
root -l -b -q 'PhotonAnalysisTree/workflows/bdt_training/ExtractBdtTrainingClusters.C("PhotonAnalysisTree/output/photon_candidate_selection/productions/20261125_ana565/jet_cluster0p12_partner0p5/maps","jet5","PhotonAnalysisTree/output/bdt_training/jet5.root")'
root -l -b -q 'PhotonAnalysisTree/workflows/bdt_training/ExtractBdtTrainingClusters.C("PhotonAnalysisTree/output/photon_candidate_selection/productions/20261126_ana565/photonjet_cluster0p12_partner0p5/maps","photonjet5","PhotonAnalysisTree/output/bdt_training/photonjet5.root")'
```

Use the corresponding sample name for each Jet and PhotonJet subsample. The optional fourth and fifth arguments are a zero-based shard index and shard count (defaults 0 and 1). Each shard reads a disjoint range of maps after validating metadata and computing the weight denominator over the complete subsample. Use a different output filename per shard. Output is published after a successful close and requires a fresh file path. A ROOT `TChain("clusters")` can read all outputs directly, or use `hadd` with each subsample's shards exactly once. The `metadata` tree records the input map root, schema, sample, normalization and half-open map range for each shard.

The extraction checks the truth isolation cone radius in every PhotonJet map of each input subsample, including maps outside the requested shard. All PhotonJet inputs to one dataset must use the same radius. The output metadata records `truth_isolation_radius` per sample and shard. Earlier schema-11 maps with `truth_isolation_definition` and no explicit radius use 0.4. Jet metadata records −1 because its training labels use no truth isolation.

## Columns and labels

- `sample_name`, `event_uid`, `cluster_id` identify the original cluster.
- `event_eligible` requires a valid event weight and valid, passing sample stitching. For eligible events `weight_pb` divides the event’s `weight_numerator_pb` by the sum of `sum_generator_weight_processed` over all maps of its subsample, before vertex selection. Ineligible events remain available with zero physical weight and three excluded labels.
- Inputs in the existing model order are `cluster_Et`, `cluster_weta_cogx`, `cluster_wphi_cogx`, `vertexz`, `cluster_Eta`, `e11_over_e33`, `cluster_et1`, `cluster_et2`, `cluster_et3`, `cluster_et4`, `e32_over_e35`. The vertex and cluster kinematics use MBD z. `inputs_valid` requires a valid shower shape and finite inputs.
- `common_main` and `common_origin` store the map's common composition class and overview (0 prompt, 1 π⁰, 2 η, 3 Other). `topology_applicability` keeps the map code. `topology` is −1 unless the corresponding π⁰ or direct η→γγ anchor is valid and evaluable. The other codes are 0 other, 1 separated, 2 merged, 3 missing and 4 single_contaminated.
- `truth_iso_valid` and `truth_pass_isolated` use the PhotonJet map flags. Both are false for Jet.

| Eligible sample and category | label_v1 | label_v2 | label_v3 |
|---|---:|---:|---:|
| PhotonJet isolated prompt with valid truth isolation | 1 | 1 | 1 |
| Jet π⁰/η missing | 1 | 1 | 1 |
| Jet π⁰/η merged or single_contaminated | 0 | 0 | 0 |
| Jet decay_photon | −1 | 1 | 1 |
| Jet π⁰/η separated | −1 | −1 | 1 |
| Jet common Other excluding decay_photon | 0 | 0 | 0 |
| All remaining clusters | −1 | −1 | −1 |

For version 1, `inputs_valid && label_v1 >= 0` selects usable labeled rows. Use `label_v2` or `label_v3` for the other versions. The training framework can group rows by `(sample_name, event_uid)` for event-level partitioning. Use the 11 listed input columns for scoring.

To cross-check an output against its source maps, run the validation script with uproot, Awkward Array and NumPy.

```sh
python3 PhotonAnalysisTree/validation/bdt_training/check_training_shard.py PhotonAnalysisTree/output/bdt_training/training.root
```
