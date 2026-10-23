# Pythia photon-candidate selection

Main map → reduce → merge workflow. Current production is map schema **8**, pi0 topology **12**, eta topology **2**, supplemental decay-origin **1**, and partial schema **10**. Physics definitions and diagnostic details are in [definition notes](../../docs/photon_candidate_definitions.md).

## Results and inputs

```text
output/photon_candidate_selection/
  current/jet/  -> validated current results
  productions/
    topology_v12_eta_v2_decay_v1/ClusterE0p12_Pi0Partner0p5/
      maps/<sample>/map_<chunk>.root
      partials/jet/<sample>/shard_<index>/photon_candidate_selection.root
      results/jet/
    topology_v11/<configuration>/  # retained schema-5 reference
  qa/
```

The existing map directories link to cold storage. Do not remove the underlying map inputs when reorganizing diagnostics. Saved ROOT metadata retains its original paths; `output/layout_manifest.json` records where files moved. See [layout](../../docs/layout.md).

Each map reads synchronized `DST_CALO_CLUSTER`, `DST_MBD_EPD`, `DST_TRUTH_JET`, and `G4Hits` basenames from `input/<sample>/segments.list`. Fun4All resolves the logical filenames through the file catalog. Jet submit defaults use 50 segments per map, 200 maps per ordinary sample and 2,000 for Jet12. Exact selection settings are in each submit file and map metadata. Jet defaults are cluster E>0.12 GeV, pi0/eta partner E>0.5 GeV, pi0 0.10<m<0.20 GeV, eta 0.45<m<0.65 GeV, missing band 0.12<E<=0.5 GeV, and recovery disabled. PhotonJet defaults are a separate configuration; do not mix them with Jet outputs.

## Build

```bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
cmake -S PhotonAnalysisTree/src -B PhotonAnalysisTree/build/candidate \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="$PWD/PhotonAnalysisTree/install"
cmake --build PhotonAnalysisTree/build/candidate --parallel 4
cmake --install PhotonAnalysisTree/build/candidate
```

Finish jobs using the installed library before installing another build.

## Map

Review the production/configuration and choose new map outputs when rerunning a production. Existing map files are not overwritten. No repository script automatically submits jobs.

```bash
mkdir -p PhotonAnalysisTree/output/condor/photon_candidate_selection
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet5.job
```

Other samples have matching submit files. The explicit local interface is:

```text
run_map.sh JOB_INDEX CHUNK_OFFSET TOTAL_FILES FILES_PER_JOB INPUT_MANIFEST SAMPLE OUTPUT_DIRECTORY [N_EVENTS] [MIN_CLUSTER_ENERGY] [PI0_PARTNER_MIN_ENERGY] [ETA_PARTNER_MIN_ENERGY] [PI0_MASS_MIN] [PI0_MASS_MAX] [ETA_MASS_MIN] [ETA_MASS_MAX] [MISSING_ENERGY_MIN] [MISSING_ENERGY_MAX]
```

[`check_sample_map_outputs.sh`](../../validation/photon_candidate_selection/check_sample_map_outputs.sh) checks expected file coverage and optionally performs deep ROOT validation; pass the output directory and exact cuts for a non-default production. [`run_small_sample.sh`](../../validation/photon_candidate_selection/run_small_sample.sh) runs limited-input end-to-end QA under `output/photon_candidate_selection/qa/`; this is not a physics normalization. Jet12 QA needs at least ten maps so every fixed shard is nonempty.

## Reduce and merge

One reducer fills all six selections: `kinematic`, `preselection`, `preselection_tight`, `preselection_isolation`, `region_a`, and `region_a_tagging_veto`. The last is Region A after the pi0-or-eta veto. The input maps also retain Region B/C/D and gap flags for later ABCD studies.

Use full-sample generator-weight normalization, including processed events later rejected by the vertex cut. Never normalize individual map files or add per-shard normalization denominators. Jet12 has ten contiguous shards; each other Jet sample has one. Schema/settings, exact sample/shard coverage, category closure, weighted Sumw2 and subset relations are checked before final plots are produced.

```bash
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_jet_samples.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_jet12_shards.job
```

Interfaces:

```text
run_reduce.sh FAMILY MAP_ROOT OUTPUT_BASE REQUIRE_COMPLETE N_BINS ET_MAX SAMPLE SHARD_INDEX
run_merge.sh FAMILY PARTIAL_ROOT OUTPUT_DIRECTORY
```

After all 16 partials complete, use a new final output directory:

```bash
bash PhotonAnalysisTree/workflows/photon_candidate_selection/run_merge.sh jet \
  PhotonAnalysisTree/output/photon_candidate_selection/productions/topology_v12_eta_v2_decay_v1/ClusterE0p12_Pi0Partner0p5/partials/jet \
  /tmp/photon_candidate_merge_replay
```

Merged outputs include candidate composition, pi0/eta anchor topology, eta cascade composition, Other-origin and supplemental decay-origin histograms, spectra/fractions/survival PDFs, ET-binned BDT distributions and summaries. ROOT histogram metadata records source map/partial versions, cuts, weights and uncertainty conventions. Statistical errors use the documented candidate/pair-level Sumw2 and subset covariance; event-level correlations are not estimated.

## Retained v11/schema-5 results

Four configurations under `productions/topology_v11/` retain maps, 16 partials and final results. Partner thresholds are 0.15 or 0.5 GeV and pi0 mass windows are 0.10–0.20 or 0.08–0.18 GeV. Eta topology extension is absent, but eta tagging remains enabled. Use `run_reference.sh` with the matching frozen Git version, rather than passing old maps to the schema-8 reducer. See [reference replay](../../docs/reproduction.md).

`ReplotPi0BdtScores.C` / `ReplotEtaBdtScores.C` redraw stored BDT ROOT distributions where present. [`CheckEtaDecayFlight.C`](../diagnostics/eta_decay_flight/README.md) is an independent DST diagnostic and does not change composition/topology selections; arguments and interpretation are in the definition notes.

## Implementation and validation

The executable reduce / merge macros delegate to `analysis/photon_candidate_selection/`. Shared histogram, metadata and plotting helpers live there;
shared settings and plotting style live in `analysis/common/`. See [internal implementation](../../analysis/README.md).
Production tree/map validators and small-sample QA live under [validation](../../validation/README.md). `run_map.sh` still validates each temporary output before publishing it.

## Checks

```bash
bash PhotonAnalysisTree/tests/run_candidate_selection_tests.sh "$PWD/PhotonAnalysisTree/build/candidate"
```

The checks cover shared topology/eta/decay-origin behavior, branch roundtrips, metadata, closure, category covariance, veto survival and a synthetic schema-8 map → reduce → merge pipeline. Logs use one fixed `output/condor/photon_candidate_selection/<stage>_<ProcId>.out/.err` per job position and are reused across productions. HTCondor event logs use fixed `/tmp/` paths.
