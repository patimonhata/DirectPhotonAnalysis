# Pythia photon-candidate selection

The development map writes schema **11**, common composition **1**, pi0 topology **14**, eta topology **4**, and supplemental decay-origin **1**. The existing reducer reads schema-10 maps and writes schema-13 partials. Use matching frozen software for retained productions. Physics definitions and diagnostic details are in [definition notes](../../docs/photon_candidate_definitions.md).

## Results and inputs

```text
output/photon_candidate_selection/
  productions/
    <production>/<configuration>/
      maps/<sample>/map_<chunk>.root
      partials/jet/<sample>/shard_<index>/photon_candidate_selection.root
      results/jet/
    topology_v11/<configuration>/  # retained schema-5 reference
  qa/
```

Some retained map directories link to cold storage. Do not remove the underlying map inputs when reorganizing diagnostics. Saved ROOT metadata retains its original paths; `output/layout_manifest.json` records where files moved. See [layout](../../docs/layout.md).

Each map reads synchronized `DST_CALO_CLUSTER`, `DST_MBD_EPD`, `DST_TRUTH_JET`, and `G4Hits` basenames from `input/<sample>/segments.list`. Fun4All resolves the logical filenames through the file catalog. Jet submit defaults use 50 segments per map, 200 maps per ordinary sample and 2,000 for Jet12. Exact selection settings are in each submit file and map metadata. Jet defaults are cluster E>0.12 GeV, pi0/eta partner E>0.5 GeV, pi0 0.10<m<0.20 GeV, eta 0.45<m<0.65 GeV, missing band 0.12<E<=0.5 GeV, and recovery disabled. PhotonJet defaults are a separate configuration; do not mix them with Jet outputs.

## Build

See [software/build provenance](../../docs/software.md) for the managed build,
versioned build and install paths, and selection through `install/current`.
For formal productions, commit the execution code, rebuild clean source, freeze a software ID, and initialize a new production manifest.

For maps on tg with manifests, partials, results and logs on `/sphenix/user`, initialize a new production with `software.py production ... --map-root /absolute/tg/path`.
The storage directory must be new or empty; the command records its canonical path and creates `configuration/maps` as a symlink.
This requires newly frozen software containing the external-map-root feature. Existing submit arguments remain valid, and `log_directory` stays independent.
HTCondor event logs (`Log`) must stay under `/tmp/` by site policy; retain the submit defaults. `log_directory` controls stdout/stderr only.
See [external map storage](../../docs/software.md#map-だけを-tg-に保存する) for the full setup.

## Map

Review the production/configuration and choose new map outputs when rerunning a production. Existing map files are not overwritten. No repository script automatically submits jobs.

The default model is `models/model_ppg15_nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_single_tmva.root`, loaded as `myBDT` by `TMVA::Experimental::RBDT`. `run_map.sh` passes the model path explicitly from `PHOTON_SOURCE_ROOT`; the macro default uses the same model. The existing tight/non-tight thresholds are retained intentionally. See [model and input order](../../docs/photon_candidate_definitions.md#bdt-model-and-input-order). Model changes require new DST maps; reduce, merge and replot do not recompute scores or selection flags.

Submit files default to `production_directory = UNCONFIGURED`. Initialize a new production manifest, then pass its configuration directory to every map/reduce submission. The submit file sets `PHOTON_PRODUCTION_MANIFEST`; runtime dispatch selects the bound frozen scripts and libraries. `software.py production` creates the manifest without exporting shell variables, and normal production submissions do not need a separate `PHOTON_SOFTWARE_ROOT`.

```bash
mkdir -p PhotonAnalysisTree/output/condor/photon_candidate_selection
condor_submit -append "production_directory = /absolute/path/to/new_configuration" PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet5.job
```

Other samples have matching submit files. The four final arguments default to 3 GeV, 60 cm, 0.7, and 0.7. Set the MBD z limit to 100 cm and the candidate eta limit to 1.0 for wide reconstructed acceptance. The topology anchor limit uses truth-vertex eta and can be configured independently.

The explicit local interface is:

```text
run_map.sh JOB_INDEX CHUNK_OFFSET TOTAL_FILES FILES_PER_JOB INPUT_MANIFEST SAMPLE OUTPUT_DIRECTORY [N_EVENTS] [MIN_CLUSTER_ENERGY] [PI0_PARTNER_MIN_ENERGY] [ETA_PARTNER_MIN_ENERGY] [PI0_MASS_MIN] [PI0_MASS_MAX] [ETA_MASS_MIN] [ETA_MASS_MAX] [MISSING_ENERGY_MIN] [MISSING_ENERGY_MAX] [SHOWER_SHAPE_MIN_CLUSTER_ET] [MAX_ABS_MBD_VERTEX_Z] [CANDIDATE_ABS_ETA_MAX] [TOPOLOGY_ANCHOR_ABS_ETA_MAX]
```

[`check_sample_map_outputs.sh`](../../validation/photon_candidate_selection/check_sample_map_outputs.sh) checks expected file coverage and optionally performs deep ROOT validation; pass the output directory and exact cuts for a non-default production. [`run_small_sample.sh`](../../validation/photon_candidate_selection/run_small_sample.sh) runs the schema-10 map and reduce chain with matching frozen software on a small input sample under `output/photon_candidate_selection/qa/`; this is not a physics normalization. Jet12 QA needs at least ten maps so every fixed shard is nonempty.

## Reduce and merge

One reducer fills all six selections: `kinematic`, `preselection`, `preselection_tight`, `preselection_isolation`, `region_a`, and `region_a_tagging_veto`. The last is Region A after the pi0-or-eta veto. The input maps also retain Region B/C/D and gap flags for later ABCD studies.

Use full-sample generator-weight normalization, including processed events later rejected by the vertex cut. Never normalize individual map files or add per-shard normalization denominators. Jet12 has ten contiguous shards; each other Jet sample has one. Schema/settings, exact sample/shard coverage, category closure, weighted Sumw2 and subset relations are checked before final plots are produced.

```bash
condor_submit -append "production_directory = /absolute/path/to/new_configuration" PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_jet_samples.job
condor_submit -append "production_directory = /absolute/path/to/new_configuration" PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_jet12_shards.job
```

Interfaces:

```text
run_reduce.sh FAMILY MAP_ROOT OUTPUT_BASE REQUIRE_COMPLETE N_BINS ET_MAX SAMPLE SHARD_INDEX
run_merge.sh FAMILY PARTIAL_ROOT OUTPUT_DIRECTORY
```

After all 16 Jet partials complete, select the same production manifest and write results within its configuration directory:

```bash
python3 PhotonAnalysisTree/tools/software.py run \
  --production-manifest /absolute/path/to/new_configuration/manifest.json \
  merge jet /absolute/path/to/new_configuration/partials/jet /absolute/path/to/new_configuration/results/jet
```

Merged outputs include all-origin common composition, transition CSVs and matching/ancestry quality, candidate composition, pi0/eta anchor topology, eta cascade composition, Other-origin and supplemental decay-origin histograms, spectra/fractions/survival PDFs, ET-binned BDT distributions and summaries. ROOT histogram metadata records source map/partial versions, cuts, weights and uncertainty conventions. Statistical errors use the documented candidate/pair-level Sumw2 and subset covariance; event-level correlations are not estimated.

## Retained v11/schema-5 results

Four configurations under `productions/topology_v11/` retain maps, 16 partials and final results. Partner thresholds are 0.15 or 0.5 GeV and pi0 mass windows are 0.10–0.20 or 0.08–0.18 GeV. Eta topology extension is absent, but eta tagging remains enabled. Use `run_reference.sh` with the matching frozen Git version, rather than passing old maps to the current schema-10 map reducer. See [reference replay](../../docs/reproduction.md).

`ReplotPi0BdtScores.C` / `ReplotEtaBdtScores.C` redraw stored BDT ROOT distributions where present. [`CheckEtaDecayFlight.C`](../diagnostics/eta_decay_flight/README.md) is an independent DST diagnostic and does not change composition/topology selections; arguments and interpretation are in the definition notes.

## Implementation and validation

The executable reduce / merge macros delegate to `analysis/photon_candidate_selection/`. Shared histogram, metadata and plotting helpers live there;
shared settings and plotting style live in `analysis/common/`. See [internal implementation](../../analysis/README.md).
Production tree/map validators and small-sample QA live under [validation](../../validation/README.md). `run_map.sh` still validates each temporary output before publishing it.

## Checks

See [the software test commands](../../docs/software.md#検証) for the build and install paths to pass.

The checks cover shared topology/eta/decay-origin behavior, branch roundtrips, metadata, closure, category covariance, veto survival and a synthetic schema-10 map → reduce → merge pipeline. By default, stdout and stderr use `output/condor/photon_candidate_selection/<stage>_<ProcId>.out/.err` across productions. Set `log_directory` at submit time to use a production-specific directory. HTCondor event logs use fixed `/tmp/` paths.
