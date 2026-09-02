# Pythia photon-candidate selection

This workflow covers the map, reduce, and merge stages of the Pythia photon-candidate analysis. A map job reads a contiguous range of synchronized DST segments and writes one event-wise ROOT TTree containing split-cluster shower shapes, BDT and isolation selections, neutral-meson tags, and detailed pi0 truth topology.

## Build

Use ana.565.

~~~bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
cmake -S PhotonAnalysisTree/src -B PhotonAnalysisTree/build/candidate \
  -DCMAKE_INSTALL_PREFIX="$PWD/PhotonAnalysisTree/install" \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build PhotonAnalysisTree/build/candidate --parallel 4
cmake --install PhotonAnalysisTree/build/candidate
~~~

Do not rebuild or reinstall while Condor jobs using the installed library are running.

## Map input and output

The input unit is a range in `PhotonAnalysisTree/input/<sample>/segments.list`. For every suffix in the range, the Fun4All macro adds these four synchronized logical file names, in the same order:

- `DST_CALO_CLUSTER_<suffix>`
- `DST_MBD_EPD_<suffix>`
- `DST_TRUTH_JET_<suffix>`
- `G4Hits_<suffix>`

The files do not need to be copied into the repository or current directory. `Fun4AllDstInputManager::AddFile` resolves correctly named production files through the normal sPHENIX file-catalog mechanism.

A production map job processes all events in several DST segments and writes one file:

~~~text
PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_<threshold>/<sample>/map_<chunk>.root
~~~

The current output schema is version 4. In addition to the schema-3 double-precision isolation branches, schema 4 stores independent pi0-partner alignment, truth-pair taggability, and observed veto-result axes. Schema versions must not be mixed in one reduce input.

The default Condor configuration uses 10 DST segments per ROOT file. This is deliberately configurable through `files_per_job`; after measuring the first jobs, change both `files_per_job` and `n_chunks = ceil(total_files / files_per_job)` together if a different file size is preferable.

For a local one-chunk test:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_map.sh \
  0 0 10000 10 \
  PhotonAnalysisTree/input/jet5/segments.list \
  jet5 \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p2/jet5 \
  10 \
  0.2
~~~

The penultimate `10` limits the test to ten events, and the final `0.2` sets the strict topology threshold to `E_cluster > 0.2 GeV`. An optional tenth argument sets the strict tagging-partner threshold; when omitted it equals the topology threshold. `N_EVENTS` defaults to zero and the topology threshold defaults to 0.1 GeV. `run_map.sh` writes to a temporary file, verifies both thresholds and the fixed 0.1 GeV diagnostic floor in metadata, runs the ROOT validator, and only then atomically publishes `map_<chunk>.root`. It refuses to overwrite an existing output.

All map production, including single-chunk tests, uses this manifest-based interface.

## Small-sample end-to-end QA

`run_small_sample.sh` runs the production map code over the first contiguous part of one sample manifest and then runs the matching Jet- or PhotonJet-family reduce in partial-production mode. For example, this command processes the first 100 Jet12 segments in ten maps of 10 segments each and writes all ten integrated partials:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_small_sample.sh jet12 100 10 0.2
~~~

Its interface is:

~~~text
run_small_sample.sh SAMPLE_NAME N_SEGMENTS [FILES_PER_MAP] [MIN_CLUSTER_ENERGY_GEV] [OUTPUT_ROOT] [N_EVENTS_PER_MAP] [TAGGING_PARTNER_MIN_ENERGY_GEV]
~~~

The default output root is `PhotonAnalysisTree/output/qa/photon_candidate_selection/cluster_e_gt_<threshold>/<sample>_<N>segments`. It contains `maps/<sample>/map_*.root` and integrated partials under `reduce/<sample>/shard_<index>/`. Existing map files are never overwritten. `MIN_CLUSTER_ENERGY_GEV` defaults to 0.1, `TAGGING_PARTNER_MIN_ENERGY_GEV` defaults to the topology threshold, and `N_EVENTS_PER_MAP` defaults to zero, meaning every event in each selected map range.

This mode requires the selected manifest range to start at row zero and remain contiguous, matching the reducer's map-completeness checks. It uses only the selected sample and normalizes with only the maps present, so its products are for code, schema, and plot QA only—not a physics result. A non-Jet12 sample can use a handful of segments for a smoke test. Jet12 requires at least ten generated map files so every fixed shard is non-empty.

## Region A/B/C/D content

Every stored split cluster has these independent flag branches:

- `split_cluster_pass_region_a`: isolated and tight
- `split_cluster_pass_region_b`: non-isolated and tight
- `split_cluster_pass_region_c`: isolated and non-tight
- `split_cluster_pass_region_d`: non-isolated and non-tight
- `split_cluster_pass_final_photon`: Region A after the pi0-or-eta tag veto

The event tree also stores `region_a_count`, `region_b_count`, `region_c_count`, `region_d_count`, and `final_photon_count`. The one-entry metadata tree stores their map-file totals as `n_clusters_region_a/b/c/d` and `n_clusters_final_photon`. Clusters in the isolation gap or BDT gap have no A/B/C/D flag, which is intentional.

Thus the future ABCD purity calculation can be performed entirely from these intermediate TTrees. No Region-B/C/D cluster is discarded at the map stage.

## Selection definitions

- The signal HepMC collision vertex must satisfy `|z| < 60 cm`; missing or rejected vertices are not written.
- `TOPOCLUSTER_ALLCALO` is reconstructed in the job with the requested EMCal+HCal topological-clustering configuration.
- Only split clusters from `CLUSTERINFO_CEMC` with `E > min_cluster_energy` are stored; the threshold is strict and configured per map production.
- Shower shapes use the existing 7x7 calculator and a 70 MeV per-tower threshold. Tower patches, constituent towers, and the all-pairs table are not populated.
- Raw isolation is `sum(TopoCluster ET, deltaR < 0.4) - candidate ET`.
- Corrected simulation isolation is `1.2 * raw + 0.1 GeV`.
- Isolated means `corrected_iso < 0.490 + 0.037 * candidate_ET`.
- Non-isolated means `corrected_iso > 0.490 + 0.037 * candidate_ET + 0.8 GeV`.
- The common ABCD candidate requirement is `5 < ET < 35 GeV`, `|eta| < 0.7`, and the shower-shape preselection.
- Tight means `score > 0.8156 - 0.00156 ET`.
- Non-tight is strictly `0.7333 - 0.01333 ET < score < 0.6844 + 0.00156 ET`.

Meson tagging considers every other stored split cluster with `E > tagging_partner_min_energy`, without a partner eta cut. This tagging threshold is configured per production and defaults to the topology threshold in the wrappers. The stored best partner is the tagged pair closest to the nominal meson mass. The definition of `split_cluster_pass_final_photon` remains `region_a && !pi0_tag && !eta_tag`.

## Detailed pi0 topology

The map runs `Pi0AnchorTopologyEvaluator` with its Pythia defaults, the same strict `min_cluster_energy` used for stored split clusters, `|eta_anchor| < 0.7`, `|z_vertex| < 60 cm`, and missing-partner diagnostics enabled. It saves per-cluster truth contributors and prompt/pi0-anchor classifications, all pi0 candidate and daughter recovery diagnostics, and the flattened anchor table. The reduce stage can reproduce the detailed `pi0_anchor_topology` classification without reopening the DST.

The topology axis remains `separated/merged/single_contaminated/missing/other` and is intentionally evaluated with the production topology threshold. Independently, the truth partner is the cluster above the fixed strict diagnostic floor `E > 0.1 GeV` with the largest absolute direct energy deposit from the partner photon. `pi0_anchor_partner_alignment` records `near` or `displaced` using the existing projection criterion `deltaR > 0.15`, with explicit invalid/unavailable states.

`pi0_anchor_truth_partner_tag_status` records whether that representative truth pair is taggable at the production tagging threshold and pi0 mass window. `pi0_anchor_tag_result` then separates a surviving anchor, a veto where the truth pair itself was taggable, and a combinatorial-only veto where some reconstructed pair tagged the anchor although the representative truth pair was not taggable. `pi0_anchor_selected_tag_partner_matches_truth_partner` additionally records whether the selected pi0-tag partner cluster is the representative truth partner. These axes are diagnostic only and do not change topology or `final_photon`.

The integer encodings are:

- `partner_alignment`: 0 not applicable, 1 near, 2 displaced, 3 projection invalid, 4 cluster unavailable;
- `truth_partner_tag_status`: 0 not applicable, 1 taggable, 2 cluster unavailable, 3 same as anchor, 4 below tagging threshold, 5 mass outside window, 6 invalid mass;
- `tag_result`: 0 not applicable, 1 survived, 2 truth pair taggable veto, 3 combinatorial-only veto.

## Condor production

There is one submit file per sample, so each `condor_submit` invocation creates a separate Condor cluster:

- `submit_jet3.job`: 10,000 DST segments, 200 jobs
- `submit_jet5.job`: 10,000 DST segments, 200 jobs
- `submit_jet8.job`: 10,000 DST segments, 1,000 jobs
- `submit_jet12.job`: 100,000 DST segments, 10,000 jobs
- `submit_jet20.job`: 10,000 DST segments, 1,000 jobs
- `submit_jet30.job`: 10,000 DST segments, 1,000 jobs
- `submit_jet40.job`: 10,000 DST segments, 1,000 jobs
- `submit_photonjet3.job`: 10,000 DST segments, 1,000 jobs
- `submit_photonjet5.job`: 10,000 DST segments, 1,000 jobs
- `submit_photonjet10.job`: 10,000 DST segments, 1,000 jobs
- `submit_photonjet20.job`: 10,000 DST segments, 1,000 jobs

Create the shared log directory once, review the paths and counts, then submit manually:

~~~bash
mkdir -p PhotonAnalysisTree/output/condor/photon_candidate_selection

condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet3.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet5.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet8.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet12.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet20.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet30.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet40.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_photonjet3.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_photonjet5.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_photonjet10.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_photonjet20.job
~~~

No repository script submits jobs automatically.
Each submit file defines `min_cluster_energy` and `tagging_partner_min_energy`; the latter defaults to the former but can be changed independently per production. The wrapper likewise defaults an omitted tenth argument to the topology threshold. The fixed partner diagnostic floor is always 0.1 GeV. Change the threshold component of `map_output_directory` consistently and never mix topology/tagging configurations in one sample directory.

Each Condor job validates its output before publication. After a sample finishes, check that every expected map file exists:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/check_sample_map_outputs.sh jet5
~~~

Use a third argument of `true` to rerun the full ROOT validator on every file:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/check_sample_map_outputs.sh jet5 10 true
~~~

For non-default thresholds, pass the output directory and expected topology threshold; an optional sixth argument supplies a distinct expected tagging threshold:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/check_sample_map_outputs.sh \
  jet5 10 true \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p2/jet5 \
  0.2
~~~

## Stitching and normalization

PhotonJet windows use the leading terminal prompt HepMC photon pT (classifier category 1 or 2). Jet windows use the leading `AntiKt_Truth_r04` jet pT. Events are retained regardless of the window and carry `sample_stitching_valid` and `sample_stitching_pass`; downstream aggregation must require both.

For each complete sample, the reduce stage must calculate

~~~text
sumw_sample = sum(metadata.sum_generator_weight_processed)
event_weight_pb = event.weight_numerator_pb / sumw_sample
~~~

where `weight_numerator_pb = sample_cross_section_pb * generator_weight`. The denominator includes events rejected later by the vertex cut. Do not normalize map files independently. Cross sections and half-open stitching windows are stored in metadata; `jet40` has no upper bound.

## Unified reduce and merge

`ReducePythiaPhotonCandidateSelection.C` is the only production reducer. One job reads one sample shard once and fills both analyses:

- photon-candidate composition for the requested selection;
- pi0-anchor topology for all five comparison selections.

The supported composition selections are `kinematic`, `preselection`, `preselection_tight`, `preselection_isolation`, `region_a`, and `region_a_tagging_veto`. `final_photon` is an alias of `region_a_tagging_veto`. The topology comparison always includes `kinematic`, `preselection`, `preselection_tight`, `preselection_isolation`, and `region_a_tagging_veto`, independently of the requested composition selection.

The composition categories partition every selected candidate into prompt photon, pi0 topology, eta, or other. The detailed pi0 categories are separated, merged, single contaminated, missing, and other. Truth-origin majorities use a strict contribution `> 0.5`; an exact contribution of `0.5` is other. An overlap is recorded in metadata and makes the reducer return code 9 after writing its output.

The topology prompt reference requires both the selected flag and `split_cluster_truth_prompt_cluster`. The pi0 categories partition selected valid pi0 anchors. Weighted fraction errors are computed from the weighted numerator/denominator subset covariance; weighted fractions are never added directly.

### Shards and partial files

Every non-Jet12 sample has exactly one shard, index 0. Jet12 has exactly ten contiguous, non-overlapping shards, indices 0 through 9. A Jet12 shard reads only its map range but uses the full Jet12 generator-weight sum for cross-section normalization.

Run one partial locally with the mandatory interface:

~~~text
run_reduce.sh FAMILY MAP_ROOT OUTPUT_BASE SELECTION REQUIRE_COMPLETE N_BINS ET_MAX_GEV SAMPLE_NAME SHARD_INDEX
~~~

For example:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_reduce.sh \
  jet \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p5 \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/reduce/cluster_e_gt_0p5/region_a/jet/partial/jet5/shard_0/photon_candidate_selection \
  region_a true 200 40.0 jet5 0
~~~

Each partial is:

~~~text
partial/<sample>/shard_<index>/photon_candidate_selection.root
├── composition/
├── anchor_topology/<selection>/
├── metadata
├── topology_summary
└── shard_metadata
~~~

The partial contains only count and cross-section-weighted histograms plus validation metadata; final fractions and PDFs are created by merge.

### Condor reduce

For the Jet family, submit six one-shard samples and ten Jet12 shards as 16 ordinary Condor jobs:

~~~bash
mkdir -p PhotonAnalysisTree/output/condor/photon_candidate_selection/reduce
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_jet_samples.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_jet12_shards.job
~~~

The submit files default to `configuration = cluster_e_gt_0p5`, `selection = region_a`, `require_complete = true`, 200 bins, and 40 GeV. Edit the configuration and selection consistently in both submit files before submitting. They write distinct partial and log paths for every sample/shard.

Retry only selected jobs with:

~~~bash
condor_submit -append 'sample_names = jet3,jet5' \
  PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_jet_samples.job
condor_submit -append 'shard_indices = 3,7' \
  PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_jet12_shards.job
~~~

There is no DAG or automatic dependency. Confirm that all 16 jobs exited successfully before merge.

### Merge

`MergePythiaPhotonCandidateSelection.C` is the only merger. Its wrapper interface is:

~~~text
run_merge.sh FAMILY PARTIAL_ROOT COMPOSITION_OUTPUT_BASE TOPOLOGY_OUTPUT_BASE SELECTION
~~~

For example:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_merge.sh \
  jet \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/reduce/cluster_e_gt_0p5/region_a/jet/partial \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/candidate_composition/cluster_e_gt_0p5/region_a/jet/photon_candidate_composition \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/region_a_pi0_anchor_topology/cluster_e_gt_0p5/jet \
  region_a
~~~

The merger requires shard 0 for each non-Jet12 sample and shards 0--9 for Jet12. It rejects missing or duplicate coverage, invalid shard ranges, inconsistent full-sample normalization, unexpected sample metadata, incompatible analysis/configuration metadata, incompatible axes, and invalid category partitions. It adds only counts and weighted spectra, then recomputes all fractions and errors.

Composition output:

~~~text
photon_candidate_composition.root
photon_candidate_composition_category_fraction_stack.pdf
photon_candidate_composition_category_fraction_stack_detailed.pdf
~~~

The ROOT file includes unweighted and weighted spectra, photon purity, every category fraction, normalization inputs, and classification QA counters.

Anchor-topology output:

~~~text
<TOPOLOGY_OUTPUT_BASE>/
├── selection_comparison.root
├── kinematic/region_a_pi0_anchor_topology*.pdf
├── preselection/region_a_pi0_anchor_topology*.pdf
├── preselection_tight/region_a_pi0_anchor_topology*.pdf
├── preselection_isolation/region_a_pi0_anchor_topology*.pdf
└── region_a_tagging_veto/region_a_pi0_anchor_topology*.pdf
~~~

The topology ROOT file stores counts, weighted spectra in pb, bin-width-normalized spectra in pb/GeV, fractions, metadata, and per-selection sample summaries.

Production results require `REQUIRE_COMPLETE=true`. `false` is for incomplete QA productions only; normalization then uses only available maps and is not a physics result. Jet and PhotonJet families remain separate and are never mixed.

ABCD purity extraction remains a separate future reduce stage.
