# Pythia photon-candidate selection

This workflow is the map stage of the Pythia photon-candidate analysis. A map job reads a contiguous range of synchronized DST segments and writes one event-wise ROOT TTree containing split-cluster shower shapes, BDT and isolation selections, neutral-meson tags, and detailed pi0 truth topology.

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

`run_small_sample.sh` runs the production map code over the first contiguous part of one sample manifest and then runs the matching Jet- or PhotonJet-family reduce in partial-production mode. For example, this command processes the first 30 Jet12 segments in three maps of 10 segments each and makes the final plots:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_small_sample.sh jet12 30 10 0.2
~~~

Its interface is:

~~~text
run_small_sample.sh SAMPLE_NAME N_SEGMENTS [FILES_PER_MAP] [MIN_CLUSTER_ENERGY_GEV] [OUTPUT_ROOT] [N_EVENTS_PER_MAP] [TAGGING_PARTNER_MIN_ENERGY_GEV]
~~~

The default output root is `PhotonAnalysisTree/output/qa/photon_candidate_selection/cluster_e_gt_<threshold>/<sample>_<N>segments`. It contains `maps/<sample>/map_*.root` and the reduce products under `plots/`. Existing map files are never overwritten. `MIN_CLUSTER_ENERGY_GEV` defaults to 0.1, `TAGGING_PARTNER_MIN_ENERGY_GEV` defaults to the topology threshold, and `N_EVENTS_PER_MAP` defaults to zero, meaning every event in each selected map range.

This mode requires the selected manifest range to start at row zero and remain contiguous, matching the reducer's map-completeness checks. It uses only the selected sample and normalizes with only the maps present, so its products are for code, schema, and plot QA only—not a physics result. A handful of segments can be enough for a smoke test, but 30--100 segments is more likely to populate the Region-A topology plots.

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

## Region-A pi0-anchor topology reduce

`ReducePythiaRegionAAnchorTopology.C` reads the completed map production, validates every sample's schema, metadata compatibility, and contiguous manifest coverage, calculates the full-sample normalization above, and makes the Region-A topology plots without reopening DSTs. Jet and PhotonJet are separate families and are never mixed.

The event loop disables every branch by default, activates only the 11 scalar/vector branches used by this reduce, and gives those branches a 64 MiB `TTreeCache` per open map file. Keep future event-branch bindings on the `bind_active` path so new inputs are both enabled and cached.

Run the compiled Jet-family reduce locally with:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_reduce.sh \
  jet \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p1 \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/region_a_pi0_anchor_topology/cluster_e_gt_0p1/jet \
  true
~~~

Run the PhotonJet family independently by replacing `jet` with `photonjet` in both the family argument and output path:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_reduce.sh \
  photonjet \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p1 \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/region_a_pi0_anchor_topology/cluster_e_gt_0p1/photonjet \
  true
~~~

Production mode requires every expected map from every sample in the selected family. The optional fourth argument `require_complete=false` is only for development with partial productions; its cross-section normalization uses only the maps present and is not a physics result.

For a disconnect-safe batch run, submit the Jet and PhotonJet reducers as separate one-job Condor clusters:

~~~bash
mkdir -p PhotonAnalysisTree/output/condor/photon_candidate_selection/reduce
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_jet.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_photonjet.job
~~~

The submit files default to `configuration = cluster_e_gt_0p1` and `require_complete = true`. Edit both values as needed before submission. Each output base must have at most one active reducer; the repository never submits these jobs automatically.

The cluster selections are:

- pi0 anchor: `sample_stitching_valid && sample_stitching_pass && split_cluster_pass_region_a && split_cluster_pi0_anchor_valid`;
- prompt reference: `sample_stitching_valid && sample_stitching_pass && split_cluster_pass_region_a && split_cluster_truth_prompt_cluster`.

The stored topology classification uses the map-time strict `E_cluster > min_cluster_energy` threshold. The reduce validates that all maps and samples use the same topology threshold, fixed diagnostic floor, tagging threshold, and topology algorithm version, and does not reclassify topology. Its default binning matches the current `pi0_anchor_topology` production: 200 bins over `0 <= ET < 40 GeV`.

Each family produces one selection-comparison ROOT file and one directory per selection under

~~~text
output/plots/photon_candidate_selection/region_a_pi0_anchor_topology/cluster_e_gt_<threshold>/<family>/
├── selection_comparison.root
├── kinematic/region_a_pi0_anchor_topology*.pdf
├── preselection/region_a_pi0_anchor_topology*.pdf
├── preselection_tight/region_a_pi0_anchor_topology*.pdf
├── preselection_isolation/region_a_pi0_anchor_topology*.pdf
└── region_a_tagging_veto/region_a_pi0_anchor_topology*.pdf
~~~

The five selections are `kinematic`, `preselection` (`kinematic && pass_preselection`), `preselection_tight`,
`preselection_isolation`, and `region_a_tagging_veto` (`region_a && !pi0_tag && !eta_tag`). Each selection directory
contains the six spectrum/fraction PDFs listed below, with the selection label included in the annotations.

The ROOT file stores one directory per selection, each containing unweighted counts, cross-section-weighted spectra in pb,
bin-width-normalized spectra in pb/GeV, and category fractions. The `selection_summary` tree stores per-sample raw selected
cluster and selected anchor counts. Its metadata records the exact selection definitions and compact reduce provenance. The
topology categories are checked to partition the selected Region-A anchor denominator in every ET bin, including underflow
and overflow.

## Photon-candidate purity and background composition

`ReducePythiaPhotonCandidateComposition.C` partitions every selected candidate into exactly one of prompt photon, pi0 topology, eta, or other. Every truth-origin requirement uses a strict contribution `> 0.5`; an exact contribution of `0.5` is therefore other. The detailed pi0 categories are separated, merged, single contaminated, missing, and other.

Prompt uses the stored prompt flag plus the stored dominant-contributor fraction. Pi0 uses the stored anchor plus its same-parent main fraction and topology. Eta sums signal-embedding contributors identified either as a G4 eta or as a generator eta-decay photon. If more than one strict majority is ever found, the candidate is assigned to other and the reducer returns nonzero after writing the overlap count to metadata.

Run the Region-A Jet-family composition with:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_reduce_composition.sh \
  jet \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p5 \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/candidate_composition/cluster_e_gt_0p5/region_a/jet/photon_candidate_composition \
  region_a \
  true
~~~

Use `final_photon` instead of `region_a` to use the stored Region-A selection after the pi0-or-eta tag veto:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_reduce_composition.sh \
  jet \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p5 \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/candidate_composition/cluster_e_gt_0p5/final_photon/jet/photon_candidate_composition \
  final_photon \
  true
~~~

For a disconnect-safe full Jet-family production, create the shared log directory and submit one independent reduce job per sample:

~~~bash
mkdir -p PhotonAnalysisTree/output/condor/photon_candidate_selection/reduce
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_composition_jet_samples.job
~~~

To retry only a subset without rerunning completed samples, override `sample_names` on submission. For example:

~~~bash
condor_submit -append 'sample_names = jet3,jet5' \
  PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_composition_jet_samples.job
~~~

This creates seven ordinary Condor jobs, one each for `jet3`, `jet5`, `jet8`, `jet12`, `jet20`, `jet30`, and `jet40`. The concurrency limit is four. Each job writes only

~~~text
partial/<sample>/photon_candidate_composition.root
~~~

The partial jobs do not create PDFs. After all seven jobs have exited successfully, merge them manually:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_merge_composition.sh \
  jet \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/candidate_composition/cluster_e_gt_0p5/region_a/jet/partial \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/candidate_composition/cluster_e_gt_0p5/region_a/jet/photon_candidate_composition \
  region_a
~~~

There is deliberately no automatic Condor dependency or DAG. The merge requires exactly one partial for every sample in the family and rejects missing files, unexpected sample metadata, incompatible analysis/configuration metadata, incompatible histogram axes, and invalid category partitions. It adds only count and weighted spectra, then recomputes every fraction and both stack PDFs from the merged weighted numerator and denominator.

The sample submit defaults to `configuration = cluster_e_gt_0p5`, `selection = region_a`, and `require_complete = true`. Change `selection` to `final_photon` for the tag-veto result. Submit and merge the two selections separately because they use different partial and final directories and different log names.

Omitting the optional final `SAMPLE_NAME` argument from `run_reduce_composition.sh` retains the original all-sample behavior. The original one-job submit file also remains available:

~~~bash
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_reduce_composition_jet.job
~~~

Do not run the one-job and sample-split reducers with the same output paths.

Each run writes:

- `photon_candidate_composition.root`, containing unweighted and cross-section-weighted spectra, `h_photon_candidate_purity`, every category fraction, counts, normalization inputs, and classification QA metadata;
- `photon_candidate_composition_category_fraction_stack.pdf`, whose summary categories partition the selected-candidate denominator in every cluster-ET bin;
- `photon_candidate_composition_category_fraction_stack_detailed.pdf`, with the individual pi0 topology categories.

The reducer validates both the unweighted and weighted partition, including underflow and overflow. `overlap_cluster_count`, `half_boundary_cluster_count`, and `invalid_truth_cluster_count` are retained as QA counters. Production results require `REQUIRE_COMPLETE=true`; partial mode is only for development.

ABCD purity extraction remains a separate future reduce stage.
