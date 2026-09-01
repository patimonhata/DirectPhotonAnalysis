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
PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/<sample>/map_<chunk>.root
~~~

The current output schema is version 3. The four isolation branches `split_cluster_iso_raw_et`, `split_cluster_iso_corrected_et`, `split_cluster_iso_boundary`, and `split_cluster_noniso_boundary` are stored as `std::vector<double>`. Schema-2 files stored these branches as floats and must not be mixed with schema-3 files in one reduce input.

The default Condor configuration uses 10 DST segments per ROOT file. This is deliberately configurable through `files_per_job`; after measuring the first jobs, change both `files_per_job` and `n_chunks = ceil(total_files / files_per_job)` together if a different file size is preferable.

For a local one-chunk test:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_map.sh \
  0 0 10001 10 \
  PhotonAnalysisTree/input/jet5/segments.list \
  jet5 \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/jet5 \
  10
~~~

The final `10` limits the test to ten events. Omit it, or use `0`, for all events in the selected DST range. `run_map.sh` writes to a temporary file, runs the ROOT validator, and only then atomically publishes `map_<chunk>.root`. It refuses to overwrite an existing output.

The older `run.sh` and `Fun4All_PythiaPhotonCandidateTree.C` remain useful for a single-segment debug run; production should use the manifest-based map interface.

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
- Only split clusters from `CLUSTERINFO_CEMC` with `E > 0.1 GeV` are stored.
- Shower shapes use the existing 7x7 calculator and a 70 MeV per-tower threshold. Tower patches, constituent towers, and the all-pairs table are not populated.
- Raw isolation is `sum(TopoCluster ET, deltaR < 0.4) - candidate ET`.
- Corrected simulation isolation is `1.2 * raw + 0.1 GeV`.
- Isolated means `corrected_iso < 0.490 + 0.037 * candidate_ET`.
- Non-isolated means `corrected_iso > 0.490 + 0.037 * candidate_ET + 0.8 GeV`.
- The common ABCD candidate requirement is `5 < ET < 35 GeV`, `|eta| < 0.7`, and the shower-shape preselection.
- Tight means `score > 0.8156 - 0.00156 ET`.
- Non-tight is strictly `0.7333 - 0.01333 ET < score < 0.6844 + 0.00156 ET`.

Meson partners are every other split cluster with `E > 0.5 GeV`, without a partner eta cut. The stored best partner is the tagged pair closest to the nominal meson mass.

## Detailed pi0 topology

The map runs `Pi0AnchorTopologyEvaluator` with its Pythia defaults, `min_cluster_energy = 0.1 GeV` with a strict threshold, `|eta_anchor| < 0.7`, `|z_vertex| < 60 cm`, and missing-partner diagnostics enabled. It saves per-cluster truth contributors and prompt/pi0-anchor classifications, all pi0 candidate and daughter recovery diagnostics, and the flattened anchor table. The future reduce stage can reproduce the detailed `pi0_anchor_topology` classification without reopening the DST.

## Condor production

There is one submit file per requested jet sample, so each `condor_submit` invocation creates a separate Condor cluster:

- `submit_jet5.job`: 10,001 DST segments, 1,001 jobs
- `submit_jet8.job`: 10,000 DST segments, 1,000 jobs
- `submit_jet12.job`: 100,000 DST segments, 10,000 jobs
- `submit_jet20.job`: 10,000 DST segments, 1,000 jobs
- `submit_jet30.job`: 10,000 DST segments, 1,000 jobs
- `submit_jet40.job`: 10,000 DST segments, 1,000 jobs

`submit_jet8.job` is configured for the full schema-3 regeneration and writes to `output/intermediate_files/photon_candidate_selection/jet8/schema3_iso_double`. This preserves the existing schema-2 Jet8 maps and prevents accidental mixed-schema reduction.

Create the shared log directory once, review the paths and counts, then submit manually:

~~~bash
mkdir -p PhotonAnalysisTree/output/condor/photon_candidate_selection

condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet5.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet8.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet12.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet20.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet30.job
condor_submit PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet40.job
~~~

No repository script submits jobs automatically.

Each Condor job validates its output before publication. After a sample finishes, check that every expected map file exists:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/check_sample_map_outputs.sh jet5
~~~

Use a third argument of `true` to rerun the full ROOT validator on every file:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/check_sample_map_outputs.sh jet5 10 true
~~~

For the regenerated Jet8 output, pass its explicit output directory as the fourth argument:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/check_sample_map_outputs.sh \
  jet8 10 true \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/jet8/schema3_iso_double
~~~

## Stitching and normalization

PhotonJet windows use the leading terminal prompt HepMC photon pT (classifier category 1 or 2). Jet windows use the leading `AntiKt_Truth_r04` jet pT. Events are retained regardless of the window and carry `sample_stitching_valid` and `sample_stitching_pass`; downstream aggregation must require both.

For each complete sample, the reduce stage must calculate

~~~text
sumw_sample = sum(metadata.sum_generator_weight_processed)
event_weight_pb = event.weight_numerator_pb / sumw_sample
~~~

where `weight_numerator_pb = sample_cross_section_pb * generator_weight`. The denominator includes events rejected later by the vertex cut. Do not normalize map files independently. Cross sections and half-open stitching windows are stored in metadata; `jet40` has no upper bound.

Plotting, cross-map normalization finalization, ABCD purity extraction, and final Region-A topology aggregation are intentionally separate reduce/plot stages and are not part of this map workflow.
