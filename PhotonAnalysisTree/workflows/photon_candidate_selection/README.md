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

The current map schema is 5 and the topology algorithm version is 10. This version disables the photon-energy recovery requirement (`min_photon_energy_recovery = 0`). The direct-deposit matching quality and dominance conditions remain unchanged. Older maps and partials must be regenerated; they cannot be mixed with this production. All partner thresholds, mass windows, missing-energy boundaries, and the recovery requirement are stored and checked in map/reduce/merge metadata.

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

The penultimate `10` limits the test to ten events, and the final `0.2` sets the strict stored-cluster/anchor threshold `E_cluster > 0.2 GeV`. Optional trailing arguments are listed below. `N_EVENTS` defaults to zero and the stored-cluster threshold defaults to 0.1 GeV. Outputs are validated and published atomically; existing maps are never overwritten.

All map production, including single-chunk tests, uses this manifest-based interface.

## Small-sample end-to-end QA

`run_small_sample.sh` runs the production map code over the first contiguous part of one sample manifest and then runs the matching Jet- or PhotonJet-family reduce in partial-production mode. For example, this command processes the first 100 Jet12 segments in ten maps of 10 segments each and writes all ten integrated partials:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_small_sample.sh jet12 100 10 0.2
~~~

Its interface is:

~~~text
run_small_sample.sh SAMPLE_NAME N_SEGMENTS [FILES_PER_MAP] [MIN_CLUSTER_ENERGY_GEV] [OUTPUT_ROOT] [N_EVENTS_PER_MAP] [PI0_PARTNER_MIN_ENERGY_GEV] [ETA_PARTNER_MIN_ENERGY_GEV] [PI0_MASS_MIN] [PI0_MASS_MAX] [ETA_MASS_MIN] [ETA_MASS_MAX] [MISSING_ENERGY_MIN] [MISSING_ENERGY_MAX]
~~~

The default output root is `PhotonAnalysisTree/output/qa/photon_candidate_selection/cluster_e_gt_<threshold>/<sample>_<N>segments`. It contains `maps/<sample>/map_*.root` and integrated partials under `reduce/<sample>/shard_<index>/`. Existing map files are never overwritten. `MIN_CLUSTER_ENERGY_GEV` defaults to 0.1, the pi0 threshold defaults to the stored-cluster threshold and the eta threshold defaults to the pi0 threshold, and `N_EVENTS_PER_MAP` defaults to zero, meaning every event in each selected map range.

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

Meson tagging considers every other valid reconstructed split cluster in `CLUSTERINFO_CEMC`, without a partner eta cut. Pi0 and eta each require their own strict `E > threshold` and `mass_min < mass < mass_max`. Any qualifying pair sets the tag; among multiple pairs, the saved veto partner is the pair closest to the nominal meson mass. Truth matching is not used in this veto. `split_cluster_pass_final_photon` remains `region_a && !pi0_tag && !eta_tag`.

## Partner and missing definitions

The anchor selection keeps `min_cluster_energy`. Independently, topology partner lookup uses the strict pi0 tagging threshold. Each lookup selects the cluster with the largest absolute direct daughter deposit within its own eligible pool. A usable deposit-selected cluster in the respective energy pool is sufficient for recovery; no minimum recovered-energy fraction is required. The calibrated photon energy estimate `cluster_energy * daughter_deposit / total_deposit` and its recovery ratio remain stored as diagnostics. Thus changing the pi0 partner threshold can change separated/merged/single-contaminated/missing/other totals. The saved candidate daughter `best_cluster`/`recovered` fields describe the anchor pool; the anchor topology additionally evaluates the independent partner pool.

For missing diagnostics, all valid positive-energy clusters are searched, including those at or below 0.1 GeV. The representative truth partner maximizes direct daughter deposit among usable truth matches, independently of the tagging threshold. Its full reconstructed cluster energy (not the truth-photon energy or the attributed fraction) determines the missing energy category. The reconstruction's own clustering thresholds remain unchanged.

The six exclusive missing categories are:

| Code | Key | Definition |
|---|---|---|
| 1 | `energy_band_taggable` | `L < E <= U`, representative pair inside pi0 mass window |
| 4 | `energy_band_not_taggable` | `L < E <= U`, finite representative pair mass outside pi0 window |
| 5 | `low_energy` | `E <= L` |
| 6 | `unclustered_or_no_cemc_deposit` | No associated cluster found, or no CEMC deposit |
| 2 | `acceptance` | Partner projection outside CEMC acceptance |
| 3 | `other` | Remaining cases, including unresolved matching, invalid mass, same cluster as anchor, and `E > U` |

`L=0.2` and `U=0.5 GeV` by default. The mass-only taggability in this table asks whether lowering the partner energy threshold could tag the pair; it does not impose the production tagging threshold. Invalid projection goes to other; acceptance takes priority over energy classification. Incomplete matching is not labelled unclustered when it prevents identifying a partner. Displacement no longer defines a missing category. Legacy detailed reason codes remain diagnostic information, not additional histogram categories.

The independent diagnostic axes remain:

- `partner_alignment`: 0 not applicable, 1 near, 2 displaced (`deltaR > 0.15`), 3 projection invalid, 4 cluster unavailable;
- `truth_partner_tag_status`: 0 not applicable, 1 taggable, 2 unavailable, 3 same as anchor, 4 below the production pi0 tagging threshold, 5 mass outside pi0 window, 6 invalid mass;
- `tag_result`: 0 not applicable, 1 survived, 2 truth pair taggable veto, 3 combinatorial-only veto.

These axes use the deposit-selected truth partner, which need not be the reconstructed partner selected by the actual veto.

## Configurable parameters

Map submit files define `pi0_partner_min_energy`, `eta_partner_min_energy`, `pi0_mass_min`, `pi0_mass_max`, `eta_mass_min`, `eta_mass_max`, `missing_energy_min`, and `missing_energy_max` (GeV). Initial pi0/eta thresholds preserve each submit file's old common threshold; the mass windows default to 0.10/0.20 and 0.45/0.65. Missing boundaries default to 0.2/0.5. Bounds must be finite, non-negative, and ordered.

The optional tail of `run_map.sh` after `MIN_CLUSTER_ENERGY_GEV` is:

```text
PI0_PARTNER_MIN_ENERGY_GEV ETA_PARTNER_MIN_ENERGY_GEV PI0_MASS_MIN PI0_MASS_MAX ETA_MASS_MIN ETA_MASS_MAX MISSING_ENERGY_MIN MISSING_ENERGY_MAX
```

`run_small_sample.sh` and `check_sample_map_outputs.sh` accept the same tail after their stored-cluster threshold / existing pi0 threshold position. Omitted pi0 energy defaults to the stored-cluster threshold, omitted eta energy to pi0, and omitted windows/bounds to the values above. All settings are preserved in ROOT metadata and compared when reducing and merging. Use a separate output directory for every configuration.

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
Each submit file defines the stored-cluster/anchor threshold and the independent parameters listed above. The diagnostic floor metadata is fixed to 0 GeV, meaning all positive-energy reconstructed clusters are searched. Choose a new output directory for this schema and never mix configurations in one sample directory.

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

- photon-candidate composition for all six comparison selections;
- pi0-anchor topology for the same six comparison selections.

Both analyses produce `kinematic`, `preselection`, `preselection_tight`, `preselection_isolation`, `region_a`, and `region_a_tagging_veto`. The `region_a` result is before the neutral-meson tag veto; `region_a_tagging_veto` corresponds to the stored `split_cluster_pass_final_photon` flag.

The composition categories partition every selected candidate into prompt photon, pi0 topology, eta, or other. The detailed pi0 categories are separated, merged, single contaminated, missing, and other. Truth-origin majorities use a strict contribution `> 0.5`; an exact contribution of `0.5` is other. An overlap is recorded per selection in metadata and makes the reducer return code 9 after writing its output.

The topology prompt reference requires both the selected flag and `split_cluster_truth_prompt_cluster`. The pi0 categories partition selected valid pi0 anchors. Weighted fraction errors are computed from the weighted numerator/denominator subset covariance; weighted fractions are never added directly.

The merge also computes category-conditional survival fractions relative to the kinematic selection. In each cluster-ET bin, every curve is the weighted yield for one category after the current selection divided by the weighted yield for the same category after the kinematic selection. The composition summary contains all candidates, prompt photons, aggregated pi0-origin candidates, eta-origin candidates, and other-origin candidates; its detailed version resolves the pi0 topology. The anchor-topology summary contains all pi0-main anchors and the five topology categories; its detailed version resolves the missing category.

These survival curves have independent category denominators and therefore are not stacked and do not sum to one. Their uncertainties use the same weighted subset-covariance calculation as the existing fractions. A zero-denominator bin is stored as zero in the ROOT histogram and omitted from the PDF graph. The merger validates the unweighted subset relation and the weighted sum-of-squared-weights subset relation for every selection, category, and bin.

The partial schema is version 5. The merged candidate-composition schema is version 7 and the merged anchor-topology schema is version 5. Both merged metadata trees record the source partial schema, the kinematic denominator, the survival-fraction definition, uncertainty prescription, and zero-denominator convention.

### Shards and partial files

Every non-Jet12 sample has exactly one shard, index 0. Jet12 has exactly ten contiguous, non-overlapping shards, indices 0 through 9. A Jet12 shard reads only its map range but uses the full Jet12 generator-weight sum for cross-section normalization.

Run one partial locally with the mandatory interface:

~~~text
run_reduce.sh FAMILY MAP_ROOT OUTPUT_BASE REQUIRE_COMPLETE N_BINS ET_MAX_GEV SAMPLE_NAME SHARD_INDEX
~~~

For example:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_reduce.sh \
  jet \
  PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p5 \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/reduce/cluster_e_gt_0p5/jet/partial/jet5/shard_0/photon_candidate_selection \
  true 200 40.0 jet5 0
~~~

Each partial is:

~~~text
partial/<sample>/shard_<index>/photon_candidate_selection.root
├── composition/<selection>/
├── anchor_topology/<selection>/
├── metadata
├── composition_summary
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

The submit files default to `configuration = cluster_e_gt_0p5`, `require_complete = true`, 200 bins, and 40 GeV. Every job produces all six selections. Edit the configuration consistently in both submit files before submitting. They write distinct partial and log paths for every sample/shard.

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
run_merge.sh FAMILY PARTIAL_ROOT OUTPUT_BASE
~~~

For example:

~~~bash
PhotonAnalysisTree/workflows/photon_candidate_selection/run_merge.sh \
  jet \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/reduce/ClusterE0p2_Pi0Partner0p2/jet/partial \
  PhotonAnalysisTree/output/plots/photon_candidate_selection/ClusterE0p2_Pi0Partner0p2/jet
~~~

The merger requires shard 0 for each non-Jet12 sample and shards 0--9 for Jet12. It rejects missing or duplicate coverage, invalid shard ranges, inconsistent full-sample normalization, unexpected sample metadata, incompatible analysis/configuration metadata, incompatible axes, and invalid category partitions. It adds only counts and weighted spectra, then recomputes all fractions and errors.

Both plot families share one output base and the same selection directories:

~~~text
<OUTPUT_BASE>/
├── candidate_composition.root
├── pi0_anchor_topology.root
├── kinematic/*.pdf
├── preselection/*.pdf
├── preselection_tight/*.pdf
├── preselection_isolation/*.pdf
├── region_a/*.pdf
└── region_a_tagging_veto/*.pdf
~~~

Each selection contains five composition PDFs (summary, detailed and superdetailed category-fraction stacks, plus summary and detailed survival fractions) and eight topology PDFs. The superdetailed composition stack splits pi0 missing into the six diagnostic categories, normalized to all selected candidates. The six additional fraction histograms are saved in `candidate_composition.root`.

Each selection also contains `workinprogress/region_a_candidate_truth_origin.pdf`. This isolated plot compares the weighted prompt-photon-, pi0-, eta-, and other-candidate spectra after that selection; the five pi0 topology subclasses are summed into one pi0 curve.

Existing schema-5 partials can be reused. Before reusing their topology missing breakdown, the merger checks that missing counts, weighted yields and errors agree with composition in every bin, including underflow and overflow. A mismatch stops the merge rather than displaying a breakdown from a different population.

Both plot families use the same caption layout: collaboration label, sample family, truth vertex cut, candidate kinematics, selection, stored/anchor cluster threshold, and separate pi0/eta tagging energy thresholds and mass windows. Values come from validated production metadata; energies are shown to two decimal places. “Candidate cluster” does not imply an event-leading selection. “Stored/anchor clusters” describes `min_cluster_energy`: lower-energy clusters can still participate in tagging or truth-partner diagnostics, so a global “ignored” label would be inaccurate.

Each ROOT file has one directory per selection. The composition file includes unweighted and weighted spectra, photon purity, category fractions, normalization inputs, and classification QA counters. Both families include summary and detailed `*_survival_fraction_relative_to_kinematic*.pdf` plots.

The weighted anchor-topology spectra use the same logarithmic y-axis range, from `1e-2` to `5e6 pb/GeV`, for every selection so their absolute changes can be compared directly.

The topology ROOT file stores counts, weighted spectra in pb, bin-width-normalized spectra in pb/GeV, fractions, metadata, and per-selection sample summaries.

Production results require `REQUIRE_COMPLETE=true`. `false` is for incomplete QA productions only; normalization then uses only available maps and is not a physics result. Jet and PhotonJet families remain separate and are never mixed.

ABCD purity extraction remains a separate future reduce stage.

## Regression checks

After sourcing ana.565 and building, run `PhotonAnalysisTree/tests/run_candidate_selection_tests.sh`. The tests cover missing energy and mass boundaries, acceptance and incomplete matching, recovery enabled at 50% and disabled, merged/separated/single-contaminated topology, and rejection of mismatched map/partial settings. `python PhotonAnalysisTree/tests/check_candidate_veto.py MAP.root` independently recomputes veto flags and best partners for QA maps whose pi0 and eta thresholds are at least the stored-cluster threshold.

On the first 200 Jet5 events, the previous implementation took 24.08 s wall / 22.32 s user and the recovery-preserving commit `b7eb4d6` took 25.43 s wall / 23.48 s user (one run each, including ROOT startup and I/O). Both used stored-cluster threshold 0.2 GeV and tagging threshold 0.5 GeV. This is a small-sample estimate for the entire change, not an isolated diagnostic-floor benchmark. The new map processed 200 events, wrote 133, rejected 67 vertices, and had no invalid events. Map validation and six-selection reduce passed. A separate 20-event map verified unequal pi0/eta thresholds 0.5/0.7 GeV.

## Recovery comparison checkpoints

Commit `b7eb4d6` contains the configurable partner cuts and new missing categories while preserving the 50% recovery requirement in photon_candidate_selection. The subsequent commit disables that requirement and records topology algorithm 10, so the two productions cannot be mixed. To inspect the first version, use `git switch --detach b7eb4d6`; return to the current version with `git switch codex/photon-selection-partner-settings`. Rebuild the library for the chosen revision.

On the identical 200-event Jet5 input, removing recovery changed 15 of 436 anchor topologies. All 22 compared reconstructed kinematics, tag/veto, partner-selection, and selection-flag branches were identical across the 133 written events. The independent veto check passed all 4666 decisions (1974 tags). `compare_recovery_maps.py BEFORE.root AFTER.root` reproduces this comparison. The direct-deposit match-coverage cut of 50% and the pi0-main contributor requirement are distinct cuts and remain enabled.
