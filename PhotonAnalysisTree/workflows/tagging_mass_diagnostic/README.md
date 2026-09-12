# Tagging invariant-mass diagnostic

Independent reduce/merge workflow for existing photon-candidate maps (schema 5, topology algorithm 10, recovery requirement disabled). No map rebuild is needed. All candidates use the saved **Region A before meson veto** flag; candidate ET is not an event-leading requirement.

## Populations

| ROOT directory | Entries |
|---|---|
| `separated_truth_pair` | Every valid reconstructed representative truth pair of a Region-A truth-pi0 separated candidate, without a partner-energy or mass-window cut |
| `separated_truth_pair_below_threshold` | Subset whose representative partner has E <= the map's pi0 partner threshold |
| `prompt_pi0_all_pairs` | Every Region-A prompt candidate paired with every other stored cluster above the pi0 partner threshold |
| `prompt_pi0_window_pairs` | All such pairs inside the strict pi0 mass window |
| `prompt_eta_all_pairs` | Every Region-A prompt candidate paired with every other stored cluster above the eta partner threshold |
| `prompt_eta_window_pairs` | All such pairs inside the strict eta mass window |

Truth-pi0 means a valid pi0 anchor with main fraction > 0.5; separated is topology code 1. Prompt means the saved prompt flag and dominant fraction > 0.5, matching the candidate-composition analysis. Overlapping origins are rejected.

The representative truth partner is the usable daughter truth match with maximum absolute direct daughter deposit across **all positive-energy clusters**. It is not chosen by reconstructed mass. The saved representative mass can therefore be used even when the partner is below the stored-cluster threshold. The threshold-qualified partner that establishes the separated topology can differ from this representative. Map pi0 partner thresholds consequently affect the **selected separated population**, although no additional energy cut is applied to its representative pair here. Invalid/unavailable/same-cluster representative pairs are counted separately, not silently put into mass underflow.

Prompt pairs are recomputed from stored cluster E, eta and phi using massless four-vectors. Partner selection is strictly E > the corresponding map threshold; it does not impose candidate BDT, isolation or shower cuts. Self-pairs are excluded by cluster identity. These are directed candidate–partner entries: two prompt candidates can each contribute the same unordered pair with their own candidate ET. Multiple window partners all contribute; there is no best-partner selection. The existence of each prompt candidate's pi0/eta window pairs must agree with its saved tag flags, otherwise reduce fails before publishing a result.

## Configuration and input validation

`configuration` selects an input directory, not numerical cut overrides. Thresholds and mass windows are read from map metadata, checked across all files, retained in outputs and used in plot captions/window lines. The pi0 **and** eta partner thresholds must be at least the map stored-cluster threshold, so all eligible combinatorial partners are available. Missing-energy boundaries are retained but not used in these separated/prompt selections. Region-A kinematic/isolation settings, release and BDT model hash are also compared. No thresholds are emulated by this workflow.

Map files must have contiguous manifest coverage starting at zero, sequential chunk IDs, compatible schemas/settings and consistent per-sample cross section/stitching bounds. `require_complete=true` additionally requires the whole sample manifest. Stitching validity/pass and event-weight validity are required. Weights are `weight_numerator_pb / sum(metadata.sum_generator_weight_processed)` for the entire available sample, never per shard. Jet and PhotonJet are processed separately.

## Output and counting

Each population stores raw and cross-section-weighted mass spectra, mass versus candidate ET, and raw/weighted one-dimensional projections for every ET bin. Merge creates four PDFs per population: weighted mass, weighted mass versus ET, and raw/weighted unit-normalized ET-bin overlays.

Mass uses 240 bins from 0 to 1.2 GeV (5 MeV per bin). Candidate ET bins are 0–5, 5–6, 6–8, 8–10, 10–15, 15–20, 20–35, 35–50, 50–100 GeV. Mass overflow is retained and included in shape-normalization denominators; only 0–1.2 GeV is displayed. An ET bin with zero signed normalization is displayed as an empty shape. These fixed ET bins define the saved 2D axes.

`*_window_count` and `*_window_pb` store exact outside/inside classification versus candidate ET, with strict mass boundaries read from metadata. Their `*_fraction` projections are calculated from these counters, including mass overflow, rather than by integrating approximate mass bins. For truth pairs the window is pi0; for prompt groups it is the named meson window. Fractions and uncertainties are recomputed after merging. Histogram Sumw2 and fraction errors use pair-level weight sums and subset covariance; pairs sharing a candidate/event are correlated, so these are not event-clustered uncertainty estimates. Photon veto efficiencies should use candidate counts, not pair totals.

Top-level `candidate_flow_count` and `candidate_flow_pb` separately count separated candidates, valid representative pairs, below-threshold representative partners, invalid pairs, prompt candidates, pi0-vetoed prompt candidates, eta-vetoed prompt candidates, and prompt candidates vetoed by either meson. Pi0/eta overlaps count once in the last category. Every counter is also binned in candidate ET.

Partial metadata records entry ranges, full-sample normalization, settings and QA limits. Merged outputs preserve each source's complete metadata under `sources/partial_<index>/metadata`; normalization denominators are not added across shards. `selection_settings` uses the order documented in `PhotonCandidateSettings.h`. `region_settings` is ordered as shower tower threshold, candidate ET min/max, candidate |eta| max, vertex |z| max, isolation radius/scale/offset and nonisolation gap.

## Condor production

The default configuration is `ClusterE0p12_Pi0Partner0p15`. Both submit files must use the same configuration. Outputs/logs include the configuration name, so the four existing configurations can be reduced separately.

```bash
workflow=PhotonAnalysisTree/workflows/tagging_mass_diagnostic
mkdir -p PhotonAnalysisTree/output/condor/tagging_mass_diagnostic/reduce
condor_submit "$workflow/submit_reduce_jet_samples.job"
condor_submit "$workflow/submit_reduce_jet12_shards.job"
```

There are six ordinary Jet sample jobs and ten contiguous Jet12 entry shards. No script submits jobs automatically. To select another existing configuration, pass the same override to both submissions:

```bash
condor_submit -append 'configuration = ClusterE0p2_Pi0Partner0p5' PhotonAnalysisTree/workflows/tagging_mass_diagnostic/submit_reduce_jet_samples.job
condor_submit -append 'configuration = ClusterE0p2_Pi0Partner0p5' PhotonAnalysisTree/workflows/tagging_mass_diagnostic/submit_reduce_jet12_shards.job
```

Once all 16 jobs succeed:

```bash
PhotonAnalysisTree/workflows/tagging_mass_diagnostic/run_merge.sh jet \
  PhotonAnalysisTree/output/plots/tagging_mass_diagnostic/reduce/ClusterE0p12_Pi0Partner0p15/jet/partial \
  PhotonAnalysisTree/output/plots/tagging_mass_diagnostic/ClusterE0p12_Pi0Partner0p15/jet
```

Merge requires complete map inputs, every expected sample/shard, exact entry coverage, matching per-sample normalization/settings, and matching histogram axes. Existing ROOT outputs are never overwritten. Reduce and merge write temporary ROOT files and publish only after successful writing.

## Local and QA usage

```text
run_reduce.sh FAMILY MAP_ROOT OUTPUT_BASE REQUIRE_COMPLETE SAMPLE_NAME SHARD_INDEX SHARD_COUNT
run_merge.sh FAMILY PARTIAL_ROOT OUTPUT_BASE
```

For limited-event QA, call the C++ reducer directly (the final argument is the event limit; zero means no limit):

```cpp
.L PhotonAnalysisTree/workflows/tagging_mass_diagnostic/ReduceTaggingMassDiagnostic.C
ReduceTaggingMassDiagnostic("jet", "MAP_ROOT", "NEW_OUTPUT", true, "jet8", 0, 1, 1000);
```

For a complete reduction of a partial input sample, use `require_complete=false`. The resulting normalization uses only available maps and is **QA-only**. Targeted QA merge accepts trailing arguments `sample_filter`, `shard_count_override`, `require_complete_partials`, `make_plots`. It still rejects event-limited partials and requires coverage of all entries in the available maps. Production wrappers deliberately use strict full-production merge defaults.

Run the synthetic integration checks after sourcing ana.565:

```bash
python PhotonAnalysisTree/workflows/tagging_mass_diagnostic/test_tagging_mass.py
```

These exercise multiple prompt window partners, independent pi0/eta energy cuts and strict boundaries, the unthresholded representative pair, candidate-versus-pair counts, negative weights, split-versus-serial merge equivalence, PDF generation, overwrite protection, and rejection of incomplete coverage or mismatched settings.

Implementation QA used the first Jet8 map from `ClusterE0p12_Pi0Partner0p15` and `ClusterE0p2_Pi0Partner0p5`: both completed reduce, merge and all 24 PDFs with no prompt veto mismatch. These single-map outputs under `output/qa/tagging_mass_diagnostic/jet8_first_map` are not full-production physics results. Condor production has not been submitted by implementation/testing.
