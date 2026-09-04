# Displaced-partner invariant-mass diagnostic

This workflow reads schema-4 photon-candidate maps produced with the strict topology threshold `E_cluster > 0.5 GeV`. It selects the representative truth-associated pair classified as `Missing: displaced partner cluster`, requires the saved representative partner energy to be strictly above an emulated tagging threshold (0.2 GeV for this study), and plots its saved reconstructed invariant mass.

It does not recompute the full combinatorial pi0 veto: maps retain one representative diagnostic partner per truth photon, not every reconstructed split cluster below 0.5 GeV.

## Condor production

The Jet production is split into 16 Reduce jobs: one job each for Jet3, 5, 8, 20, 30, and 40, plus ten entry-range shards for Jet12. Each worker writes ROOT histograms only; plotting happens once after Merge.

~~~bash
workflow=/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/workflows/displaced_partner_mass_diagnostic
logs=/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/condor/displaced_partner_mass_diagnostic/reduce
mkdir -p "$logs"
cd "$workflow"
condor_submit submit_reduce_jet_samples.job
condor_submit submit_reduce_jet12_shards.job
~~~

After all 16 jobs finish:

~~~bash
partial_root=/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/plots/displaced_partner_mass_diagnostic/reduce/cluster_e_gt_0p5/jet/partial
output_base=/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/plots/displaced_partner_mass_diagnostic/cluster_e_gt_0p5/jet
./run_merge.sh jet "$partial_root" "$output_base"
~~~

Merge refuses incomplete, overlapping, or mismatched shards. It requires complete input maps, exact entry coverage, schema 4, topology threshold 0.5 GeV, production tagging threshold 0.5 GeV, diagnostic floor 0.1 GeV, emulated tagging threshold 0.2 GeV, a common release/model hash, and consistent per-sample normalization.

To retry a subset, override `sample_names` or `shard_indices`, for example:

~~~bash
condor_submit -append 'sample_names = jet8,jet20' submit_reduce_jet_samples.job
condor_submit -append 'shard_indices = 3,7' submit_reduce_jet12_shards.job
~~~

## Serial and QA runs

The serial entry point remains useful for one sample or a limited event check:

~~~text
run.sh FAMILY MAP_ROOT OUTPUT_BASE REQUIRE_COMPLETE TAGGING_PARTNER_MIN_ENERGY_GEV [SAMPLE_FILTER] [MAX_EVENTS_PER_SAMPLE]
~~~

Example with one thousand Jet8 events:

~~~bash
PhotonAnalysisTree/workflows/displaced_partner_mass_diagnostic/run.sh \
  jet \
  /sphenix/tg/tg01/coldqcd/ryotaro/DirectPhotonAnalysis/photon_candidate_selection/cluster_e_gt_0p5 \
  PhotonAnalysisTree/output/qa/displaced_partner_mass_diagnostic/jet8_1000events \
  false 0.2 jet8 1000
~~~

`MAX_EVENTS_PER_SAMPLE > 0` is QA-only. Normalization still uses the full matched sample sum of generator weights, so a limited-event yield is not a physics result.

`run_reduce.sh` accepts a sample name plus `SHARD_INDEX SHARD_COUNT`; it partitions the sample TChain into contiguous entry ranges. `MergeDisplacedPartnerMassDiagnostic.C` also has optional sample/shard arguments for targeted QA, while `run_merge.sh` deliberately uses the strict full-production defaults.

## Selection and output

Selected pairs must have topology `missing` (3), missing category `displaced_partner_cluster` (4), alignment `displaced` (2), and production truth-pair tag status `below_energy_threshold` (4). The saved partner energy must satisfy `E_partner > 0.2 GeV`. The macro also requires the saved diagnostic-pair and representative-truth-pair masses to agree.

`displaced_partner_mass_diagnostic.root` stores raw and cross-section-weighted invariant-mass spectra, partner-energy spectra, mass versus anchor ET, mass versus truth-pi0 pT, all one-dimensional projections, metadata, and raw/weighted pi0-window fractions. Directories for inclusive, kinematic, preselection, preselection plus TightBDT, preselection plus isolation, and Region A selections are included. No production meson-tag veto is applied.

PDFs show the inclusive weighted spectrum, weighted two-dimensional distributions, and raw-count unit-normalized projections. Red dashed lines mark the strict `0.10 < mass < 0.20 GeV` window. Anchor-ET bins are `0-5, 5-6, 6-8, 8-10, 10-15, 15-20, 20-35, 35-50, 50-100 GeV`; truth-pi0-pT bins are `0-3, 3-5, 5-6, 6-8, 8-10, 10-15, 15-20, 20-35, 35-50, 50-100 GeV`. The saved two-dimensional histograms allow arbitrary later rebinning.
