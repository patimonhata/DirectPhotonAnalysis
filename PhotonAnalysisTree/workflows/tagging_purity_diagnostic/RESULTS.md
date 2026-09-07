# Photon tagging / purity diagnostic results

## Scope and definitions

This study uses the existing schema-4 Jet map trees for Jet3, 5, 8, 12, 20, 30, and 40. Candidate yields are stitched with
`cross_section_pb * generator_weight / full_sample_sum_generator_weights`. A prompt candidate has a strict prompt truth contribution above 0.5.

The denominator in the threshold study is Region A before the meson veto: kinematic selection, preselection, tight shower selection, and isolation have
already been applied. The final selection rejects the union of pi0 and eta tags. Thus, the prompt and background survival values below are veto-only
survivals, not final-over-kinematic survivals.

The paired diagnostic joins the 0.2 and 0.5 GeV map trees by event UID and cluster ID. Across 16 partials it processed 49,377,040 stitched events and
75,829 Region-A anchors over all ET. Event UID mismatches, missing anchor IDs, and Region-A flag mismatches were all zero.

## ET = 5--6 GeV result

All weighted yields are pb and uncertainties are weighted statistical errors.

| pi0 partner threshold | eta partner threshold | final raw | prompt yield | background yield | purity | prompt survival | background survival | truth-pi0 tag fraction |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.5 | 12,816 | 10,285.95 +/- 531.96 | 38,025.88 +/- 924.09 | 0.2129 +/- 0.0096 | 0.9922 +/- 0.0030 | 0.8852 +/- 0.0058 | 0.0863 +/- 0.0058 |
| 0.2 | 0.2 | 5,987 | 7,433.87 +/- 407.92 | 16,193.65 +/- 603.62 | 0.3146 +/- 0.0143 | 0.7171 +/- 0.0262 | 0.3770 +/- 0.0109 | 0.4476 +/- 0.0123 |
| **0.2** | **0.5** | **10,127** | **10,118.05 +/- 529.95** | **25,829.35 +/- 781.36** | **0.2815 +/- 0.0122** | **0.9760 +/- 0.0054** | **0.6013 +/- 0.0106** | **0.4476 +/- 0.0123** |
| 0.5 | 0.2 | 7,405 | 7,580.44 +/- 410.29 | 24,577.99 +/- 719.61 | 0.2357 +/- 0.0111 | 0.7312 +/- 0.0262 | 0.5722 +/- 0.0111 | 0.0863 +/- 0.0058 |

The common pre-veto Region-A yield is 53,322.73 +/- 1,096.71 pb (14,601 raw), split into 10,367.11 pb prompt and 42,955.62 pb background.
Its purity is 0.1944 +/- 0.0088.

The pi0=0.2, eta=0.5 hybrid is the preferred tested working point:

- Relative to 0.5/0.5, prompt yield is 0.9837 +/- 0.0045, background yield is 0.6793 +/- 0.0110, and purity increases by
  0.0686 +/- 0.0040. These errors include the candidate overlap covariance between working points.
- Relative to 0.2/0.2, the hybrid recovers a prompt-yield factor of 1.3611 +/- 0.0496 while accepting a background factor of
  1.5950 +/- 0.0378; purity is lower by 0.0332 +/- 0.0092.
- Using the existing kinematic denominator, hybrid final-over-kinematic survival is 0.7138 +/- 0.0185 for prompt,
  0.1188 +/- 0.0033 for background, and 0.1552 +/- 0.0036 for all candidates.

## Why the common 0.2 GeV threshold loses prompt photons

At 0.2/0.2 the prompt veto loss is 2,933.24 pb. Eta-only tags account for 2,720.90 pb (92.8%) and pi0+eta overlap adds another 37.98 pb.
Therefore the large prompt loss is driven by low-energy eta combinatorics, not by the low pi0 threshold itself.

For the preferred pi0=0.2, eta=0.5 hybrid, the ET=5--6 GeV flow is:

| population | pi0 only veto | eta only veto | both veto | survives |
|---|---:|---:|---:|---:|
| prompt | 209.66 pb (246 raw) | 36.71 pb (132 raw) | 2.68 pb (13 raw) | 10,118.05 pb (5,198 raw) |
| background | 15,006.86 pb (3,130 raw) | 1,672.34 pb (656 raw) | 447.08 pb (297 raw) | 25,829.35 pb (4,929 raw) |

The low pi0 threshold raises the truth-pi0 tag fraction from 8.63% to 44.76%, while the high eta threshold avoids most of the prompt loss.

## Exploratory selected-partner quality scan

The stored selected pi0 partner IDs join to the cluster arrays with zero failures. A one-variable scan on the 0.2 GeV output, integrated over anchor ET,
found the following cuts at approximately 90% weighted true-partner efficiency:

| variable and accepted side | true-partner efficiency | prompt-combinatorial rejection |
|---|---:|---:|
| mass <= 0.165 GeV | 0.9126 | 0.5511 |
| deltaR <= 0.13 | 0.9039 | 0.4555 |
| energy asymmetry <= 0.93 | 0.9083 | 0.4200 |
| partner ET >= 0.20 GeV | 0.9453 | 0.2983 |
| partner E >= 0.22 GeV | 0.9200 | 0.2960 |

These are exploratory, histogram-binned, one-dimensional results among already selected partners. They are not an ET=5--6 GeV retagging result and
do not include correlations. The stored BDT score is descriptive only because the low-energy partners may be outside its training range.

## Recommendation and remaining work

Use separate partner thresholds and take pi0=0.2 GeV, eta=0.5 GeV as the current working point. Before changing the production default, recompute
pairing from the stored cluster four-vectors to scan intermediate eta thresholds (for example 0.3 and 0.4 GeV) and ET-specific mass/deltaR/asymmetry
cuts. Validate the chosen point on an independent sample.

The map trees are sufficient for cluster-level re-pairing at thresholds at or above 0.2 GeV. Their generic all-pair arrays are empty, and tower-key
detail is absent. A detector-quality study involving bad/dead-tower proximity requires additional production.

## Outputs

- `hybrid_threshold/jet/hybrid_threshold_diagnostic.root`: merged paired diagnostic, including inter-working-point overlap histograms.
- `hybrid_threshold/jet/hybrid_threshold_diagnostic_flow.tsv`: raw and weighted count flow.
- `hybrid_threshold/jet/hybrid_threshold_diagnostic_metrics.tsv`: purity, survival, and truth-pi0 metrics.
- `hybrid_threshold/jet/hybrid_threshold_diagnostic_comparisons.tsv`: correlated working-point comparisons.
- `cluster_e_gt_0p2/jet/selected_partner_feature_scan_90pct.tsv`: exploratory selected-partner scan.
- `cluster_e_gt_0p2/jet/selected_partner_feature_roc.tsv`: all one-dimensional ROC threshold points.
