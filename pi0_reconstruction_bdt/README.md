# Pi0Reconstruction BDT adapter

This package adds the PPG15/PPG12 `base_v3E` SPLIT photon-identification BDT
score to the vector-based `event_tree` written by `Pi0Reconstruction`.

The input ROOT file is never modified. A separate output ROOT file is created
with the original tree and histograms plus BDT branches and provenance
metadata.

## Quick start

The sPHENIX ROOT environment must provide ROOT 6.32 or newer with TMVA `RBDT`.

```bash
tar -xzf pi0_reconstruction_bdt_adapter.tar.gz
cd pi0_reconstruction_bdt
./run_add_bdt.sh INPUT.root OUTPUT.root
```

The wrapper uses this model by default:

```text
/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root
```

A different model file can be passed as the third argument:

```bash
./run_add_bdt.sh INPUT.root OUTPUT.root /path/to/model.root
```

The wrapper checks the input and model paths, refuses to overwrite either the
input or an existing output, runs the BDT adapter, and validates the resulting
branch sizes and score metadata. Remove an old output explicitly before
re-running the command.

## Output branches

Both branches are vectors with exactly `ncluster` elements in each event.

- `cluster_bdt_base_v3E_split` (`std::vector<float>`): photon-ID BDT score.
- `cluster_bdt_base_v3E_split_valid` (`std::vector<unsigned char>`): `1` when
  `cluster_shower_valid` is set and all 11 BDT inputs are finite; otherwise `0`.

The score branch is still filled when a non-finite feature is replaced by zero,
matching the behavior of the existing PhotonAna score-adder. Such a score has
`valid=0` and should normally be excluded from analysis.

This is a photon-identification BDT score. It is not an NPB score.

## BDT definition

- Model: `model_base_v3E_split_single_tmva.root`
- Model key: `myBDT`
- Feature order:

```text
ET
weta_cogx
wphi_cogx
vertex_z
cluster_eta
e11_over_e33
et1
et2
et3
et4
e32_over_e35
```

The model path, model key, feature order, processing counts, compatibility
warning, and kinematic warning are stored as top-level ROOT objects in every
output file.

## Required input schema

The input must contain `event_tree` with these branches and types:

| Branch | Type |
|---|---|
| `ncluster` | `UInt_t` |
| `vertex_z` | `Double_t` |
| `cluster_et` | `std::vector<double>` |
| `cluster_eta` | `std::vector<double>` |
| `cluster_shower_valid` | `std::vector<unsigned char>` |
| `cluster_shower_w_eta_cogx` | `std::vector<float>` |
| `cluster_shower_w_phi_cogx` | `std::vector<float>` |
| `cluster_shower_e11_over_e33` | `std::vector<float>` |
| `cluster_shower_e32_over_e35` | `std::vector<float>` |
| `cluster_shower_et1` through `cluster_shower_et4` | `std::vector<float>` |

Other branches and top-level histograms are copied unchanged.

## Test data

- `inputs/ForJaein_pi0_reconstruction_SPLIT_000000.root`: preserved 20-event
  SPLIT-cluster smoke input.
- `output/ForJaein_pi0_reconstruction_SPLIT_000000_with_bdt.root`: reference
  output produced by this adapter.

Reference validation result:

```text
Events: 20
Clusters: 35
Valid scores: 35
Malformed events: 0
Valid score range: 0.000140803 to 0.932545
```

Run the checker independently with:

```bash
root -l -b -q 'macros/check_bdt_output.C("OUTPUT.root")'
```

## Physics-use limitations

1. The test input reports shower-shape algorithm version 1, a 7x7 tower patch,
   and a 0.070 GeV tower threshold. However, the producer's
   `ShowerShapeCalculator.cc` and `.h` were owner-only (`0600`) during this
   check. Their formulas could not be compared line-by-line with the variables
   used by PhotonAna. Confirm this compatibility with the producer before using
   the scores in a final result.
2. The 35 clusters in the supplied smoke file have
   `0.050 < cluster ET < 4.844 GeV`. The model's documented performance bins
   begin at 6 GeV. These test-file scores therefore validate the software and
   branch format, not the BDT performance in this low-ET range.
3. The input macro used `VertexMode::Origin`, so `vertex_z` is zero in this test
   sample. Confirm that this is the intended vertex treatment for the target
   production.

## Package contents

```text
README.md
SHARE_MESSAGE.md
run_add_bdt.sh
macros/add_ppg15_bdt_to_pi0_tree.C
macros/check_bdt_output.C
inputs/ForJaein_pi0_reconstruction_SPLIT_000000.root
output/ForJaein_pi0_reconstruction_SPLIT_000000_with_bdt.root
```
