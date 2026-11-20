# 過去 study の再実行

以下のコマンドは repository root から実行する。元の図・ROOT は移動前後で照合して保持している。再実行先は必ず別の場所を明示し、保存結果を上書きしない。

## Conditional efficiency / reweighted / pair matching

各 workflow の `reproduction.json` に入力 tree と引数を記録している。conditional 系は保存 ROOT の TNamed/TParameter を使用する。reweighted の重みは `pow(truth_pt / 3.0, -8.122)` であり、保存された指数・参照 pT・閾値を configuration ごとに渡す。

```bash
bash PhotonAnalysisTree/workflows/studies/conditional_efficiency/run.sh ET5to7GeV_threshold200MeV --check
bash PhotonAnalysisTree/workflows/studies/conditional_efficiency/run.sh ET5to7GeV_threshold200MeV --output-base /tmp/conditional_replay/result
bash PhotonAnalysisTree/workflows/studies/conditional_efficiency_reweighted/run.sh ET3to5GeV_threshold100MeV_widepT --output-base /tmp/reweighted_replay/result
bash PhotonAnalysisTree/workflows/studies/pair_matching_efficiency/run.sh default --output-base /tmp/pair_matching_replay/result
```

`--check` は入力とレシピの存在を確認する。再実行は ana.565 を使う。pair matching は旧 ROOT に入力/cut metadata がないため、既定の pi0 schema-3 dataset と `DeltaR=0.03`、`0.10<m<0.18`、`E_cluster>0.1 GeV`、asymmetry 20 bins で保存ヒストグラム112個との一致を確認した。代表例の比較で、conditional 1,800個、reweighted 36個、pair matching 112個のヒストグラムが一致した。検証結果は `build/layout_migration/validation.json` に保存する。

## Photon candidate の v11/schema 5

現行 reducer は map schema 10 を入力に取り、partial schema 13 を出力する。v11 の再実行には、partner 0.15 GeV の partial schema 5 用に `5c1a09ba54b091e2f3f4419d8106e4027cc9160f`、partner 0.5 GeV の partial schema 6 用に `f9813670f083bbcf6e3d8116a7d87aae1f868601` を使う。wrapper が必要な ROOT マクロと描画スタイルだけを build 内に展開する。現行のソースや install は変更しない。

```bash
bash PhotonAnalysisTree/workflows/photon_candidate_selection/run_reference.sh merge ClusterE0p12_Pi0Partner0p5 /tmp/v11_reference_replay
bash PhotonAnalysisTree/workflows/photon_candidate_selection/run_reference.sh reduce ClusterE0p12_Pi0Partner0p5 jet5 0 /tmp/v11_partial_replay/result
```

4 configuration とも map と16 partial を保存している。通常は merge だけで再描画できる。元の recovery_manifest.json 内の旧パスも、layout_manifest.json の対応で読み替える。

## Tagging mass / truth pT

v11 tagging mass の入力 map は candidate workflow が所有する。最新定義の map を旧定義用 reducer へ渡さない。truth pT の過去結果は `schema3/prompt_eta07_unweighted_inclusive_pi0_decay` にある。詳細な実行引数は各 workflow の README を参照する。

## Cluster ET と candidate selection の比較

入力 sample、頂点 cut、reco/truth eta、cluster/partner 閾値、起源 majority、topology matching、recovery、weight/正規化を揃えてから比較する。cluster_et は minimum-bias の unweighted study で、過去の contributor 境界や幾何 matching も使う。candidate selection の Jet stitching・断面積重み・厳密な多数寄与・現行 topology と同一条件の図ではない。過去結果はそのまま保存し、名称だけで同一視しない。
