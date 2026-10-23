# PhotonAnalysisTree の配置

`workflows/` は人間が実行する入口、`analysis/` は ROOT による集計・描画の内部実装、
`src/` は Fun4All module と共有の物理処理です。生成ファイルの検証は `validation/`、自動テストは `tests/` に分けます。

```text
src/                                      # DST 読み込み・truth matching・topology・推論
analysis/
  photon_candidate_selection/             # 共通型・histogram・描画・reduce/merge 実装
  common/                                 # 解析条件検証・描画スタイル
  tagging_mass/                           # diagnostic 内部実装
  tagging_working_point_scan/
workflows/
  photon_candidate_selection/             # 主解析: submit・run・実行用 ROOT マクロ
  pythia_tree/                            # Pythia tree 生成・score 追加
  single_particle/                        # gun tree 生成・score 追加
  common/                                 # manifest・scoring・replay の共通入口
  studies/<study>/                        # 独立解析・過去 study の実行入口
  diagnostics/<topic>/                    # 主解析を調べる補助解析の実行入口
validation/<workflow>/                    # 出力検証・比較・小規模 QA
tests/                                    # 自動テスト（synthetic 入力を含む）
```

実行入口の一覧は [workflows](../workflows/README.md)、
内部部品の依存関係は [analysis](../analysis/README.md)、検証の使い方は [validation](../validation/README.md) を参照してください。

## 成果物の配置

コードを `studies/`・`diagnostics/` に分類しても、成果物の workflow ID と保存先は変更しません。
例えば `workflows/diagnostics/tagging_mass/` の成果物は引き続き `output/tagging_mass_diagnostic/` にあります。

```text
output/<workflow>/
  productions/<production>/<configuration>/
    manifest.json
    maps/       # tree 生成では個別 tree、candidate 解析では外部 map へのリンク
    partials/   # 必要な workflow のみ
    results/    # 最終 ROOT・PDF・CSV。tree 生成では merged dataset
  qa/
output/condor/<workflow>/  # configuration/production に依存しない固定ログ
```

旧 `macro` と共通の `output/plots`・`root`・`merged`・`intermediate_files`・`qa` は前回の成果物整理で廃止済みです。
旧コードパスへの互換リンクは作りません。保存 ROOT 内の過去のパスは履歴として維持します。

小規模な過去 study は production を `historical` とする。`preserved_results` は元の結果配置を保管した集合であり、単一の物理的 cut を意味しない。個別 partial の設定・sample と対応させる際は ROOT metadata を確認する。

## 主解析と過去参照

- 正式結果の入口: `output/photon_candidate_selection/current/jet/`。この新しい入口は `productions/topology_v12_eta_v2_decay_v1/ClusterE0p12_Pi0Partner0p5/results/jet/` を指す。
- 過去参照: `output/photon_candidate_selection/productions/topology_v11/` の4 configuration。map schema 5、eta topology 拡張前。eta tagging 自体は存在する。
- 各 configuration の `maps` は cold-storage の実体へのリンク。実体の移動・削除は行っていない。
- 保存 ROOT 内の入力パス・metadata は当時の記録として維持する。移動対応は `output/layout_manifest.json` に記録する。
- 新しい設定で実行するときは、新しい configuration/output を明示する。過去結果を上書きしない。

## Workflow の扱い

| Workflow | 方針 |
|---|---|
| photon_candidate_selection | 主解析。現行と v11/schema 5 を保存 |
| conditional_efficiency / conditional_efficiency_reweighted | 重要な過去 study。保存 ROOT から再現レシピを抽出 |
| pair_matching_efficiency | 重要な過去 study。旧 ROOT に入力/cut metadata がないため、旧マクロの既定値との再実行比較で確認 |
| cluster_et | 過去 study を保存。candidate selection の類似図と同一条件とは扱わない |
| tagging_mass_diagnostic | v11 本番を保存。最新定義への対応は今回行わない |
| truth_pt | 過去結果を保存。最新定義への対応は今回行わない |
| pi0_cluster_multiplicity | 小さな独立解析として保存 |
| score_comparison | PDF 6枚と専用マクロを保存。入力の探索・再現性の整備は対象外 |
| tagging_purity_diagnostic / tagging_working_point_scan | 低優先の調査結果。旧入力が現存するとは限らず、再実行可能性は未保証 |
| single_particle / pythia_tree | 共有入力 dataset の生成と保存 |

`pi0_anchor_topology` と `displaced_partner_mass_diagnostic` の専用 workflow・専用成果物・ログ、および tagging mass の旧 QA は削除対象。共有 `Pi0AnchorTopologyEvaluator` と EventDisplay は維持する。

## 将来の新定義

`output/tagging_mass_diagnostic/productions/pending_latest_definition/results/` と `output/truth_pt/productions/pending_latest_definition/results/` は空の予約先。manifest の状態は `not_generated` であり、正式結果ではない。定義を確定したときに production/configuration を決める。

## ログ

submit ごとに同じ `output/condor/<workflow>/<stage>_<ProcId>.out/.err` を使い、別 production/configuration でも上書きする。HTCondor のイベントログはサイト規則に従って固定の `/tmp/ryotaro-PhotonAnalysisTree-<workflow>-<stage>_<ProcId>.log` を使う。

## コード配置変更の検証

今回の整理では保存成果物を移動せず、実行入口と内部実装・検証・補助解析の配置を変更した。
Reduce / Merge の計算本体は変更前と文字列単位で一致する。

- 既存の candidate テストと tagging mass の synthetic integration test を共通ランナーで実行し、すべて成功した。
- 変更前後の synthetic schema-8 map → reduce → merge を比較した。ROOT 26ファイルの histogram 12,654個について全 bin の内容・誤差、軸・label・entry 数が一致し、44 tree の branch・182行の内容も一致した。出力先に由来するパス文字列だけを正規化した。
- CSV 15ファイルは byte 単位で一致し、PDF 246ファイルは相対配置とページ数が一致した。PDF の描画画像そのものの一致検証は行っていない。
- 実 DST 1 segment・2イベント（97 cluster）で、小規模 QA の map 生成 → 移動後の validator → reduce が正常終了した。QA 出力は `/tmp/` の新規ディレクトリを使った。
- 過去 study の12 configuration で入力・再現レシピの存在を確認し、全28 submit file の dry-run を確認した。ジョブは投入していない。
- include、Markdown link、shell の参照、Condor の実行先を確認した。既存の submit / recipe の成果物パスは維持した。

テストの再実行（ana.565 を設定し、ライブラリを build 後）:

```bash
bash PhotonAnalysisTree/tests/run_candidate_selection_tests.sh "$PWD/PhotonAnalysisTree/build/candidate"
```

## 前回の成果物移行時の検証

保持データ30,314件を照合し、うち1,672件は SHA-256 も確認した。外部 map リンク5件、28 submit の dry-run、ライブラリ build、既存 candidate pipeline と tagging mass のテストを確認した。代表的な過去 study の1,948 histograms と v11 両 partial 世代の1,836 histograms は再実行と保存結果が一致した。working-point scan の8,400比較も通過した。

`qa/other_origin_categories` は、当時の schema-7 map への参照が移行前から欠けているため、QA manifest に `historical_input_missing` と記録して保存した。正式結果の入力は保持・検証している。
