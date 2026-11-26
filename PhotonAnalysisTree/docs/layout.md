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

## ビルドと実行ソフトウェア

開発用は `build/dev/<version>/` と `install/dev/<version>/`、正式実行用は `software/<software-id>/` に置きます。既定の開発版名は release 名です。
ローカルでは `install/current` が検証済みの開発 install を指します。ビルド記録と production manifest の設計は [software](software.md) を参照してください。

## 成果物の配置

コードを `studies/`・`diagnostics/` に分類しても、成果物の workflow ID と保存先は変更しません。
例えば `workflows/diagnostics/tagging_mass/` の成果物は引き続き `output/tagging_mass_diagnostic/` にあります。

```text
output/<workflow>/
  productions/<production>/<configuration>/
    manifest.json
    maps/       # 個別 tree、または外部 map へのリンク
    partials/   # 必要な workflow のみ
    results/    # 最終 ROOT・PDF・CSV。tree 生成では merged dataset
  qa/
output/condor/<workflow>/  # submit の stdout/stderr の既定先
```

小規模な過去 study は production を `historical` とする。`preserved_results` は元の結果配置を保管した集合であり、単一の物理的 cut を意味しない。個別 partial の設定・sample と対応させる際は ROOT metadata を確認する。

## 主解析と過去参照

- 正式結果は `output/photon_candidate_selection/productions/<production>/<configuration>/results/` を参照する。ローカルに残る production は `common_composition_anchor_v14_eta_parent_v1_ana565`、`newsignal_v1`、`topology_v11`。
- `topology_v11` は4 configuration の過去参照。map schema 5 で、eta topology 拡張前の eta tagging を含む。
- `topology_v11` と `common_composition_anchor_v14_eta_parent_v1_ana565` の `maps/` は tg を指す。`newsignal_v1` の map は production 配下にある。
- `common_composition_anchor_v14_eta_parent_v1_ana565` は既存の map を tg に移し、当時の manifest を保持している。移動先は manifest の `output_roots.map` に登録されていないため、その manifest を指定した通常の reduce では入力検証が止まる。
- 保存 ROOT 内の入力パス・metadata は当時の記録として維持する。前回の成果物移行の対応表は `output/layout_manifest.json` にある。
- 新しい設定で実行するときは、新しい configuration/output を明示する。過去結果を上書きしない。

## ログ

stdout/stderr の既定先は `output/condor/<workflow>/<stage>_<ProcId>.out/.err`。主解析の submit は `log_directory` で production ごとの保存先も指定できる。HTCondor のイベントログはサイト規則に従って固定の `/tmp/ryotaro-PhotonAnalysisTree-<workflow>-<stage>_<ProcId>.log` を使う。
