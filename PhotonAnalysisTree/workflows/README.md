# 解析を実行する入口

通常操作する submit、shell script、ROOT マクロと、その実行条件の README / recipe を置きます。
内部ヘッダ・集計実装は [analysis](../analysis/README.md)、生成ファイルの検証は [validation](../validation/README.md) にあります。

| 入口 | 役割 |
|---|---|
| [photon_candidate_selection](photon_candidate_selection/README.md) | 主解析。Pythia DST → map → reduce → merge |
| [bdt_training](bdt_training/README.md) | schema 11 map から inclusive なクラスターと3版の教師ラベルを軽量 ROOT に抽出 |
| [pythia_tree](pythia_tree/README.md) | 汎用 Pythia tree 生成・score 追加 |
| [single_particle](single_particle/README.md) | gun tree 生成・score 追加・dataset 作成 |
| [physics_data](../input/physics_data/README.md) | Run 2 pp 物理データの Spin QA と DST 入力準備 |
| [common](common/README.md) | 共有の manifest 生成・score 追加・過去 study replay |
| [studies](studies/README.md) | 独立解析・保存した過去 study の再実行 |
| [diagnostics](diagnostics/README.md) | topology、tagging 等の補助的な物理調査 |

主解析の map は DST を直接読むため、汎用 Pythia tree の作成は前提ではありません。
実行スクリプトは必要な出力検証を自動実行します。検証実装を `validation/` に分けても、成功判定は省略しません。
