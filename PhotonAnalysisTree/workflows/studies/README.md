# 独立解析・過去 study

主解析とは独立した解析を実行する入口です。過去 study も保存しています。
[再実行の説明](../../docs/reproduction.md) と各 README で入力・解析条件を確認してください。
コードの分類を変更しても、成果物は既存の `output/<study>/` に保存します。

| Study | 状態・入口 |
|---|---|
| [conditional_efficiency](conditional_efficiency/README.md) | 保存結果と再現レシピあり |
| [conditional_efficiency_reweighted](conditional_efficiency_reweighted/README.md) | 保存結果と再現レシピあり |
| [pair_matching_efficiency](pair_matching_efficiency/README.md) | 保存結果と再現レシピあり |
| [cluster_et](cluster_et/README.md) | 過去 study。candidate と同じ物理条件の図として扱わない |
| [truth_pt](truth_pt/README.md) | 過去結果を保存。最新定義の結果は未生成 |
| [pi0_cluster_multiplicity](pi0_cluster_multiplicity/README.md) | 独立した小規模解析 |
| [score_comparison](score_comparison/README.md) | 保存 PDF と専用マクロあり。入力再現性の整備は未完了 |

`run.sh ... --check` を持つ study は、ROOT を実行せず入力とレシピの存在を確認できます。
