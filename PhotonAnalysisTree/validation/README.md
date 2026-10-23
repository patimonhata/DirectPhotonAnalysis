# 出力検証と小規模 QA

実際の ROOT 出力の schema、branch alignment、metadata、file coverage、保存結果との一致を検証します。
`tests/` の自動テストや `workflows/diagnostics/` の物理調査とは用途を分けています。
本番の実行スクリプトは、従来どおり必要な validator を呼び出してから成功とします。

| 配置 | 対象 |
|---|---|
| `photon_candidate_selection/` | candidate tree/map、sample coverage、veto、recovery 比較、小規模実入力 QA |
| `single_particle/` | gun base/scored tree |
| `pythia_tree/` | Pythia base/scored/truth-spectrum tree |
| `cluster_et/` / `truth_pt/` / `pi0_cluster_multiplicity/` | 各 study の partial |
| `tagging_working_point_scan/` | scan と入力 map / 過去結果との比較 |

repository root からの例:

```bash
bash PhotonAnalysisTree/validation/photon_candidate_selection/check_sample_map_outputs.sh jet5
bash PhotonAnalysisTree/validation/photon_candidate_selection/run_small_sample.sh jet5 1 1
root -l -b -q 'PhotonAnalysisTree/validation/pythia_tree/check_pythia_tree.C("input.root")'
```

coverage 検証と小規模 QA は既定の cut を持ちます。本番とは異なる条件の場合は引数で入力・出力・cut を明示してください。
小規模 QA は `output/photon_candidate_selection/qa/` に map と partial を作り、物理的な full-sample 正規化の代用にはしません。

コード変更時の既存テスト（ana.565 を設定し、ライブラリを build 後）:

```bash
bash PhotonAnalysisTree/tests/run_candidate_selection_tests.sh "$PWD/PhotonAnalysisTree/build/candidate"
```
