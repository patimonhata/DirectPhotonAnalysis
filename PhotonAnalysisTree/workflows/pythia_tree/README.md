# Pythia analysis tree

Pythia 4-stream DST から tree を生成し、split score を追加する workflow。`run_tree.sh` / `submit_tree.job`、`run_add_scores.sh` / `submit_scores.job`、`check_pythia_tree.C` / `check_pythia_scored_tree.C` をまとめています。

```bash
bash PhotonAnalysisTree/workflows/common/make_pythia_input_manifest.sh
bash PhotonAnalysisTree/workflows/pythia_tree/run_tree.sh pythia8_Jet5-0000000028-000000.root 0
```

個別 tree は `output/pythia_tree/productions/historical/default/maps/`、既存 merged dataset は同 configuration の `results/` に保存します。既存 dataset の schema は ROOT metadata で確認してください。scoring は split ONNX model の明示指定が必須で、入力を上書きしません。同梱 NO_SPLIT 学習モデルの SPLIT への適用は学習条件外です。

`Fun4All_PythiaTruthSpectrumTree.C` と対応 checker は旧 tree-based truth study の補助として保持しています。通常の truth-pT workflow は DST を直接読む `workflows/truth_pt/` です。
