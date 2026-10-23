# Prompt efficiency / purity projection

`estimate_prompt_efficiency_purity.py` は保存した candidate composition から、preselection 後の prompt TightBDT 効率を
100% と仮定した場合の purity を計算します。学習済みモデルの性能予測や統計誤差の評価ではありません。
背景収量は固定し、ET bin ごとに名目上の prompt isolation / veto acceptance を使います。

既定入力は保存 v11 の `ClusterE0p12_Pi0Partner0p5` 結果です。ana.565 の PyROOT 環境で、repository root から実行します。

```bash
python3 PhotonAnalysisTree/workflows/diagnostics/prompt_efficiency_projection/estimate_prompt_efficiency_purity.py --help
```

入力・出力・ET 範囲を明示する場合は `--help` の引数を使ってください。
