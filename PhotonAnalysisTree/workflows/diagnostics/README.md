# 主解析を調べる補助解析

特定の物理的な疑問を調べるための入口です。主解析の map → reduce → merge には含まれません。
ファイル形式の検証は [validation](../../validation/README.md) を使います。

| 診断 | 役割・状態 | 既存の成果物 workflow ID |
|---|---|---|
| [eta_decay_flight](eta_decay_flight/README.md) | DST から eta の飛行距離・崩壊を独立調査 | `photon_candidate_selection/qa/eta_decay_flight` |
| [tagging_mass](tagging_mass/README.md) | 保存した v11 の tagging mass 調査 | `tagging_mass_diagnostic` |
| [tagging_purity](tagging_purity/README.md) | 過去の低優先調査。旧 map の再実行可能性は未保証 | `tagging_purity_diagnostic` |
| [tagging_working_point_scan](tagging_working_point_scan/README.md) | 過去の閾値・質量窓比較。旧入力の再実行可能性は未保証 | `tagging_working_point_scan` |
| [prompt_efficiency_projection](prompt_efficiency_projection/README.md) | 保存 v11 結果に基づく仮定上の効率・purity 予測 | `photon_candidate_selection/.../results/jet/prompt_efficiency_projection` |

v11 用の診断に現行 schema-8 map を渡さないでください。最新定義への対応とは別に、今回の整理では既存の入力条件を維持しています。
