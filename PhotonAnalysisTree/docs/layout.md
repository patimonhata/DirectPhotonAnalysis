# PhotonAnalysisTree の配置

コードは `workflows/<workflow>/`、成果物は `output/<workflow>/` に置く。旧 `macro` と、共通の `output/plots`・`root`・`merged`・`intermediate_files`・`qa` は廃止する。旧パスへの互換リンクは作らない。

```text
workflows/<workflow>/
output/<workflow>/
  productions/<production>/<configuration>/
    manifest.json
    maps/       # tree 生成では個別 tree、candidate 解析では外部 map へのリンク
    partials/   # 必要な workflow のみ
    results/    # 最終 ROOT・PDF・CSV。tree 生成では merged dataset
  qa/
output/condor/<workflow>/  # configuration/production に依存しない固定ログ
```

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

## 移行時の検証

保持データ30,314件を照合し、うち1,672件は SHA-256 も確認した。外部 map リンク5件、28 submit の dry-run、ライブラリ build、既存 candidate pipeline と tagging mass のテストを確認した。代表的な過去 study の1,948 histograms と v11 両 partial 世代の1,836 histograms は再実行と保存結果が一致した。working-point scan の8,400比較も通過した。

`qa/other_origin_categories` は、当時の schema-7 map への参照が移行前から欠けているため、QA manifest に `historical_input_missing` と記録して保存した。正式結果の入力は保持・検証している。
