# 出力検証と小規模 QA

実際の ROOT 出力の schema、branch alignment、metadata、file coverage、保存結果との一致を検証します。
`tests/` の自動テストや `workflows/diagnostics/` の物理調査とは用途を分けています。
本番の実行スクリプトは、従来どおり必要な validator を呼び出してから成功とします。

| 配置 | 対象 |
|---|---|
| `calorimeter_acceptance/` | tower geometry 境界、独立角度計算、vertex scan、PDF ページ数 |
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
bash PhotonAnalysisTree/tests/run_candidate_selection_tests.sh "$PWD/PhotonAnalysisTree/build/dev/ana.565/photon"
```

## Common composition の移行・matching 診断

`photon_candidate_selection/diagnose_common_composition.py OLD_PRODUCTION NEW_PRODUCTION NEW_QA_DIRECTORY`
は既存 production を読み取り、全 partial の分類移行を sample・selection・ET 別に集計します。
merged CSV との一致、選択分母の新旧比較、指定 map のイベント identity と再構成 branch の一致を確認します。
既定では Jet3/12/40 の map 0 の先頭1,000 written events を比較します。
`--samples jet3 jet5 jet8 jet12 jet20 jet30 jet40 --map-indices 0 100` で対象を広げられます。
`migration_rates.csv` の割合の分母は **旧 Other** です。`representative_clusters.jsonl` は source file/event/cluster ID と ancestry を含みます。
全 production の移行は新 map 内で並行計算した legacy origin に対する値で、実際の旧 map との直接比較は指定した標本に限ります。
map の標本集計には stitching と kinematic selection を適用し、全 sample の規格化には既存 partial の重みを使います。
標本の先頭イベントは無作為抽出ではなく、全 production の物理割合の推定には使いません。

`photon_candidate_selection/run_common_matching.sh FROZEN_SOFTWARE_ROOT INPUT_MANIFEST MANIFEST_ROW NEW_QA_DIRECTORY [N_EVENTS=30]`
は固定ライブラリの evaluator を使う独立した DST 診断です。診断条件は今回の Jet production と同じ
cluster E>0.12 GeV、partner E>0.5 GeV、|vertex z|<=60 cm、recovery 無効です。
通常の map producer、BDT、topocluster isolation は実行せず、matching が使用する CEMC truth cells/towers を構築します。
各評価済み候補について、閾値以上・未満それぞれの最初の失敗、失敗件数、anchor 自身の成否、保存 topology を CSV に出力します。
閾値未満の再計算も含め、診断で再構成した AND 条件が production evaluator の値と一致しなければ失敗します。
DST 診断には stitching/kinematic selection を適用しないため、map 標本・最終ヒストグラムとは母集団が異なります。

manifest row 0 の診断は、以下で保存済み map 0 と全 candidate identity / matching flag / anchor topology を照合できます。

```text
python3 PhotonAnalysisTree/validation/photon_candidate_selection/summarize_common_matching.py \
  NEW_PRODUCTION NEW_SUMMARY_JSON jet3=JET3_DIAGNOSTIC_DIRECTORY jet40=JET40_DIAGNOSTIC_DIRECTORY
```

出力する代替 gate の件数は保存済み判定に対する診断値です。特に anchor のみを確認する条件は partner の未検出を保証せず、
両 photon の代表 cluster を要求する条件は本来の missing を落とします。これらは当初の調査時点の比較であり、その後の採用条件は以下の schema-10 検証を参照してください。
各スクリプトは新しい出力先を要求し、production や解析ライブラリを書き換えません。


## Anchor matching 適用条件の回帰検証（map schema 10）

`photon_candidate_selection/compare_anchor_matching.py OLD_MAP NEW_MAP DST_DIAGNOSTIC_DIRECTORY NEW_SUMMARY_JSON`
は schema 9 / pi0 v13 / eta v3 と schema 10 / pi0 v14 / eta v4 の同一 source-file/event を比較します。
新 map に保存された全イベントについて、適用可否・理由を除く既存の全 event branch を完全一致で検証します。
元の DST 診断 CSV の anchor matching・projection・interaction から独立に期待値を計算し、新しい適用結果と照合します。
入力は同じ manifest row 0 と Jet 設定で、DST 診断と新 map は同じイベント範囲を使ってください。
新しい anchor 診断3枝も CSV と照合し、非評価 cluster の適用結果が変わっていないことを確認します。

今回の小規模検証は `output/photon_candidate_selection/qa/anchor_matching_v14/` に保存しています。
Jet3/Jet40 各30入力イベント（各21 written events）で、既存464枝が完全一致し、独立診断の423 anchor のうち422件が評価可能、
coverage 約0.452の1件が unknown のままでした。これは少数イベントの動作検証であり、全 production の割合の推定ではありません。
開発ビルド・テスト・map/reduce のログ、map と partial、比較 JSON を同じ QA directory に置いています。
正式 production は clean-source software を新規 freeze し、新規出力先で map から再実行してください。
