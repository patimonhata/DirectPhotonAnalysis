# Tagging invariant-mass diagnostic

保存する本番は map schema 5 / topology v11 の `ClusterE0p12_Pi0Partner0p5`。v10 本番・旧 QA は削除済み。最新定義への更新や新しい物理結果の生成は今回行っていない。

```text
output/tagging_mass_diagnostic/productions/topology_v11/ClusterE0p12_Pi0Partner0p5/
  partials/jet/<sample>/shard_<index>/tagging_mass_diagnostic.root
  results/jet/tagging_mass_diagnostic.root
```

元 map は `output/photon_candidate_selection/productions/topology_v11/ClusterE0p12_Pi0Partner0p5/maps/<sample>/` を参照する。元の shared input を削除しない。current candidate の schema-8 map を渡してはいけない。

## Populations and counting

Region A before meson veto を使用する。truth pi0 は有効 anchor かつ単一 parent の main fraction >0.5。prompt は保存 flag かつ dominant fraction >0.5。origin overlap は拒否する。

| Population | Entries |
|---|---|
| separated_truth_pair | separated 候補の代表 truth pair。追加の partner E / mass cut はかけない |
| missing_energy_band_truth_pair | saved missing category が energy-band mass-window in/out の候補 |
| separated_or_missing_energy_band_truth_pair | 上記2群の和。1候補に代表 pair 1個 |
| prompt_pi0_all_pairs / prompt_pi0_window_pairs | prompt と保存 cluster の全組合せ / strict pi0 window 内 |
| prompt_eta_all_pairs / prompt_eta_window_pairs | 同様の eta 組合せ |

代表 partner は正の direct daughter deposit が最大の usable cluster。mass で最適化しない。prompt の組合せは directed pair で、自己 pair を除き、複数の window partner をすべて数える。pair 数を候補の veto 効率と同一視しない。window pair の有無は保存 tag flag と一致させる。

閾値・mass window・missing band は map metadata から読み、全 map の一致、manifest 完全性、sample/shard coverage、full-sample generator-weight 正規化を確認する。v11 の保存設定は cluster E>0.12、pi0/eta partner E>0.5 GeV、pi0 0.10<m<0.20、eta 0.45<m<0.65 GeV、missing band 0.12<E<=0.5 GeV。recovery requirement は無効。

raw count と cross-section-weighted mass / mass-vs-ET を保存する。mass は0–1.2 GeVを240 bins、ET は0,5,6,8,10,18,30,40 GeV。shape overlay は5–6、6–8、8–10 GeV。overflow を保ち、正規化は merge 後に行う。誤差は pair-level Sumw2 / subset covariance で、event-level の相関は推定しない。7群×4 PDFを出力する。

## Replay

既存16 partial から再描画する場合は、新しい出力先を指定する。

```bash
bash PhotonAnalysisTree/workflows/diagnostics/tagging_mass/run_merge.sh jet \
  PhotonAnalysisTree/output/tagging_mass_diagnostic/productions/topology_v11/ClusterE0p12_Pi0Partner0p5/partials/jet \
  /tmp/tagging_mass_v11_replay
```

map からの再集計は `run_reduce.sh FAMILY MAP_ROOT OUTPUT_BASE REQUIRE_COMPLETE SAMPLE SHARD_INDEX SHARD_COUNT`。非Jet12は1 shard、Jet12は10 shard。既存 ROOT は上書きしない。production は `require_complete=true` を使う。

```bash
python3 PhotonAnalysisTree/tests/test_tagging_mass.py
```

全 sample/shard・設定・軸を検証する。各 partial の metadata は merged ROOT の sources に保持する。有限 event 数や不完全 map による QA は本番結果として扱わない。

将来の新定義用は `output/tagging_mass_diagnostic/productions/pending_latest_definition/results/` に予約し、`not_generated` として明示する。
