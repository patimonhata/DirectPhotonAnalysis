# マクロ入出力一覧

パスと既定値は現行コードを基準にしています。ROOT 出力内の詳しい Tree / branch は
[ROOT TTree ブランチ一覧](root-tree-schema.md) を参照してください。

## 共通事項

- 現行プロジェクトルート: `/sphenix/user/ryotaro/DirectPhotonAnalysis`
- `processID` は通常6桁ゼロ埋めされます（例: `0` → `000000`）。
- `n_events=0` を Fun4All に渡すマクロでは入力の全イベントを処理します。
- ROOT file を `RECREATE` で開くマクロは同名ファイルを上書きします。
- `CLUSTERINFO_CEMC` は split、`CLUSTERINFO_CEMC_NO_SPLIT` は no-split cluster node です。

## SinglePi0GunSimulation

### `Fun4All_SingleParticlePi0.C`

| 項目 | 内容 |
| --- | --- |
| 目的 | 単一 pi0 を生成し、Geant4、calorimeter reconstruction、追加の no-split clustering を実行する |
| 引数 | `processID=0`, `nEvents=5`, `save_tree=false` |
| 入力ファイル | なし。`PHG4SimpleEventGenerator` を使用 |
| 現行生成条件 | pi0 1個、vertex `(0,0,0)`、`eta=0.5`、`phi=[-pi,pi]`、`pt=5 GeV` |
| 主な生成 node | `TOWERINFO_CALIB_CEMC`, `CLUSTERINFO_CEMC`, `CLUSTERINFO_CEMC_NO_SPLIT` |
| 出力 | `SinglePi0GunSimulation/output/DST_pi0_5GeV_eta05_towerinfo/DST_single_pi0_reconstructedInfo_<pid>.root` |
| 注意 | no-split builder の tower threshold は `0.070 GeV`。`save_tree` を使う dumper は現在コメントアウト |

実行 wrapper:

```bash
./run_job.sh <processID> <nEvents> <save_tree>
```

`Calo_Calib_ryotaro.C` は calibration と cluster setup の補助マクロで、通常は単独実行しません。

## Pi0Reconstruction

### Fun4All 実行マクロ

| マクロ | 引数 | 必要な入力 | 出力 |
| --- | --- | --- | --- |
| `Fun4All_Pi0Reconstruction.C` | `processID`, `n_events`, `cluster_node`, `input_directory`, `output_tag` | `DST_single_pi0_reconstructedInfo_<pid>.root` と指定 cluster node | `output/root/pi0_reconstruction_[<tag>_]<pid>.root` |
| `Fun4All_Pi0ReconstructionSplitComparison.C` | `processID`, `n_events`, `input_directory`, `output_directory` | 同じDST内の `CLUSTERINFO_CEMC` と `CLUSTERINFO_CEMC_NO_SPLIT` | `pi0_reconstruction_SPLIT_<pid>.root` と `...NO_SPLIT_<pid>.root` |
| `Fun4All_TowerClusterEnergyAudit.C` | `processID`, `n_events`, `tower_node`, `split_cluster_node`, `no_split_cluster_node`, `input_directory`, `output_tag` | tower node と両cluster node | `output/tower_cluster_energy_audit_[<tag>_]<pid>.root` |

`Pi0Reconstruction` 出力には `event_tree` と `h_m_gg`, `h_ncluster`, `h_cluster_e`,
`h_pair_e_asym` が入ります。vertex は現行マクロでは原点固定です。

`TowerClusterEnergyAudit` の現行 threshold は tower / cluster ともに `0.070 GeV` です。

実行 wrapper:

```bash
./run_job.sh <processID> <nEvents>
./run_split_comparison_job.sh <processID> <nEvents>
```

### 後処理マクロ

| マクロ | 入力 | 主な引数・条件 | 出力 |
| --- | --- | --- | --- |
| `MakePi0HistogramsFromTree.C` | `event_tree` | `min_cluster_energy=0.07`, mass window `[0.10,0.18] GeV`, `processID` | histogram ROOT file。入力・出力名はマクロ内で固定 |
| `CompareSplitClusterCounts.C` | split / no-split の `event_tree`。`event_uid`, `cluster_e` 必須 | optional `min_cluster_energy` | 比較 histogram、summary (`TNamed`) を含むROOTとPDF |
| `CompareSplitClusterMass.C` | split / no-split の `event_tree`。`event_uid` とcluster 4-momentum必須 | NO_SPLIT側 `ncluster==2` を起点に比較 | 比較 histogram、summaryを含むROOTと `_mass_overlay.pdf` |

`CompareSplitClusterCounts.C` と `CompareSplitClusterMass.C` は `event_uid` の重複をエラーにし、
entry 番号ではなく UID で同一イベントを結合します。

## EventDisplay

### `Fun4All_Pi0EventDisplayDump.C`

| 項目 | 内容 |
| --- | --- |
| 引数 | `processID=0`, `n_events=100`, `cluster_node="CLUSTERINFO_CEMC"`, `tower_node="TOWERINFO_CALIB_CEMC"` |
| 必要な node | `G4TruthInfo`, cluster node, tower node, `TOWERGEOM_CEMC`; hit保存時は `G4HIT_CEMC` |
| 現行入力 | 別プロジェクトの `SingleGammaGunSimulation` 出力。ファイル名はマクロ内で固定 |
| 現行出力 | `EventDisplay/output/root/0_event_display_<pid>.root` |
| 出力Tree | `events`, `truth_particles`, `truth_segments`, `cemc_clusters`, `cemc_cluster_towers`, `cemc_towers`, `cemc_hits` |
| threshold | cluster / tower ともに現行 `0.0 GeV`、hits保存は有効 |

wrapper は `./run_job.sh <processID>` です。

### 描画マクロ

| マクロ | 引数 | 入力Tree | 出力 |
| --- | --- | --- | --- |
| `draw_event_xy.C` | `filename`, `event_id=0`, `draw_other_truth=false`, `draw_hits=false` | EventDisplayの複数Tree | `output/image/event_<event>_xy.pdf` |
| `draw_event_zr.C` | 同上 | 同上 | `output/image/event_<event>_zr.pdf` |
| `draw_event_decay_zoomed.C` | `filename`, `event_id=0`, `zoom_cm=5.0` | 主にtruth Tree | `output/image/event_<event>_decay_zoomed.pdf` |

## TruthAnalysis

### `MakeTruthPi0HistogramsFromEventDisplayTree.C`

| 項目 | 内容 |
| --- | --- |
| 引数 | `processID=0`, `acceptance_eta_max=1.1`, `progress_interval=100000` |
| 入力 | EventDisplay ROOT内の `truth_particles` と `truth_segments` |
| 選別 | primary pi0 と、その direct gamma daughters。投影位置で `abs(eta1) < acceptance_eta_max` |
| 出力 | `TruthAnalysis/output/root/...root`。ファイル名はマクロ内で固定 |
| ROOT内容 | truth histogram、`truth_pi0_tree`、acceptance・event count・fractionの `TParameter` |

wrapper は `./run_job.sh <processID>` です。EventDisplay側とファイル名が一致しているか実行前に確認します。

## EmcalEtViewer

### dump マクロ

| マクロ | 引数 | 入力 | 出力 |
| --- | --- | --- | --- |
| `Fun4All_CemcTowerDumper.cc` | `process_id`, `run`, `n_events`, `save_tree=false` | マクロ内で指定したDST、`TOWERINFO_CALIB_CEMC`, `TOWERGEOM_CEMC` | `CemcTowerDumper/output/<run>/<run>-<pid>.root`、Tree `cemc_towers` |
| `Fun4All_CemcClusterDumper.cc` | 上記 + `cluster_node_name="CLUSTERINFO_CEMC"` | マクロ内で指定したDSTとcluster node | `CemcClusterDumper/output/<run>/<run>-<pid>.root`、Tree `cemc_clusters` |

両方とも `save_tree=true` の場合だけ TTree を保存します。入力ファイルと出力先は現状マクロ内に固定されています。

### 描画マクロ

| マクロ | 引数 | 入力 | 出力・動作 |
| --- | --- | --- | --- |
| `DrawCemcLego.C` | `file_name`, `tree_name="cemc_towers"`, `mode="vtx"` | `cemc_towers` | eventごとに `event.png` を更新し、Enterで次へ進む |
| `DrawCemcClusterLego.C` | `file_name`, `tree_name="cemc_clusters"`, `mode="vtx"`, `zmode="et"`, `cluster_index=-1` | `cemc_clusters` | `cluster_event.png`。cluster指定または対話表示 |

画像名が固定なので、複数イベントや複数条件を残す場合は実行後に名前を変更します。
