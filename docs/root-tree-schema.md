# ROOT TTree ブランチ一覧

この文書は `event_uid` 追加後の現行コードを対象とします。

## 共通規約

- energy と momentum: GeV
- position、radius、vertex、zoom: cm
- phi と角度: rad
- eta: 無次元
- `-999` または `-999.0`: 値が取得できない場合の sentinel（主に EventDisplay / TruthAnalysis / EmcalEtViewer）
- `is_*`, `has_*`: `0=false`, `1=true`
- vector branch は、同じグループ内の同じ添字が同じ object を表します。

Treeごとに event の数え方が異なります。`Pi0Reconstruction/event_tree` と EventDisplay は0始まり、
`CemcTowerDumper` と `CemcClusterDumper` は1始まりです。別ファイル間の対応付けに `event` だけを
使わないでください。

## `event_uid`

`Pi0Reconstruction/event_tree` の `event_uid` は次の64 bit整数です。

```cpp
event_uid = (static_cast<unsigned long long>(process_id) << 32U) | event;
```

- 上位32 bit: `process_id`
- 下位32 bit: その process 内の0始まり `event`

復元例:

```cpp
unsigned int process_id = static_cast<unsigned int>(event_uid >> 32U);
unsigned int event = static_cast<unsigned int>(event_uid & 0xffffffffULL);
```

同じ production 内で各ジョブの `process_id` が重複せず、1ジョブのイベント数が \(2^{32}\) 未満なら
UID は一意です。split / no-split の同一イベントには、同じ `process_id` と `event` を設定することで
同じ UID が入ります。

`CompareSplitClusterCounts.C` と `CompareSplitClusterMass.C` はこの UID で2本の Tree を結合し、
各入力 Tree 内の重複 UID をエラーにします。

現時点で `event_uid` を持つのは `Pi0Reconstruction/event_tree` です。EventDisplay、TruthAnalysis、
EmcalEtViewer の Tree は `event` のみなので、複数 process の file を結合する際は process 情報を
別途保持する必要があります。

## Pi0Reconstruction: `event_tree`

1 entry = 1 event です。

### Event情報

| Branch | C++型 | 説明 |
| --- | --- | --- |
| `process_id` | `unsigned int` | マクロから module に渡したprocess ID |
| `event` | `unsigned int` | process内の0始まりevent番号 |
| `event_uid` | `unsigned long long` | `process_id` と `event` を結合した64 bit UID |
| `ncluster` | `unsigned int` | `min_cluster_energy` 以上の有効cluster数 |
| `ncluster_all` | `unsigned int` | finiteなenergy / positionからphoton candidateを構成できた全cluster数 |
| `min_cluster_energy` | `double` | `ncluster` とpair作成に使ったenergy threshold [GeV] |
| `vertex_x`, `vertex_y`, `vertex_z` | `double` | momentum方向の計算に使ったvertex [cm]。現行マクロでは原点 |

### Cluster情報

以下はすべて `std::vector<double>` で、長さは `ncluster_all` です。
energy cut 前の有効clusterも入る点に注意してください。

| Branch | 説明 |
| --- | --- |
| `cluster_e` | cluster energy [GeV] |
| `cluster_x`, `cluster_y`, `cluster_z` | cluster position [cm] |
| `cluster_px`, `cluster_py`, `cluster_pz` | vertexからcluster位置へ向かうmassless photon仮定のmomentum [GeV] |

cluster vector の index は、入力 `RawClusterContainer` の cluster ID ではなく、このTree内で
有効candidateを追加した順番です。

### Cluster pair情報

pairは `min_cluster_energy` 以上のclusterの全組合せです。以下のvectorはすべて同じ長さです。

| Branch | C++型 | 説明 |
| --- | --- | --- |
| `pair_cluster_i`, `pair_cluster_j` | `std::vector<unsigned int>` | pairを構成する `cluster_*` vector のindex |
| `pair_m_gg` | `std::vector<double>` | 2 clusterの不変質量 [GeV] |
| `pair_e_asym` | `std::vector<double>` | `abs(E1-E2)/(E1+E2)`。通常 `[0,1]` |

同じ ROOT file には次の histogram も保存されます。

| Object | 内容 |
| --- | --- |
| `h_m_gg` | 選別後clusterの全pair不変質量 |
| `h_ncluster` | eventごとの `ncluster` |
| `h_cluster_e` | threshold以上のcluster energy |
| `h_pair_e_asym` | pair energy asymmetry |

## Pi0Reconstruction: `energy_audit_tree`

`TowerClusterEnergyAudit` の出力で、1 entry = 1 event です。`event` は0始まりです。

各 prefix `tower`, `split_cluster`, `no_split_cluster` に対して同じ種類の branch があります。

| Branch pattern | C++型 | 説明 |
| --- | --- | --- |
| `event` | `unsigned int` | process内event番号。UIDではない |
| `<prefix>_count` | `unsigned int` | finite energyを持つ全object数 |
| `<prefix>_count_above_threshold` | `unsigned int` | threshold以上のobject数 |
| `<prefix>_energy_sum` | `double` | finiteな全objectのenergy和 [GeV] |
| `<prefix>_energy_sum_above_threshold` | `double` | threshold以上のenergy和 [GeV] |
| `<prefix>_max_energy` | `double` | 最大energy [GeV] |
| `<prefix>_energies` | `std::vector<double>` | threshold以上の各object energy [GeV] |

実際の branch 名は次のとおりです。

| Prefix | Branch |
| --- | --- |
| `tower` | `tower_count`, `tower_count_above_threshold`, `tower_energy_sum`, `tower_energy_sum_above_threshold`, `tower_max_energy`, `tower_energies` |
| `split_cluster` | `split_cluster_count`, `split_cluster_count_above_threshold`, `split_cluster_energy_sum`, `split_cluster_energy_sum_above_threshold`, `split_cluster_max_energy`, `split_cluster_energies` |
| `no_split_cluster` | `no_split_cluster_count`, `no_split_cluster_count_above_threshold`, `no_split_cluster_energy_sum`, `no_split_cluster_energy_sum_above_threshold`, `no_split_cluster_max_energy`, `no_split_cluster_energies` |

現行 Fun4All マクロでは tower / cluster threshold ともに `0.070 GeV` です。

## EventDisplay

出力には7本の Tree が入ります。すべての `event` は EventDisplay module 内の0始まり番号で、
`event_uid` ではありません。

### `events`

1 entry = 1 event のsummaryです。

| Branch | 型 | 説明 |
| --- | --- | --- |
| `event` | `int` | event番号 |
| `n_truth_particles` | `int` | `truth_particles` に保存したparticle数 |
| `n_truth_pi0` | `int` | PID 111 のtruth particle数 |
| `n_truth_gamma` | `int` | PID 22 のtruth particle数 |
| `n_clusters` | `int` | cluster threshold以上のcluster数 |
| `n_towers_above_threshold` | `int` | tower threshold以上のtower数 |
| `n_cemc_hits` | `int` | 保存したCEMC G4 hit数 |

### `truth_particles`

1 entry = 1 truth particle です。

| Branch group | 型 | 説明 |
| --- | --- | --- |
| `event` | `int` | event番号 |
| `track_id`, `parent_id`, `primary_id`, `vtx_id` | `int` | Geant4 truth ID |
| `pid` | `int` | PDG particle ID |
| `px`, `py`, `pz`, `e`, `pt`, `p` | `double` | momentum / energy [GeV] |
| `vx`, `vy`, `vz` | `double` | production vertex [cm] |
| `eta`, `phi` | `double` | momentumから計算した方向 |
| `is_primary`, `is_pi0`, `is_gamma` | `int` | particle分類flag |
| `is_pi0_daughter` | `int` | direct parent が PID 111 なら1 |
| `ancestor_pi0` | `int` | ancestryを遡って最初に見つかったPID 111のtrack ID。なければ `-999` |
| `ancestor_gamma` | `int` | ancestryを遡って最初に見つかったPID 22のtrack ID。なければ `-999` |
| `generation` | `int` | 自身から `ancestor_gamma` まで遡った世代数。なければ `-999` |

### `truth_segments`

1 entry = 1 truth particle の表示用線分です。

| Branch group | 型 | 説明 |
| --- | --- | --- |
| `event`, `track_id`, `pid`, `parent_id` | `int` | event / particle識別子 |
| `ancestor_pi0`, `ancestor_gamma` | `int` | 対応するancestor track ID |
| `x0`, `y0`, `z0` | `double` | 線分始点（particle vertex）[cm] |
| `x1`, `y1`, `z1` | `double` | 線分終点（pi0は最初のdaughter vertex、それ以外はCEMC半径への投影）[cm] |
| `r0`, `r1` | `double` | 始点・終点の横半径 [cm] |
| `eta0`, `phi0`, `eta1`, `phi1` | `double` | 原点から見た始点・終点の方向 |
| `segment_type` | `int` | `1=pi0`, `2=pi0のdirect gamma daughter`, `3=その他` |

### `cemc_clusters`

1 entry = threshold以上の1 clusterです。

| Branch group | 型 | 説明 |
| --- | --- | --- |
| `event` | `int` | event番号 |
| `cluster_node` | `std::string` | 読み出したcluster node名 |
| `cluster_id` | `unsigned int` | `RawCluster::get_id()` |
| `energy`, `ecore` | `double` | cluster energy [GeV] |
| `chi2`, `prob`, `merged_cluster_prob` | `double` | cluster shape / probability量 |
| `x`, `y`, `z`, `r` | `double` | cluster位置 [cm] |
| `phi`, `eta_vtx0` | `double` | 原点からcluster位置を見た方向 |
| `ntowers` | `int` | cluster member tower数 |
| `lead_tower_key` | `unsigned int` | contributionが最大のmember tower key |
| `lead_tower_ieta`, `lead_tower_iphi` | `int` | lead tower index |
| `lead_tower_energy` | `double` | tower nodeから得たlead tower energy [GeV] |
| `nearest_truth_gamma_track_id` | `int` | CEMC投影位置が最も近いdirect pi0-daughter gammaのtrack ID |
| `nearest_truth_gamma_delta_eta`, `nearest_truth_gamma_delta_phi`, `nearest_truth_gamma_delta_r` | `double` | cluster − truth projection の差 |
| `nearest_truth_gamma_delta_tower` | `double` | 予約済みbranch。現行実装では計算されず `-999` |

### `cemc_cluster_towers`

1 entry = clusterとmember towerの1組です。

| Branch group | 型 | 説明 |
| --- | --- | --- |
| `event` | `int` | event番号 |
| `cluster_id` | `unsigned int` | cluster ID |
| `cluster_node` | `std::string` | cluster node名 |
| `tower_key` | `unsigned int` | member tower raw key |
| `ieta`, `iphi` | `int` | tower index |
| `cluster_tower_value` | `double` | `RawCluster` のtower mapに保存されたcontribution |
| `tower_energy` | `double` | tower nodeでのcalibrated energy [GeV] |
| `eta`, `phi` | `double` | 予約済みbranch。現行実装では `-999` |

### `cemc_towers`

1 entry = threshold以上の1 towerです。

| Branch group | 型 | 説明 |
| --- | --- | --- |
| `event`, `channel` | `int` | event番号とTowerInfo channel |
| `tower_node` | `std::string` | tower node名 |
| `key` | `unsigned int` | CEMC raw tower key |
| `ieta`, `iphi` | `int` | tower index |
| `eta`, `phi` | `double` | tower centerの方向 |
| `x`, `y`, `z` | `double` | tower center [cm] |
| `energy`, `time` | `double` | tower energy [GeV] と time |
| `is_good` | `int` | `TowerInfo::get_isGood()` |
| `status` | `int` | `TowerInfo::get_status()` |

### `cemc_hits`

`set_write_hits(true)` のとき、1 entry = 1 CEMC G4 hitです。

| Branch group | 型 | 説明 |
| --- | --- | --- |
| `event`, `trkid` | `int` | event番号とhitのtruth track ID |
| `ancestor_pi0`, `ancestor_gamma` | `int` | truth ancestry上のtrack ID |
| `x0`, `y0`, `z0`, `x1`, `y1`, `z1` | `double` | hit入口・出口 [cm] |
| `r0`, `r1` | `double` | 入口・出口の横半径 [cm] |
| `edep`, `eion`, `light_yield` | `double` | Geant4 hitのenergy deposit関連量 [GeV] |

## TruthAnalysis: `truth_pi0_tree`

1 entry = 1 primary pi0 です。入力 EventDisplay の `truth_particles` と `truth_segments` を
`(event, track_id)` で対応付けます。

| Branch group | 型 | 説明 |
| --- | --- | --- |
| `event`, `pi0_track_id` | `int` | EventDisplay event番号とprimary pi0 track ID |
| `n_direct_gamma` | `int` | direct gamma daughter数 |
| `is_pi0_to_2gamma` | `int` | direct gammaがちょうど2個なら1 |
| `gamma1_track_id`, `gamma2_track_id` | `int` | 2 gammaのtrack ID。2gammaでなければ `-999` |
| `gamma1_e`, `gamma2_e` | `double` | truth gamma energy [GeV] |
| `gamma1_px`, `gamma1_py`, `gamma1_pz` | `double` | gamma 1 momentum [GeV] |
| `gamma2_px`, `gamma2_py`, `gamma2_pz` | `double` | gamma 2 momentum [GeV] |
| `gamma1_eta1`, `gamma1_phi1`, `gamma2_eta1`, `gamma2_phi1` | `double` | EventDisplay `truth_segments` の終点方向 |
| `gamma1_in_acceptance`, `gamma2_in_acceptance` | `int` | projectionがあり `abs(eta1)<acceptance_eta_max` なら1 |
| `both_gamma_in_acceptance` | `int` | 両方acceptance内なら1 |
| `at_least_one_gamma_out_acceptance` | `int` | 両方acceptance内でなければ1 |
| `missing_gamma_projection` | `int` | 少なくとも一方のprojectionがなければ1 |
| `m_gg` | `double` | truth 2 gamma invariant mass [GeV] |
| `pair_e_asym` | `double` | `abs(E1-E2)/(E1+E2)` |

2gammaでない entry では、gamma固有量、`m_gg`、`pair_e_asym` は `-999` です。

## EmcalEtViewer: `cemc_towers`

1 entry = 1 event、各vectorの1要素 = 1 towerです。`event` は1始まりです。

| Branch | 型 | 説明 |
| --- | --- | --- |
| `event` | `int` | module内event番号 |
| `vtx_x`, `vtx_y`, `vtx_z` | `float` | `GlobalVertexMap`から選んだvertex [cm]。取得失敗時 `-999` |
| `eta0`, `phi0` | `std::vector<float>` | 原点からtower centerを見た方向 |
| `etavtx`, `phivtx` | `std::vector<float>` | event vertexからtower centerを見た方向 |
| `et0` | `std::vector<float>` | `energy/cosh(eta0)` [GeV] |
| `etvtx` | `std::vector<float>` | `energy/cosh(etavtx)` [GeV] |
| `energy` | `std::vector<float>` | tower energy [GeV] |
| `ieta`, `iphi` | `std::vector<int>` | geometry上のtower index |

全tower vectorは同じ長さ・同じ添字対応です。

## EmcalEtViewer: `cemc_clusters`

1 entry = 1 eventです。`event` は1始まりです。`cluster_*` vectorの1要素が1 cluster、
`member_*` vectorの1要素が1 cluster-member tower対応を表します。

### Event / cluster vector

| Branch group | 型 | 説明 |
| --- | --- | --- |
| `event`, `has_vertex` | `int` | event番号とvertex取得flag |
| `vtx_x`, `vtx_y`, `vtx_z` | `float` | event vertex [cm]。取得失敗時 `-999` |
| `cluster_id` | `std::vector<unsigned int>` | `RawCluster` ID |
| `cluster_energy`, `cluster_ecore` | `std::vector<float>` | energy [GeV] |
| `cluster_chi2`, `cluster_prob`, `cluster_merged_prob` | `std::vector<float>` | cluster probability関連量 |
| `cluster_et_iso`, `cluster_mean_time` | `std::vector<float>` | isolation ET と mean time |
| `cluster_x`, `cluster_y`, `cluster_z`, `cluster_r` | `std::vector<float>` | cluster位置 [cm] |
| `cluster_phi_det` | `std::vector<float>` | `RawCluster::get_phi()` |
| `cluster_eta0`, `cluster_phi0`, `cluster_et0` | `std::vector<float>` | 原点基準の方向とtransverse energy |
| `cluster_etavtx`, `cluster_phivtx`, `cluster_etvtx` | `std::vector<float>` | event vertex基準。vertexなしでは `-999` |
| `cluster_n_towers` | `std::vector<int>` | member tower数 |
| `cluster_lead_tower_ieta`, `cluster_lead_tower_iphi` | `std::vector<int>` | lead tower index |
| `cluster_x_tower_raw`, `cluster_y_tower_raw`, `cluster_x_tower_corr`, `cluster_y_tower_corr` | `std::vector<float>` | tower座標系でのraw / corrected cluster位置 |

### Member tower vector

| Branch | 型 | 説明 |
| --- | --- | --- |
| `member_cluster_index` | `std::vector<int>` | `cluster_*` vectorへのindex |
| `member_cluster_id` | `std::vector<unsigned int>` | 親cluster ID |
| `member_tower_key` | `std::vector<unsigned int>` | tower key |
| `member_tower_caloid` | `std::vector<int>` | keyからdecodeしたcalorimeter ID |
| `member_tower_ieta`, `member_tower_iphi` | `std::vector<int>` | keyからdecodeしたtower index |
| `member_energy` | `std::vector<float>` | cluster tower mapのmember energy / contribution |
| `member_energy_fraction` | `std::vector<float>` | `member_energy / cluster_energy` |

## 簡単な確認例

```cpp
_file0->ls();
event_tree->Print();
event_tree->Scan("process_id:event:event_uid:ncluster:ncluster_all", "", "", 10);
event_tree->Draw("pair_m_gg");
```

vectorの対応を確認する例:

```cpp
event_tree->Scan("event:cluster_e.size():cluster_px.size():pair_m_gg.size()", "", "", 10);
```
