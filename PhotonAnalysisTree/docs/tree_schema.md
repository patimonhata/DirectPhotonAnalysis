# TTree schema

Schema version: 4. Energy/momentumはGeV、positionはcm、angleはradianです。invalid scalar/scoreは`-999`、valid flagは`0`です。

## Index contract

- `event_tree`: 1 entry = 1 input DST event
- `event_uid = (source_file_id << 32) | event_in_file`
- SPLIT/NO_SPLIT clusterはそれぞれenergy降順、同energyならcluster ID昇順
- `<collection>_pair_cluster_i/j`は同じcollectionのcluster vector index
- `<collection>_tower_cluster_index`は同じcollectionのcluster vector index
- score vectorは対応するcluster vectorと同じ長さ・順序
- shower patchはclusterごとに7x7 = 49要素をcluster順にflatten

## Common event and truth

| Branch group | Meaning |
|---|---|
| `source_file_id`, `event_in_file`, `event_uid` | file/event identity |
| `vertex_x/y/z` | reconstruction vertex; current production is origin |
| `label` | gamma=1, pi0=0, invalid=-999 |
| `truth_valid` | expected single primary truth was decoded successfully |
| `truth_primary_pdg_id`, `truth_primary_track_id` | primary identity |
| `truth_e/px/py/pz/pt/eta/phi` | primary kinematics |
| `truth_vx/vy/vz` | primary production vertex |
| `truth_n_direct_daughter`, `truth_is_pi0_to_2gamma` | decay summary |
| `truth_daughter_*` | direct daughters, track ID ascending |
| `truth_daughter_projection_eta/phi/valid` | straight projection from daughter vertex to CEMC tower radius |
| `truth_daughter_in_acceptance` | valid projection and `abs(eta) < acceptance_eta_max` |
| `truth_both_gamma_in_acceptance` | both direct gamma daughters accepted |
| `truth_at_least_one_gamma_out_acceptance` | inverse of previous flag for a two-gamma decay |
| `truth_missing_gamma_projection` | at least one direct gamma projection failed |
| `truth_m_gg`, `truth_pair_e_asym` | truth two-gamma mass and energy asymmetry |

Truth projection and acceptance reproduce the definitions used by`TruthAnalysis/MakeTruthPi0HistogramsFromEventDisplayTree.C`; direct daughters are stored once as vectors instead of duplicated `gamma1_*`/`gamma2_*` scalars.

`require_truth_node=false`で実行し、入力DSTにtruth nodeがない場合もeventは保存されます。その場合は`truth_valid=0`、truth scalarはinvalid値、truth vectorは空です。

## Cluster collections

`<c>` is `split` or `nosplit`.

| Branch | Type/meaning |
|---|---|
| `<c>_ncluster` | `UInt_t`, length of every ordinary `<c>_cluster_*` vector |
| `<c>_cluster_id` | original RawCluster ID |
| `<c>_cluster_ntower` | number of cluster tower entries |
| `<c>_cluster_e/et/eta/phi` | reconstructed kinematics from origin |
| `<c>_cluster_x/y/z`, `<c>_cluster_px/py/pz` | position and massless momentum |
| `<c>_cluster_shower_*` | `ShowerShapeCalculator` outputs used by the BDT |
| `<c>_cluster_shower_patch_*` | flattened 7x7 patch; length `49*ncluster` |
| `<c>_pair_cluster_i/j` | pair endpoint indices |
| `<c>_pair_m_gg`, `<c>_pair_e_asym` | reconstructed pair quantities |

SPLITは`CLUSTERINFO_CEMC`、NO_SPLITは`CLUSTERINFO_CEMC_NO_SPLIT`です。同名のkinematic branchを共通化しないのは、両collectionが一般に異なるcluster個数・ID・energyを持つためです。

## Single-pi0 energy-contribution matching

`<c>`は`split`または`nosplit`です。各vectorは`<c>_cluster_*`と同じ長さ・順序です。

| Branch | Meaning |
|---|---|
| `<c>_cluster_truth_match_valid` | constituent towerからcell/hit provenanceを完全に取得できた |
| `<c>_cluster_truth_total_edep` | cluster constituent towersに属する全G4 hit energy deposit |
| `<c>_cluster_truth_gamma0/1_edep` | direct daughter gamma 0/1の子孫trackによるdeposit |
| `<c>_cluster_truth_other_edep` | direct daughter gamma 0/1へ遡れないdeposit |
| `<c>_cluster_truth_gamma0/1_fraction` | gamma deposit / total deposit |
| `<c>_cluster_truth_other_fraction` | other deposit / total deposit |
| `<c>_cluster_truth_gamma0/1_recovery` | gamma deposit / direct daughter truth energy |

`TOWER_CALIB_CEMC -> G4CELL_CEMC -> G4HIT_CEMC`を辿り、各hitの
`track_id`からG4 ancestryを遡って`truth_daughter_track_id[0/1]`へ最初に
到達したgammaへ、cellのhit mapに記録されたdepositを割り当てます。
single-pi0 gunでは2本のgammaが同じprimary showerを共有するため、
shower IDだけでは両者を分離できません。

SPLIT collectionではPythia matcherと同様に、各tower depositへ
`clamp(cluster tower value / calibrated TowerInfo energy, 0, 1)`を掛けます。
NO_SPLIT collectionではtower deposit全量を使います。既存DSTに
`TOWER_CALIB_CEMC`、`G4CELL_CEMC`、`G4HIT_CEMC`が残っていればsimulation
DSTの再生成は不要ですが、このbranchを追加したPhotonAnalysisTreeは再生成が必要です。
いずれかの対応付けが欠けるclusterでは`truth_match_valid=0`になります。

`require_nosplit_cluster_node=false`で実行し、入力DSTにNO_SPLIT nodeがない場合は、`nosplit_ncluster=0`、`nosplit_ntower=0`、すべての`nosplit_*` vectorが空になります。SPLIT collection、calibrated tower、tower geometryは常に必須です。

## Constituent towers

`<c>`は`split`または`nosplit`です。

| Branch | Meaning |
|---|---|
| `<c>_ntower` | every `<c>_tower_*` vector length; SPLITではclusterごとのtower entry数の和 |
| `<c>_tower_cluster_index` | owner cluster index |
| `<c>_tower_key`, `ieta`, `iphi` | RawTower key and decoded bins |
| `<c>_tower_x/y/z/r/eta/phi` | tower geometry |
| `<c>_tower_energy` | 未分配のcalibrated TowerInfo energy |
| `<c>_tower_cluster_value` | RawClusterに保存された当該clusterへの割当energy |
| `<c>_tower_time/is_good/status` | TowerInfo metadata |

NO_SPLIT collectionではtower keyはevent内で一意です。SPLIT collectionでは同じtower keyが複数clusterに現れ得るため、`split_ntower`はunique tower数ではなくentry数です。SPLIT ONNXの`log1p(tower_energy)`とenergy fractionには、cluster energyとの整合性を保つため`split_tower_cluster_value`を使います。これらのraw quantityからONNX inputを再現し、derived featureはbase treeへ重複保存しません。

## Added score branches

| Branch | Meaning |
|---|---|
| `split_cluster_bdt_base_v3E_score` | SPLIT clusterで学習し、SPLIT clusterへ適用したphoton-ID BDT output |
| `split_cluster_bdt_base_v3E_valid` | finite input and valid shower shape |
| `split_cluster_bdt_ppg15v1_score` | low-pT sampleを加えて学習したPPG15v1 photon-ID BDT output |
| `split_cluster_bdt_ppg15v1_valid` | 上記scoreのfinite input and valid shower shape |
| `nosplit_cluster_bdt_base_v3E_score` | NO_SPLIT clusterで学習し、NO_SPLIT clusterへ適用したphoton-ID BDT output |
| `nosplit_cluster_bdt_base_v3E_valid` | 上記scoreのfinite input and valid shower shape |
| `nosplit_cluster_p_gamma` | ONNX direct-gamma probability |
| `nosplit_cluster_p_gamma_valid` | exact training-time feature integrity checks passed and inference succeeded |
| `split_cluster_p_gamma` | SPLIT学習ONNXによるSPLIT clusterのdirect-gamma probability |
| `split_cluster_p_gamma_valid` | SPLIT constituent feature integrity checks passed and inference succeeded |

Scored files intentionally contain only the `event_tree` and `metadata` TTrees. Adapter counters are written to the job log. Model paths, hashes, feature order, and domain warnings are recorded once per merged dataset in `output/merged/manifest.json`.

## Metadata tree

One entry per source DST with schema version, input/output paths, node names, required-node flags, ordering rule, acceptance, and processed/written/invalid counters. `hadd` concatenates metadata entries while concatenating event entries.
