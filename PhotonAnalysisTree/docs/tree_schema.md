# TTree schema

Schema version: 1. Energy/momentumはGeV、positionはcm、angleはradianです。invalid scalar/scoreは`-999`、valid flagは`0`です。

## Index contract

- `event_tree`: 1 entry = 1 input DST event
- `event_uid = (source_file_id << 32) | event_in_file`
- SPLIT/NO_SPLIT clusterはそれぞれenergy降順、同energyならcluster ID昇順
- `<collection>_pair_cluster_i/j`は同じcollectionのcluster vector index
- `nosplit_tower_cluster_index`はNO_SPLIT cluster vector index
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

## NO_SPLIT constituent towers

| Branch | Meaning |
|---|---|
| `nosplit_ntower` | every `nosplit_tower_*` vector length |
| `nosplit_tower_cluster_index` | owner cluster index |
| `nosplit_tower_key`, `ieta`, `iphi` | RawTower key and decoded bins |
| `nosplit_tower_x/y/z/r/eta/phi` | tower geometry |
| `nosplit_tower_energy` | calibrated TowerInfo energy |
| `nosplit_tower_cluster_value` | value stored in RawCluster |
| `nosplit_tower_time/is_good/status` | TowerInfo metadata |

These raw quantities are sufficient to reproduce the ONNX inputs; derived model features are intentionally not duplicated in the base tree.

## Added score branches

| Branch | Meaning |
|---|---|
| `split_cluster_bdt_base_v3E_score` | SPLIT clusterで学習し、SPLIT clusterへ適用したphoton-ID BDT output |
| `split_cluster_bdt_base_v3E_valid` | finite input and valid shower shape |
| `nosplit_cluster_bdt_base_v3E_score` | NO_SPLIT clusterで学習し、NO_SPLIT clusterへ適用したphoton-ID BDT output |
| `nosplit_cluster_bdt_base_v3E_valid` | 上記scoreのfinite input and valid shower shape |
| `nosplit_cluster_p_gamma` | ONNX direct-gamma probability |
| `nosplit_cluster_p_gamma_valid` | exact training-time feature integrity checks passed and inference succeeded |

Scored files intentionally contain only the `event_tree` and `metadata` TTrees. Adapter counters are written to the job log. Model paths, hashes, feature order, and domain warnings are recorded once per merged dataset in `output/merged/manifest.json`.

## Metadata tree

One entry per source DST with schema version, input/output paths, node names, ordering rule, acceptance, and processed/written/invalid counters. `hadd` concatenates metadata entries while concatenating event entries.
