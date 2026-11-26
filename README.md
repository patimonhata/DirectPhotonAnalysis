# DirectPhotonAnalysis

sPHENIX の photon / pi0 解析、入力シミュレーション、truth matching と event display をまとめたプロジェクトです。
現在の主解析は [PhotonAnalysisTree の photon candidate selection](PhotonAnalysisTree/workflows/photon_candidate_selection/README.md) です。

## 目的から探す

| 目的 | 入口 |
|---|---|
| 主解析: photon candidate の map → reduce → merge | [PhotonAnalysisTree](PhotonAnalysisTree/README.md) |
| Pythia / single-particle の解析 tree 生成と score 追加 | [workflow 一覧](PhotonAnalysisTree/workflows/README.md) |
| 独立解析・過去 study の再実行 | [studies](PhotonAnalysisTree/workflows/studies/README.md) |
| 主解析の分類・tagging を調べる補助解析 | [diagnostics](PhotonAnalysisTree/workflows/diagnostics/README.md) |
| 生成ファイルの整合性確認・小規模 QA | [validation](PhotonAnalysisTree/validation/README.md) |
| topology 判定を event display で見る | [TopologyEventDisplay](TopologyEventDisplay/README.md) |
| CEMC / HCAL の TopoCluster 比較 | [TopoClusterHCalStudy](TopoClusterHCalStudy/README.md) |
| 単一 pi0 の simulation・再構成・表示 | 以下の解析フローと [実行手順](docs/workflow.md) |

PhotonAnalysisTree のコードは、実行入口を `workflows/`、ROOT 集計・描画の内部実装を `analysis/`、
Fun4All module と共有の物理処理を `src/`、出力検証を `validation/`、自動テストを `tests/` に分けています。
詳細は [配置方針](PhotonAnalysisTree/docs/layout.md) を参照してください。

## 単一粒子シミュレーションからの解析の流れ

```text
SinglePi0GunSimulation
  └─ DST_single_pi0_reconstructedInfo_<processID>.root
       ├─ Pi0Reconstruction
       │    ├─ pi0_reconstruction_SPLIT_<processID>.root
       │    └─ pi0_reconstruction_NO_SPLIT_<processID>.root
       │          └─ split / no-split 比較、再ヒストグラム化
       ├─ EventDisplay
       │    └─ events, truth_particles, truth_segments,
       │       cemc_clusters, cemc_cluster_towers, cemc_towers, cemc_hits
       │          ├─ XY / ZR / decay zoomed PDF
       │          └─ TruthAnalysis
       └─ EmcalEtViewer
            └─ tower / cluster lego plot
```

現在の各マクロの既定入力が、常にこの順序でそのまま接続されているとは限りません。
特に energy、eta、入力ディレクトリ、cluster node、出力タグを実行前に確認してください。

## ドキュメント

* [標準的な実行手順](docs/workflow.md): 環境設定、ビルド、ローカル確認、Condor 投入の順序

* [マクロ入出力一覧](docs/macro-reference.md): 各マクロの引数、必要な node / Tree、生成ファイル

* [ROOT TTree ブランチ一覧](docs/root-tree-schema.md): Tree スキーマと branch の意味

過去の生成条件や解析手順は、対象のコミットで確認できます（例: `git show <commit>:SinglePi0GunSimulation/macro/Fun4All_SingleParticlePi0.C`）。

## ディレクトリ

| ディレクトリ                    | 役割                                              |
| ------------------------- | ----------------------------------------------- |
| `SinglePi0GunSimulation/` | 単一pi0の生成、検出器シミュレーション、DST 出力                     |
| `Pi0Reconstruction/`      | cluster pair によるpi0再構成と split / no-split 比較     |
| `EventDisplay/`           | truth、CEMC hit / tower / cluster を表示用 TTree に保存 |
| `TruthAnalysis/`          | event-display Tree から primary pi0→2gamma を集計    |
| `EmcalEtViewer/`          | CEMC tower / cluster の lego plot 用 dump と描画     |
| `PhotonAnalysisTree/`     | photon candidate 主解析、tree 生成、truth/topology の共有処理 |
| `TopologyEventDisplay/`   | 共有 topology evaluator を使った診断用 event display |
| `TopoClusterHCalStudy/`   | all-calorimeter TopoCluster の CEMC / HCAL 比較 |

従来の package では `src/` は Fun4All module、`macro/` は実行・描画マクロ、`build/` と `install/` はローカルビルド、
`output/` は生成物です。

## 重要な規約

* 現行パスは `/sphenix/user/ryotaro/DirectPhotonAnalysis` です。

* `processID` はファイル名では `000000` のような6桁表記になります。

* `CLUSTERINFO_CEMC` は通常の split cluster、`CLUSTERINFO_CEMC_NO_SPLIT` は subcluster splitting を無効にして追加した node です。

* ROOT 出力の多くは `RECREATE` で開かれるため、同名ファイルを再実行すると上書きされます。

* `Pi0Reconstruction/event_tree` のイベント対応付けには `event` 単独ではなく `event_uid` を使います。

* energy / momentum の単位は GeV、位置は cm、角度は rad です。

## 単一粒子シミュレーションの動作確認

まず [標準的な実行手順](docs/workflow.md) の環境設定とビルドを行い、5イベントだけ生成します。

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/macro
./run_job.sh 0 5 true
```

以下は、単一 pi0 を `pt=5 GeV`、`eta=0.5`、`phi=[-pi,pi]` で生成する場合の出力例です。実際の保存先は生成マクロの `outDir` で確認してください。

```text
SinglePi0GunSimulation/output/DST_pi0_5GeV_eta05_towerinfo/
  DST_single_pi0_reconstructedInfo_000000.root
```

後段へ進む前に、入力ディレクトリと node 名がこの DST の設定に合っていることを確認してください。

## 単一粒子フローの注意点

* `Pi0Reconstruction` の入力ディレクトリには、使用する DST を生成したディレクトリを指定します。

* `EventDisplay/macro/Fun4All_Pi0EventDisplayDump.C` の現行既定入力は別プロジェクトの single-gamma sample です。

* `TruthAnalysis` を EventDisplay の後に使う場合は、両マクロ内の入力・出力ファイル名が一致しているか確認します。

* `source /opt/sphenix/core/bin/sphenix_setup.sh -n ana` は rolling release です。厳密に再現したい解析では固定 release を指定します。
