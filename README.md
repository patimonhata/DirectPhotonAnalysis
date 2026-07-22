# DirectPhotonAnalysis

単一粒子シミュレーションを入力として、CEMC の tower / cluster、pi0→\$\gamma\$ gamma 再構成、
event display、truth 情報を調べるための個人用解析コードです。

この README は「しばらく離れた後に、どこから再開すればよいか」を短時間で思い出すための入口です。
詳細は `docs/` に分けています。

## 解析の流れ

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

* [ROOT TTree ブランチ一覧](docs/root-tree-schema.md): `event_uid` 追加後の Tree スキーマと branch の意味

* [再現メモ](docs/reproducibility-log.md): 解析条件と生成物を後から追うための短い記録テンプレート

## ディレクトリ

| ディレクトリ                    | 役割                                              |
| ------------------------- | ----------------------------------------------- |
| `SinglePi0GunSimulation/` | 単一pi0の生成、検出器シミュレーション、DST 出力                     |
| `Pi0Reconstruction/`      | cluster pair によるpi0再構成と split / no-split 比較     |
| `EventDisplay/`           | truth、CEMC hit / tower / cluster を表示用 TTree に保存 |
| `TruthAnalysis/`          | event-display Tree から primary pi0→2gamma を集計    |
| `EmcalEtViewer/`          | CEMC tower / cluster の lego plot 用 dump と描画     |

`src/` は Fun4All module、`macro/` は実行・描画マクロ、`build/` と `install/` はローカルビルド、
`output/` は生成物です。

## 重要な規約

* 現行パスは `/sphenix/user/ryotaro/DirectPhotonAnalysis` です。

* `processID` はファイル名では `000000` のような6桁表記になります。

* `CLUSTERINFO_CEMC` は通常の split cluster、`CLUSTERINFO_CEMC_NO_SPLIT` は subcluster splitting を無効にして追加した node です。

* ROOT 出力の多くは `RECREATE` で開かれるため、同名ファイルを再実行すると上書きされます。

* `Pi0Reconstruction/event_tree` のイベント対応付けには `event` 単独ではなく `event_uid` を使います。

* energy / momentum の単位は GeV、位置は cm、角度は rad です。

## 最短の動作確認

まず [標準的な実行手順](docs/workflow.md) の環境設定とビルドを行い、5イベントだけ生成します。

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/macro
./run_job.sh 0 5 true
```

現在の generator 設定では、主な出力は次です。

```text
SinglePi0GunSimulation/output/DST_pi0_5GeV_eta05_towerinfo/
  DST_single_pi0_reconstructedInfo_000000.root
```

後段へ進む前に、入力ディレクトリと node 名がこの DST の設定に合っていることを確認してください。

## 現在の注意点

* `SinglePi0GunSimulation` の現行 generator は `pt=5 GeV`, `eta=0.5`, `phi=[-pi,pi]` です。

* `Pi0Reconstruction` の既定入力ディレクトリは `DST_pi0_5GeV_eta0_towerinfo` なので、現行 simulation 出力 `eta05` を読む場合は引数で入力ディレクトリを指定します。

* `EventDisplay/macro/Fun4All_Pi0EventDisplayDump.C` の現行既定入力は別プロジェクトの single-gamma sample です。

* `TruthAnalysis` を EventDisplay の後に使う場合は、両マクロ内の入力・出力ファイル名が一致しているか確認します。

* `source /opt/sphenix/core/bin/sphenix_setup.sh -n ana` は rolling release です。厳密に再現したい解析では実際に使った release 名を再現メモに残します。
