# 標準的な実行手順

この文書は、環境設定から結果確認までを実行順に並べたものです。
コマンドは SDCC 上の現行配置を前提とします。

## 1. 環境設定

```bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export PROJECT_ROOT=/sphenix/user/ryotaro/DirectPhotonAnalysis
```

厳密な再現が必要な場合は `ana` の代わりに実行時の固定 release を使い、
[再現メモ](reproducibility-log.md) に記録します。

## 2. Fun4All module のビルド

それぞれのモジュールをビルドします。sPHENIX のスタンダードな手法でビルドできます。

### EmcalEtViewer

`SinglePi0GunSimulation` はこの install にある `CemcTowerDumper` header を include します。  \
`build/` がまだない場合は module の直下で `mkdir -p build install` を先に実行します。

```bash
cd "${PROJECT_ROOT}/EmcalEtViewer/build"
${PROJECT_ROOT}/EmcalEtViewer/src/autogen.sh --prefix="${PROJECT_ROOT}/EmcalEtViewer/install"
make -j4
make install
source /opt/sphenix/core/bin/setup_local.sh "${PROJECT_ROOT}/EmcalEtViewer/install"
```

### Pi0Reconstruction

```bash
cd "${PROJECT_ROOT}/Pi0Reconstruction/build"
${PROJECT_ROOT}/EmcalEtViewer/src/autogen.sh --prefix="${PROJECT_ROOT}/Pi0Reconstruction/install"
make -j4
make install
source /opt/sphenix/core/bin/setup_local.sh "${PROJECT_ROOT}/Pi0Reconstruction/install"
```

### EventDisplay

```bash
cd "${PROJECT_ROOT}/EventDisplay/build"
${PROJECT_ROOT}/EmcalEtViewer/src/autogen.sh --prefix="${PROJECT_ROOT}/EventDisplay/install"
make -j4
make install
source /opt/sphenix/core/bin/setup_local.sh "${PROJECT_ROOT}/EventDisplay/install"
```

## 3. 単一 pi0 sample の生成

実行前に `Fun4All_SingleParticlePi0.C` の次の設定を確認します。

* particle、vertex、eta、phi、pt

* 出力ディレクトリ名

* CDB global tag / timestamp

* tower calibration と cluster threshold

* split / no-split node の作り方

5イベントのローカル確認例です。

```bash
cd "${PROJECT_ROOT}/SinglePi0GunSimulation/macro"
./run_job.sh 0 5 true
```

引数は順に `processID nEvents save_tree` です。現在は DST が常に保存され、
`save_tree` を使う追加 dumper はコメントアウトされています。

現行設定で期待される出力は次です。

```text
SinglePi0GunSimulation/output/DST_pi0_5GeV_eta05_towerinfo/
  DST_single_pi0_reconstructedInfo_000000.root
```

`processID=0` が `000000` に対応します。

## 4. split / no-split pi0 再構成

現行 simulation の `eta05` sample を明示的に読むローカル例です。
`n_events=0` は入力にある全イベントを処理します。

```bash
cd "${PROJECT_ROOT}/Pi0Reconstruction/macro"
source /opt/sphenix/core/bin/setup_local.sh "${PROJECT_ROOT}/Pi0Reconstruction/install"
root.exe -q -b 'Fun4All_Pi0ReconstructionSplitComparison.C(\
                  0,\
                  0,\
                  "/sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0_5GeV_eta05_towerinfo",\
                  "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/root")'
```

期待される出力は次の2ファイルです。

```text
Pi0Reconstruction/output/root/pi0_reconstruction_SPLIT_000000.root
Pi0Reconstruction/output/root/pi0_reconstruction_NO_SPLIT_000000.root
```

両方の `event_tree` で同じ `processID` とイベント順を使うため、同一イベントは同じ
`event_uid` を持ちます。詳しくは [ROOT TTree ブランチ一覧](root-tree-schema.md#event_uid) を参照してください。

一方の cluster node だけを処理する場合は、次の形です。

```bash
root.exe -q -b 'Fun4All_Pi0Reconstruction.C(\
                  0,\
                  0,\
                  "CLUSTERINFO_CEMC",\
                  "/sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0_5GeV_eta05_towerinfo",\
                  "SPLIT")'
```

## 5. split / no-split の比較

cluster 数の比較例です。以下のような、pdf 形式の画像(と画像の元になる ROOT file)が生成されます。  
![生成される画像例](../ImageSample/Pi0Reconstruction_compare_split_cluster_counts_5GeV-1.png)

最後の `0.07` は比較時の最小 cluster energy \[GeV] です。
省略した場合は energy cut を追加せず、`cluster_e` の全要素を数えます。

```bash
cd "${PROJECT_ROOT}/Pi0Reconstruction/macro"
root.exe -q -b 'CompareSplitClusterCounts.C(\
                  "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/root/pi0_reconstruction_NO_SPLIT_000000.root",\
                  "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/root/pi0_reconstruction_SPLIT_000000.root",\
                  "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/compare_split_cluster_counts_000000.root",\
                  0.07)'
```

出力は ROOT と、同じ basename の PDF です。比較は entry 番号ではなく `event_uid` で対応付けます。

cluster pair mass の比較例です。

```bash
root.exe -q -b 'CompareSplitClusterMass.C(\
                  "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/root/pi0_reconstruction_NO_SPLIT_000000.root",\
                  "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/root/pi0_reconstruction_SPLIT_000000.root",\
                  "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/compare_split_cluster_mass_000000.root")'
```

この mass 比較では、NO\_SPLIT 側で `ncluster==2` のイベントを起点にし、同じ `event_uid` の
SPLIT 側を探します。

## 6. EventDisplay Tree の作成と描画

SingleParticleGun で生成したイベントを視覚化するための event display を独自に用意しています。

||||
|---|---|---|
|![画像例1](../ImageSample/EventDisplay_event_0_decay_zoomed-1.png) | ![画像例2](../ImageSample/EventDisplay_event_0_xy_towerinforyotaro_SPLIT-1.png) | ![画像例3](../ImageSample/EventDisplay_event_0_zr_towerinforyotaro_SPLIT-1.png) |

DST file を input として受け取ります。
マクロ内の `input_file` と必要に応じて `output_file` を先に変更します。

```bash
cd "${PROJECT_ROOT}/EventDisplay/macro"
./run_job.sh 0
```

`run_job.sh` は `processID` だけを受け取り、マクロ既定値では100イベントを処理します。
現行の既定出力名は次です。

```text
EventDisplay/output/root/0_event_display_000000.root
```

XY、ZR、decay vertex 付近の描画例です。

```bash
root.exe -q -b 'draw_event_xy.C("/sphenix/user/ryotaro/DirectPhotonAnalysis/EventDisplay/output/root/0_event_display_000000.root",0,false,false)'
root.exe -q -b 'draw_event_zr.C("/sphenix/user/ryotaro/DirectPhotonAnalysis/EventDisplay/output/root/0_event_display_000000.root",0,false,false)'
root.exe -q -b 'draw_event_decay_zoomed.C("/sphenix/user/ryotaro/DirectPhotonAnalysis/EventDisplay/output/root/0_event_display_000000.root",0,5.0)'
```

PDF は `EventDisplay/output/image/` に保存されます。

## 7. truth pi0 の集計

`MakeTruthPi0HistogramsFromEventDisplayTree.C` は EventDisplay 出力の `truth_particles` と
`truth_segments` を読みます。現状では入力ファイル名が EventDisplay の現行既定出力名
`0_event_display_...` と異なるため、マクロ内の `input_file` を合わせてから実行します。

```bash
cd "${PROJECT_ROOT}/TruthAnalysis/macro"
root.exe -q -b 'MakeTruthPi0HistogramsFromEventDisplayTree.C(0,1.1,100000)'
```

引数は `processID acceptance_eta_max progress_interval` です。出力には histogram、
`truth_pi0_tree`、acceptance と集計値の `TParameter` が入ります。

