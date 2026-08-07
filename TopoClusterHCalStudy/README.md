# TopoCluster HCAL study

高運動量のsingle gammaとpi0 decay photonについて、all-calorimeter TopoCluster内の
CEMC energyとHCAL energyの相関を比較する解析です。入力DSTに保存済みの
`TOPOCLUSTER_ALLCALO`を読み、TopoClusterの再構成や追加matchingは行いません。

## 入出力

Fun4All macroはsample名と入力directoryを引数で受け取ります。
`process_id=12`, `sample=gamma`の場合は、入力directory直下の次のfileを読みます。

```text
DST_single_gamma_reconstructedInfo_000012.root
```

pi0の場合は`DST_single_pi0_reconstructedInfo_000012.root`です。
出力は1 DSTにつき1 ROOT fileです。

```text
output/root/topocluster_hcal_gamma_000012.root
output/root/topocluster_hcal_pi0_000012.root
output/merge/
```

## Tree schema

Tree名は`topocluster_tree`、1 entryは1 eventです。energyの単位はGeVです。
同じvector index `i`は、すべて同じ`TOPOCLUSTER_ALLCALO` clusterを表します。

| branch | 型 | 内容 |
| --- | --- | --- |
| `sample_id` | `unsigned int` | `0=gamma`, `1=pi0` |
| `job_index` | `unsigned int` | sample内のCondor `process_shift` |
| `process_id` | `unsigned int` | 入力DSTの6桁file番号 |
| `event` | `unsigned int` | DST内の0始まりevent番号 |
| `event_uid` | `unsigned long long` | sample、process ID、eventから作るUID |
| `n_topocluster` | `unsigned int` | event内のall-calorimeter TopoCluster数 |
| `topocluster_id` | `vector<unsigned int>` | `RawCluster::get_id()` |
| `emcal_energy` | `vector<float>` | CEMC tower contributionの和 |
| `hcalin_energy` | `vector<float>` | inner HCAL tower contributionの和 |
| `hcalout_energy` | `vector<float>` | outer HCAL tower contributionの和 |
| `hcal_total_energy` | `vector<float>` | `hcalin_energy + hcalout_energy` |
| `topocluster_energy` | `vector<float>` | `RawCluster::get_energy()` |
| `other_calo_energy` | `vector<float>` | 想定外のcalorimeter IDを持つcontributionの和 |
| `energy_residual` | `vector<float>` | `topocluster_energy`から全component和を引いた値 |
| `<calo>_ntower` | `vector<unsigned int>` | 各componentを構成するtower数 |

`event_uid`の上位32 bitには次の`source_id`を使い、gammaとpi0で同じ
process IDを使っても衝突しないようにします。

```cpp
source_id = (sample_id << 31U) | process_id;
event_uid = (static_cast<unsigned long long>(source_id) << 32U) | event;
```

## Build

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy
./src/build.sh
```

`build.sh`は`src/autogen.sh`、`make`、`make install`を順に実行し、headerを
`install/include/topoclusterhcalstudy/`、libraryを`install/lib/`へinstallします。
既定releaseは既存ジョブと同じ`ana`です。

## Local test

`run_job.sh`の引数は順に`PROCESS_ID SAMPLE INPUT_DIRECTORY [N_EVENTS]`です。

```bash
./run_job.sh 0 gamma \
  /sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_gamma_25to35GeV_etapm1_vertexpm60 2

./run_job.sh 0 pi0 \
  /sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0_25to35GeV_etapm1_vertexpm60 2
```

## Condor

`run_condor.job`の次の値を、1回のsubmitで処理するsampleに合わせます。

```condor
sample          = gamma
input_directory = /path/to/gamma/DST/directory
job_offset      = 0
Queue 500
```

`Queue N`で作られる`Process=0 ... N-1`に`job_offset`を加え、
`process_shift = Process + job_offset`を入力file番号とします。したがって処理範囲は
`job_offset ... job_offset + N - 1`です。

gammaとpi0は別々に設定してsubmitします。

```bash
condor_submit run_condor.job
```

## Merge

gammaとpi0の必要なjobが揃った後に実行します。file数は固定していません。

```bash
./merge_outputs.sh
```

gamma、pi0、および両sampleを合わせた3 filesを`output/merge`に作ります。

## 二次元相関の例

ROOTではvector branch同士の同じindexが対応するため、次のように全TopoClusterを描けます。

```cpp
topocluster_tree->Draw("hcal_total_energy:emcal_energy", "", "colz");
```
