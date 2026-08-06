# TopoCluster HCAL study

高運動量のsingle gammaとpi0 decay photonについて、all-calorimeter TopoCluster内の
CEMC energyとHCAL energyの相関を比較する解析です。入力DSTに保存済みの
`TOPOCLUSTER_ALLCALO`を読み、TopoClusterの再構成や追加matchingは行いません。

## 入出力

入力は各500 filesです。

- `SinglePi0GunSimulation/output/DST_gamma_25to35GeV_etapm1_vertexpm60`
- `SinglePi0GunSimulation/output/DST_pi0_25to35GeV_etapm1_vertexpm60`

出力は1 DSTにつき1 ROOT fileです。

- `output/root/topocluster_hcal_gamma_000000.root` ... `000499.root`
- `output/root/topocluster_hcal_pi0_000000.root` ... `000499.root`
- `output/merge/`: `hadd`後のfile

## Tree schema

Tree名は`topocluster_tree`、1 entryは1 eventです。energyの単位はGeVです。
同じvector index `i`は、すべて同じ`TOPOCLUSTER_ALLCALO` clusterを表します。

| branch | 型 | 内容 |
| --- | --- | --- |
| `sample_id` | `unsigned int` | `0=gamma`, `1=pi0` |
| `job_index` | `unsigned int` | Condor job index。gammaは0--499、pi0は500--999 |
| `process_id` | `unsigned int` | sample内のDST番号、0--499 |
| `event` | `unsigned int` | DST内の0始まりevent番号 |
| `event_uid` | `unsigned long long` | `(job_index << 32) | event` |
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

`other_calo_energy`と`energy_residual`はenergy分解の検証用で、通常は0または浮動小数点精度程度です。

## Build

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy
./src/build.sh
```

`build.sh`は`src/autogen.sh`、`make`、`make install`を順に実行し、headerを`install/include/topoclusterhcalstudy/`、libraryを`install/lib/`へinstallします。

既定releaseは既存ジョブと同じ`ana`です。固定releaseに変更する場合は、buildと実行の両方で同じ値を指定します。

```bash
SPHENIX_RELEASE=ana.557 ./src/build.sh
```

## Local test

`JOB_INDEX=0`はgamma file 000000、`JOB_INDEX=500`はpi0 file 000000です。
第2引数は処理event数で、0ならDST末尾まで処理します。

```bash
./run_job.sh 0 2
./run_job.sh 500 2
```

## Condor

```bash
condor_submit run_condor.job
```

1000 jobsを投入します。既存の出力fileは上書きしないため、再投入前に対象jobの出力状況を確認してください。

## Merge

1000 filesが揃った後に実行します。

```bash
./merge_outputs.sh
```

gamma、pi0、および両sampleを合わせた3 filesを`output/merge`に作ります。

## 二次元相関の例

ROOTではvector branch同士の同じindexが対応するため、次のように全TopoClusterを描けます。

```cpp
topocluster_tree->Draw("hcal_total_energy:emcal_energy", "", "colz");
```
