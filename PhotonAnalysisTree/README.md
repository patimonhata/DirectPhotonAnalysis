# PhotonAnalysisTree

Single-particle gun DSTを1回だけ読み、truth、SPLIT cluster、NO_SPLIT cluster、再構成gamma pairを同じ`event_tree`へ保存するFun4All subsystemです。`event_tree`は1 entry = 1 input eventです。

## なぜ独立ディレクトリか

`Pi0Reconstruction`はcluster collectionからpairを再構成する部品です。このpackageはtruth、2種類のcluster collection、ML用tower point set、2つの推論器を束ねる上位のI/O pipelineなので、`DirectPhotonAnalysis/PhotonAnalysisTree`として分離しています。

重複しそうなcluster branchは削除せず、物理的に異なるcollectionとして`split_`と`nosplit_` prefixで分けます。event ID、vertex、truthは共通なので1組だけです。両collectionは必ずcluster energy降順、同値ならcluster ID昇順です。全cluster vector、pair index、score vectorはこの順序を共有します。

## 内容

- truth: primary、direct daughter、CEMC半径へのgamma投影、acceptance、truth mass/asymmetry
- `split_*`: `CLUSTERINFO_CEMC`、shower shape、全cluster pair
- `nosplit_*`: `CLUSTERINFO_CEMC_NO_SPLIT`、shower shape、全cluster pair、cluster constituent tower
- `split_cluster_bdt_base_v3E_*`: PPG15/PPG12 `base_v3E` BDT adapterが後段で追加
- `nosplit_cluster_p_gamma*`: 同梱ONNX modelが後段で追加

完全なbranch契約は[docs/tree_schema.md](docs/tree_schema.md)を参照してください。

## Build

先に既存`Pi0Reconstruction`をbuild/installしてから実行します。

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis
./PhotonAnalysisTree/build.sh
```

buildは次を作ります。

- `libPhotonAnalysisTree.so`: DSTから統合TTreeを作るFun4All subsystem
- `libPi0GammaOnnx.so`: ROOT/Fun4Allに依存しない単一cluster ONNX interface

## DSTからTTreeを作る

```bash
./PhotonAnalysisTree/run_tree.sh 0 0
```

第1引数はsource/process ID、第2引数はevent数です。`0` event指定は全eventです。default inputは依頼時の`newDST_pi0_5to15GeV_etapm1`、default outputは`PhotonAnalysisTree/output/root`です。

直接ROOT macroを呼ぶ場合は、input/output directoryとexpected primary PDGも変更できます。

## scoreを追加する

inputは変更せず、BDT追加fileと最終fileを別々に作ります。

```bash
./PhotonAnalysisTree/run_add_scores.sh \
  PhotonAnalysisTree/output/root/photon_analysis_tree_000000.root \
  PhotonAnalysisTree/output/root/photon_analysis_tree_000000_with_bdt.root \
  PhotonAnalysisTree/output/root/photon_analysis_tree_000000_scored.root
```

wrapperは最後に`check_scored_tree.C`を実行し、cluster、tower、pair、score vectorの長さに加え、`metadata` TTreeが実際に読めること、event数とsource file IDが一致することを検証します。scored ROOT fileのtop-level objectは`event_tree`と`metadata`だけです。adapterのcluster数、valid数、malformed数などの集計は標準出力（Condor log）だけに記録します。既存outputは上書きしません。Condorで多数jobを同時実行してもACLiCの共有build fileが競合しないよう、score macroは各process内でloadします。

## scored fileをmergeする

merged fileは個別fileと同じdirectoryへ置かず、入力に6桁のprocess IDを明示したpatternを使います。

```bash
mkdir -p PhotonAnalysisTree/output/merged
hadd -f PhotonAnalysisTree/output/merged/all.root \
  PhotonAnalysisTree/output/root/photon_analysis_tree_[0-9][0-9][0-9][0-9][0-9][0-9]_scored.root
```

正常なmergeでは`event_tree`と`metadata`がそれぞれmergeされ、`metadata`のentry数は入力file数になります。個別fileに`TNamed`や`TParameter`を保存しないため、同名objectのbackup cycleは発生しません。

merge後にdataset全体の生成条件を`manifest.json`へ記録します。scriptはROOTのkey構成、event/metadata数、model SHA-256を検証してからmanifestを作ります。

```bash
./PhotonAnalysisTree/make_dataset_manifest.sh
```

default出力は`PhotonAnalysisTree/output/merged/manifest.json`です。引数は`SCORED_DIR MERGED_ROOT MANIFEST_JSON BDT_MODEL ONNX_MODEL`の順で、省略できます。既存manifestを意図的に更新する場合だけ`FORCE=1`を指定します。

## Model provenanceと制約

model path、SHA-256、score branch、feature schema、domain warning、event数、source file範囲、codeのgit commitはdataset単位の`output/merged/manifest.json`へ1回だけ保存します。個別ROOT fileには重複保存しません。

同梱`models/best_model.onnx`のSHA-256は`194eb9a0394b5752619d324aa8a17dde547ff67713ea908dcfed869ef0931a6d`です。元は`Pi0DirectGammaSeparation/Pi0GammaClassifier/output/baseline`で、`models/onnx_metadata.json`も同梱しています。

ONNX modelは`n_cluster == 1`のNO_SPLIT eventだけで学習されています。adapterはmulti-cluster eventもclusterごとに評価しますが、そのscoreは学習条件外なので物理解析で使う前に別途validationが必要です。

BDT modelのdocumented performance binはcluster ET 6 GeVからです。それより低いETのscoreはsoftware上は計算されますが、性能保証範囲外です。

## HTCondor production

実行scriptとjob fileは`PhotonAnalysisTree`直下にあります。

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree
condor_submit run_tree.job
# tree jobsの完了後
condor_submit run_add_scores.job
```

job数やoffsetを変更する場合はsubmit時に上書きできます。

```bash
condor_submit -append "n_jobs = 100" -append "job_offset = 500" run_tree.job
condor_submit -append "n_jobs = 100" -append "job_offset = 500" run_add_scores.job
```

- `run_tree.sh`, `run_tree.job`: DSTからbase TTreeを生成
- `run_add_scores.sh`, `run_add_scores.job`: BDT/gamma score branchを追加して検証
- `make_dataset_manifest.sh`: `hadd`後のROOT fileとmodelを検証し、dataset manifestを生成
