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

wrapperは最後に`check_scored_tree.C`を実行し、cluster、tower、pair、score vectorの長さに加え、`metadata` TTreeが実際に読めること、event数とsource file IDが一致することを検証します。既存outputは上書きしません。Condorで多数jobを同時実行してもACLiCの共有build fileが競合しないよう、score macroは各process内でloadします。

## scored fileをmergeする

merged fileは個別fileと同じdirectoryへ置かず、入力に6桁のprocess IDを明示したpatternを使います。

```bash
mkdir -p PhotonAnalysisTree/output/merged
hadd -f PhotonAnalysisTree/output/merged/all.root \
  PhotonAnalysisTree/output/root/photon_analysis_tree_[0-9][0-9][0-9][0-9][0-9][0-9]_scored.root
```

正常なmergeでは`event_tree`と`metadata`がそれぞれmergeされ、`metadata`のentry数は入力file数になります。

## Model provenanceと制約

同梱`models/best_model.onnx`のSHA-256は`194eb9a0394b5752619d324aa8a17dde547ff67713ea908dcfed869ef0931a6d`です。元は`Pi0DirectGammaSeparation/Pi0GammaClassifier/output/baseline`で、`models/onnx_metadata.json`も同梱しています。

ONNX modelは`n_cluster == 1`のNO_SPLIT eventだけで学習されています。adapterはmulti-cluster eventもclusterごとに評価しますが、そのscoreは学習条件外なので物理解析で使う前に別途validationが必要です。

BDT modelのdocumented performance binはcluster ET 6 GeVからです。それより低いETのscoreはsoftware上は計算されますが、性能保証範囲外です。

## HTCondor production

`condor/`にはtree生成用、score付与用、および両段階を順番に実行するDAGがあります。

全5000 DSTを一括で実行する場合:

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree
./condor/submit_condor.sh pipeline
```

`pipeline.dag`は全5000件のtree jobが成功してからscore batchを開始します。treeが1件でも失敗・holdした場合、score段階は開始しません。

小規模testや部分再実行では段階別にsubmitできます。次はprocess ID 500--599を各DSTの先頭20 eventだけ処理する例です。

```bash
./condor/submit_condor.sh tree 100 500 20
# tree jobsの完了を確認してから
./condor/submit_condor.sh scores 100 500
```

引数は`STAGE N_JOBS JOB_OFFSET N_EVENTS`です。score stageでは`N_EVENTS`を使いません。submit helperは必要なoutput directoryを作り、対象範囲に既存outputがある場合やscore入力treeが欠ける場合は、submit前に停止します。部分的なBDT/scored fileが残ったjobを再実行するときは、内容を確認してからそのfileを明示的に整理してください。

Condor files:

- `condor/run_tree.job`: `run_tree.sh PROCESS_ID N_EVENTS`
- `condor/run_add_scores.job`: process IDからinput/BDT/final ROOT pathを構成
- `condor/pipeline.dag`: `TREE -> SCORES`のbatch依存
- `condor/submit_condor.sh`: directory作成、範囲preflight、submit
