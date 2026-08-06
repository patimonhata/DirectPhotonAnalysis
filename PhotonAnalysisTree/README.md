# PhotonAnalysisTree

Single-particle gun DSTを1回だけ読み、truth、SPLIT cluster、NO_SPLIT cluster、再構成gamma pairを同じ`event_tree`へ保存するFun4All subsystemです。`event_tree`は1 entry = 1 input eventです。

## なぜ独立ディレクトリか

`Pi0Reconstruction`はcluster collectionからpairを再構成する部品です。このpackageはtruth、2種類のcluster collection、ML用tower point set、2つの推論器を束ねる上位のI/O pipelineなので、`DirectPhotonAnalysis/PhotonAnalysisTree`として分離しています。

重複しそうなcluster branchは削除せず、物理的に異なるcollectionとして`split_`と`nosplit_` prefixで分けます。event ID、vertex、truthは共通なので1組だけです。両collectionは必ずcluster energy降順、同値ならcluster ID昇順です。全cluster vector、pair index、score vectorはこの順序を共有します。

## 内容

- truth: primary、direct daughter、CEMC半径へのgamma投影、acceptance、truth mass/asymmetry
- `split_*`: `CLUSTERINFO_CEMC`、shower shape、全cluster pair、cluster constituent tower
- `nosplit_*`: `CLUSTERINFO_CEMC_NO_SPLIT`、shower shape、全cluster pair、cluster constituent tower
- `split_cluster_bdt_base_v3E_*`: SPLIT clusterで学習した`base_v3E` BDTをSPLIT clusterへ適用
- `nosplit_cluster_bdt_base_v3E_*`: NO_SPLIT clusterで学習した`base_v3E` BDTをNO_SPLIT clusterへ適用
- `nosplit_cluster_p_gamma*`: 同梱ONNX modelが後段で追加
- `split_cluster_p_gamma*`: 将来のSPLIT学習ONNX modelを`add_split_gamma_onnx.C`で追加

完全なbranch契約は[docs/tree_schema.md](docs/tree_schema.md)を参照してください。

TruthまたはNO_SPLIT clusterを持たないDSTでは、それぞれ`set_require_truth_node(false)`、
`set_require_nosplit_cluster_node(false)`を設定できます。省略collectionのbranchは削除せず、
invalid値または空vectorとして保存します。

## Build

先に既存`Pi0Reconstruction`をbuild/installしてから実行します。

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis
./PhotonAnalysisTree/src/build.sh
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

Pythia Jet5では同じsegment IDを持つ次の4 streamを、それぞれ別のinput managerへ
basenameで登録します。

- `DST_CALO_CLUSTER`
- `DST_MBD_EPD`
- `DST_TRUTH_JET`
- `G4Hits`

```bash
root -l -b -q \
  'PhotonAnalysisTree/macro/Fun4All_PhotonAnalysisTreePythia.C(0,0)'
```

第1引数はsegment ID、第2引数はevent数（`0`は全event）です。出力名は
`photon_analysis_tree_<segment>.root`です。Truth nodeは必須、NO_SPLIT clusterはoptionalです。

## scoreを追加する

入力ROOT fileの`event_tree`へBDT/gamma score branchを直接追加し、全段と最終検証が成功した後に`_scored.root`へrenameします。途中で失敗した場合は入力ファイルを元DSTから作り直してください。

```bash
./PhotonAnalysisTree/run_add_scores.sh \
  PhotonAnalysisTree/output/root/photon_analysis_tree_000000.root \
  PhotonAnalysisTree/output/root/photon_analysis_tree_000000_scored.root
```

wrapperは最後に`check_scored_tree.C`を実行し、cluster、tower、pair、score vectorの長さに加え、`metadata` TTreeが実際に読めること、event数とsource file IDが一致することを検証します。scored ROOT fileのtop-level objectは`event_tree`と`metadata`だけです。adapterのcluster数、valid数、malformed数などの集計は標準出力（Condor log）だけに記録します。既存の最終outputは上書きせず、score/valid branchが既にある入力もエラーにします。Condorで多数jobを同時実行してもACLiCの共有build fileが競合しないよう、score macroは各process内でloadします。

NO_SPLIT clusterで学習したBDTは`run_add_scores.sh`の第2段でNO_SPLIT clusterへ適用します。default modelは`model_base_v3E_nosplit_single_tmva.root`です。必要に応じて独立macroとしても実行できます。

```bash
root -l -b \
  -e '.L PhotonAnalysisTree/macro/add_nosplit_bdt.C' \
  -e 'gSystem->Exit(add_nosplit_bdt("input.root"));'
```

このmacroが追加するのは`nosplit_cluster_bdt_base_v3E_score`と`nosplit_cluster_bdt_base_v3E_valid`です。入力featureの定義と順序は既存のSPLIT学習版macroと同一ですが、値はすべて`nosplit_cluster_*` branchから取得します。

SPLIT clusterで学習したONNX modelが用意できた後は、model pathを明示して次のmacroを使います。このmacroには意図的にdefault modelを設定していません。

```bash
root -l -b \
  -e '.L PhotonAnalysisTree/macro/add_split_gamma_onnx.C' \
  -e 'gSystem->Exit(add_split_gamma_onnx("input.root", "split_model.onnx"));'
```

入力にはschema version 2以降の`split_tower_*` branchが必要です。SPLITでは同じtowerが複数clusterへ分配され得るため、point featureのenergyには未分配の`split_tower_energy`ではなく、RawClusterに保存された割当energy `split_tower_cluster_value`を使います。

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

`run_tree.job`のjob数とoffsetはsubmit時に上書きできます。`run_add_scores.job`はdefaultで5000 jobをqueueするため、job数を変える場合はjob file末尾の`Queue 5000`を編集し、offsetだけをsubmit時に上書きします。

```bash
condor_submit -append "n_jobs = 100" -append "job_offset = 500" run_tree.job
condor_submit -append "job_offset = 500" run_add_scores.job
```

- `run_tree.sh`, `run_tree.job`: DSTからbase TTreeを生成
- `run_add_scores.sh`, `run_add_scores.job`: BDT/gamma score branchを追加して検証
- `make_dataset_manifest.sh`: `hadd`後のROOT fileとmodelを検証し、dataset manifestを生成

`run_tree.job`はPythia Jet5のsegment 0--999をdefaultでqueueします。
