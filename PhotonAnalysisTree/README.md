# PhotonAnalysisTree

Pythia DSTには別frontendの`PythiaPhotonAnalysisTree`を使います。cluster/tower/shower-shape
処理はsingle-particle版と共有しますが、event truthは複数primaryを前提としたcluster contributor
schemaです。詳細は[docs/pythia_tree_schema.md](docs/pythia_tree_schema.md)を参照してください。

Single-particle gun DSTを1回だけ読み、truth、SPLIT cluster、NO_SPLIT cluster、再構成gamma pairを同じ`event_tree`へ保存するFun4All subsystemです。`event_tree`は1 entry = 1 input eventです。

## なぜ独立ディレクトリか

`Pi0Reconstruction`はcluster collectionからpairを再構成する部品です。このpackageはtruth、2種類のcluster collection、ML用tower point set、2つの推論器を束ねる上位のI/O pipelineなので、`DirectPhotonAnalysis/PhotonAnalysisTree`として分離しています。

重複しそうなcluster branchは削除せず、物理的に異なるcollectionとして`split_`と`nosplit_` prefixで分けます。event ID、vertex、truthは共通なので1組だけです。両collectionは必ずcluster energy降順、同値ならcluster ID昇順です。全cluster vector、pair index、score vectorはこの順序を共有します。

## 内容

- truth: primary、direct daughter、CEMC半径へのgamma投影、acceptance、truth mass/asymmetry
- `split_*`: `CLUSTERINFO_CEMC`、shower shape、全cluster pair、cluster constituent tower
- `nosplit_*`: `CLUSTERINFO_CEMC_NO_SPLIT`、shower shape、全cluster pair、cluster constituent tower
- `split_cluster_bdt_base_v3E_*`: SPLIT clusterで学習した`base_v3E` BDTをSPLIT clusterへ適用
- `split_cluster_bdt_ppg15v1_*`: low-pT sampleを加えて学習したPPG15v1 BDTをSPLIT clusterへ適用
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
basenameで登録し、schema version 4の専用treeを作ります。

- `DST_CALO_CLUSTER`
- `DST_MBD_EPD`
- `DST_TRUTH_JET`
- `G4Hits`

```bash
root -l -b -q \
  'PhotonAnalysisTree/macro/Fun4All_PhotonAnalysisTreePythia.C("pythia8_Jet5-0000000028-000000.root",0)'
```

第1引数はstream prefixを除いた共通suffix、第2引数はevent数（`0`は全event）です。
macroは4つのstream prefixを付けたbasenameを各input managerへ登録し、Fun4Allの
file lookupへ渡します。出力名は
`pythia_photon_analysis_tree_<input suffix without .root>.root`です。
G4 truth、`PHHepMCGenEventMap`、SPLIT clusterは必須、NO_SPLIT clusterはoptionalです。
wrapperはsingle-particle版と分離しており、生成後に`check_pythia_tree.C`で
schema version、branch alignment、contributor offsetとfractionを検証します。

```bash
./PhotonAnalysisTree/run_tree_pythia.sh \
  pythia8_Jet5-0000000028-000000.root 0
```

`input/jet5`の4つのstream別listからCondor用の共通suffix manifestを生成するには
次を実行します。4リストの行数、順序、suffixが一致しない場合は生成に失敗します。

```bash
./PhotonAnalysisTree/make_pythia_input_manifest.sh
```

## scoreを追加する

single-particle用`run_add_scores.sh`は入力ROOT fileの`event_tree`へBDT/gamma score branchを直接追加し、全段と最終検証が成功した後に`_scored.root`へrenameします。途中で失敗した場合は入力ファイルを元DSTから作り直してください。

```bash
./PhotonAnalysisTree/run_add_scores.sh \
  PhotonAnalysisTree/output/root/photon_analysis_tree_000000.root \
  PhotonAnalysisTree/output/root/photon_analysis_tree_000000_scored.root
```

wrapperは最後に`check_scored_tree.C`を実行し、cluster、tower、pair、score vectorの長さに加え、`metadata` TTreeが実際に読めること、event数とsource file IDが一致することを検証します。scored ROOT fileのtop-level objectは`event_tree`と`metadata`だけです。adapterのcluster数、valid数、malformed数などの集計は標準出力（Condor log）だけに記録します。既存の最終outputは上書きせず、score/valid branchが既にある入力もエラーにします。Condorで多数jobを同時実行してもACLiCの共有build fileが競合しないよう、score macroは各process内でloadします。

Pythia schema version 4にはsplit専用wrapperを使います。第3引数には
split clusterで学習したONNX modelを必ず明示します。

```bash
./PhotonAnalysisTree/run_add_scores_pythia.sh \
  PhotonAnalysisTree/output/root/pythia_photon_analysis_tree_pythia8_Jet5-0000000028-000000.root \
  PhotonAnalysisTree/output/root/pythia_photon_analysis_tree_pythia8_Jet5-0000000028-000000_scored.root \
  /path/to/split_model.onnx
```

このwrapperは入力を変更せず、最終outputと同じdirectoryの一時コピーへ次のbranchだけを
追加します。

- `split_cluster_bdt_base_v3E_{score,valid}`
- `split_cluster_bdt_ppg15v1_{score,valid}`
- `split_cluster_p_gamma{,_valid}`

変更前後に`check_pythia_tree.C`、最後に`check_pythia_scored_tree.C`を実行し、
全検証成功後だけ一時ファイルを最終outputへrenameします。no-split score branchは
作りません。専用Condor jobでは、split modelがまだない間の暫定設定として、
NO_SPLIT・1 cluster条件で学習された同梱`models/best_model.onnx`を明示的に指定します。
このONNX scoreは学習domain外であり、物理解析へ使用する前に別途validationが必要です。
wrapperを直接呼ぶ場合は引き続き第3引数でmodel pathを明示します。

`add_split_bdt_ppg15v1.C`は既存SPLIT BDTと同じ11 featureを使い、low-pT sampleを追加したモデルの出力を`split_cluster_bdt_ppg15v1_score`、入力・shower shape・推論が有効かを`split_cluster_bdt_ppg15v1_valid`へ保存します。wrapperの第7引数でdefault model pathを上書きできます。

NO_SPLIT clusterで学習したBDTも`run_add_scores.sh`内でNO_SPLIT clusterへ適用します。default modelは`model_base_v3E_nosplit_single_tmva.root`です。必要に応じて独立macroとしても実行できます。

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

既存`base_v3E` BDT modelのdocumented performance binはcluster ET 6 GeVからです。それより低いETのscoreはsoftware上は計算されますが、性能保証範囲外です。PPG15v1 modelは同じfeatureを使い、よりlow-pTのsampleを含めて学習しています。

## HTCondor production

実行scriptとjob fileは`PhotonAnalysisTree`直下にあります。

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree
condor_submit run_tree.job
condor_submit run_tree_pythia.job
# tree jobsの完了後
condor_submit run_add_scores.job
# manifestの先頭から1 jobだけsubmit
condor_submit -append "manifest_slice = [0:1]" -maxjobs 1 run_add_scores_pythia.job
```

`run_tree.job`のjob数とoffsetはsubmit時に上書きできます。`run_add_scores.job`はdefaultで5000 jobをqueueするため、job数を変える場合はjob file末尾の`Queue 5000`を編集し、offsetだけをsubmit時に上書きします。

```bash
condor_submit -append "n_jobs = 100" -append "job_offset = 500" run_tree.job
condor_submit -append "job_offset = 500" run_add_scores.job
```

Pythia scoringを段階的に増やす場合はmanifestとsliceをsubmit時に上書きできます。

```bash
condor_submit -append "input_manifest = input/jet5/segments.list" -append "manifest_slice = [0:10]" -maxjobs 10 run_add_scores_pythia.job
```

- `run_tree.sh`, `run_tree.job`: single-particle DSTからbase TTreeを生成
- `run_tree_pythia.sh`, `run_tree_pythia.job`: Pythia 4-stream DSTからschema 4 treeを生成
- `make_pythia_input_manifest.sh`: 4つのstream別listを検証して共通suffix manifestを生成
- `run_add_scores.sh`, `run_add_scores.job`: BDT/gamma score branchを追加して検証
- `run_add_scores_pythia.sh`, `run_add_scores_pythia.job`: Pythia schema 4 treeへsplit scoreだけをtransactionalに追加して検証
- `make_dataset_manifest.sh`: `hadd`後のROOT fileとmodelを検証し、dataset manifestを生成

## Minimum-bias truth pT spectra

cluster reconstructionを必要としないtruth粒子数の測定には、G4Hits DSTだけを読む軽量な
`PythiaTruthSpectrumTree`を使います。prompt photon候補とpi0はsignal HepMC eventから保存します。
pi0 decay photonのpT spectrumは、HepMC内で崩壊したpi0由来のfinal photonと、
Geant4に委譲されたsignal-primary pi0由来のphotonを合算します。HepMC/G4別の
diagnostic histogramも保存します。詳細な粒子定義とbranchは
`docs/pythia_truth_spectrum_schema.md`を参照してください。

4 stream listのsuffix同期を検証してfull manifestを作ります。

```bash
cd /sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree
./make_pythia_input_manifest.sh input/minimum_bias input/minimum_bias/segments.list
```

1 fileの先頭eventだけを確認する場合:

```bash
./run_truth_spectrum_pythia.sh \
  pythia8_Detroit-0000000028-000000.root 1 output/truth_root
```

defaultのCondor jobはmanifestのhalf-open slice `[0:1000]`、つまり先頭1000 filesだけを
queueします。`-maxjobs`は切り詰めではなく誤投入防止です。

```bash
condor_submit -maxjobs 1000 run_truth_spectrum_pythia.job
```

少数の生成ROOT filesなら直接TChainへ読み、count densityを描けます。photonは
classifier category 1 (direct)または2 (fragmentation)のprompt photonです。defaultは
eta cutなしで、最後から2つ目の引数を`0.7`にすると`|eta_truth| < 0.7`を3粒子種へ適用します。
pi0 decay photonの合計とHepMC/G4内訳は同じROOT fileに保存されます。

```bash
root -l -b -q 'macro/PlotPythiaTruthPtSpectra.C("output/truth_root/pythia_truth_spectrum_tree_*.root","output/plots/minbias_truth_pt",100,20.0,-1.0,false)'
```

数万files以上では、1 processで全fileを開かずmap-reduce型で処理します。defaultは
先頭40,000 truth treesを500 filesずつ、80 partial jobsへ分割します。

```bash
condor_submit -maxjobs 80 run_truth_pt_partial_pythia.job
```

全200,000 filesを最初から処理する場合:

```bash
condor_submit \
  -append "total_files = 200000" \
  -append "n_chunks = 400" \
  -maxjobs 400 \
  run_truth_pt_partial_pythia.job
```

先頭40,000 filesのpartialが完成済みなら、chunk 80から残りだけを追加できます。

```bash
condor_submit \
  -append "total_files = 200000" \
  -append "chunk_offset = 80" \
  -append "n_chunks = 320" \
  -maxjobs 320 \
  run_truth_pt_partial_pythia.job
```

partialの範囲・解析条件を検証して統合し、最後にbin幅で規格化してplotします。
`expected_manifest_end`を指定するため、partialの欠落も検出します。

```bash
root -l -b -q 'macro/FinalizePythiaTruthPtSpectra.C("output/truth_pt_partial/prompt_eta07_unweighted_inclusive_pi0_decay/partial_*.root","output/plots/minbias_truth_pt_prompt_eta07_inclusive_pi0_decay",0,40000)'
```

詳細は`docs/pythia_truth_pt_map_reduce.md`を参照してください。

## Minimum-bias cluster E_T spectra

再構成clusterのprompt/pi0起源とseparated/merged/missing topologyは、4 stream DSTを直接読む`PythiaClusterEtSpectrum`でmap-reduceします。各G4 photonを実際の生成vertexからCEMCへ投影するため、displaced decayもcollision vertex近似には戻しません。default cut、pathway定義、実行方法は[docs/pythia_cluster_et_map_reduce.md](docs/pythia_cluster_et_map_reduce.md)を参照してください。最終density histogramの単位は`clusters/GeV`です。

`run_tree_pythia.job`は`input/jet5/segments.list`を`queue ... from`で読み、
manifestの1行につき1つの同期した4-stream DST jobを生成します。
