# Single-particle gun workflow

Single-particle gun DSTから解析用TTreeを作り、photon-ID scoreを追加してdatasetへまとめるworkflowです。branch契約は[tree schema](../../docs/tree_schema.md)を参照してください。

## Pipeline

1. run_tree.sh / submit_tree.job: DSTからscoreなしのbase TTreeを生成
2. check_tree.C: metadata、event UID、truth、cluster、tower、pair、truth-match vectorを検証
3. run_add_scores.sh / submit_scores.job: 5系統のscore branchをtransactionalに追加
4. check_scored_tree.C: 全score vectorとbase treeの主要な整合性を検証
5. hadd: scored filesのevent_treeとmetadataをmerge
6. make_dataset_manifest.sh: merge結果、event数、model hash、Git revisionをmanifestへ記録

低水準のscore macro (add_*.C) はPythia workflowとも共有するため、PhotonAnalysisTree/macroに残しています。

## Base TTree

repository rootから1 segmentだけ実行する例です。第1引数はsource/process ID、第2引数はevent数で、0は入力をEOFまで処理します。

~~~bash
PhotonAnalysisTree/workflows/single_particle/run_tree.sh 0 0
~~~

defaultの入力はSinglePi0GunSimulation/output/DST_pi0_3to15GeV_etapm1_vertexpm60、出力はPhotonAnalysisTree/output/root/photon_analysis_tree_<ID>.rootです。入力directory、出力directory、primary PDGを変える場合はFun4All_PhotonAnalysisTree.Cを直接呼びます。

Condor production:

~~~bash
cd PhotonAnalysisTree/workflows/single_particle
condor_submit submit_tree.job
~~~

job数とoffsetはsubmit時に上書きできます。

~~~bash
condor_submit -append 'n_jobs = 100' -append 'job_offset = 500' submit_tree.job
~~~

生成後にcheck_tree.Cが自動実行され、失敗したCondor jobはholdされます。

## Score追加

~~~bash
PhotonAnalysisTree/workflows/single_particle/run_add_scores.sh \
  PhotonAnalysisTree/output/root/photon_analysis_tree_000000.root \
  PhotonAnalysisTree/output/root/photon_analysis_tree_000000_scored.root
~~~

wrapperはbase fileを事前検査し、final outputと同じdirectoryに一時コピーを作ります。そのコピーへsplit/no-split BDT、PPG15v1 split BDT、split/no-split ONNX scoreを追加し、base/scored checkerがすべて成功した場合だけfinal output名へrenameします。入力fileは変更しません。失敗時は一時fileを削除し、既存のfinal outputは上書きしません。

Condor productionはtree jobs完了後に実行します。

~~~bash
cd PhotonAnalysisTree/workflows/single_particle
condor_submit submit_scores.job
~~~

score wrapperの任意引数でmodel pathを上書きできます。default modelと適用範囲の注意事項は[PhotonAnalysisTree README](../../README.md#model-provenanceと制約)に記載しています。

## Mergeとmanifest

scored fileだけを明示的な6桁ID patternでmergeします。

~~~bash
mkdir -p PhotonAnalysisTree/output/merged
hadd -f PhotonAnalysisTree/output/merged/all.root \
  PhotonAnalysisTree/output/root/photon_analysis_tree_[0-9][0-9][0-9][0-9][0-9][0-9]_scored.root
PhotonAnalysisTree/workflows/single_particle/make_dataset_manifest.sh
~~~

default manifestはPhotonAnalysisTree/output/merged/manifest.jsonです。引数はSCORED_DIR MERGED_ROOT MANIFEST_JSON BDT_MODEL ONNX_MODELの順です。既存manifestを置き換える場合だけFORCE=1を指定します。

copy_analysis_trees.Cは、古いROOT fileに余分なtop-level objectやkey cycleがある場合に、event_treeとmetadataだけをclean copyする補助macroです。通常の新規productionでは必須ではありません。

## Files

- Fun4All_PhotonAnalysisTree.C: single-particle DST frontend
- run_tree.sh, submit_tree.job: base tree production
- check_tree.C: base schema checker
- run_add_scores.sh, submit_scores.job: transactional scoring
- check_scored_tree.C: scored schema checker
- copy_analysis_trees.C: optional clean-copy helper
- make_dataset_manifest.sh: merged dataset manifest
