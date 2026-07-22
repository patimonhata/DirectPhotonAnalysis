# 再現メモ

これは厳密な実験ノートではなく、「このROOT fileを何の設定で作ったか」を思い出すための短い索引です。
本番生成や比較条件を変えたときだけ1エントリ追加します。細かすぎて負担になる場合は、不要な欄を削除してください。

## 記録するタイミング

- generator の energy / eta / particle を変えた
- CDB、calibration、cluster threshold、split設定を変えた
- 大量のCondor jobを投入した
- 論文・スライド・比較図に使うROOT/PDFを作った
- 後から同じsampleを作り直す可能性がある

## テンプレート

以下をコピーして、この文書の先頭側へ追加します。

```markdown
## YYYY-MM-DD: 短い名前

- 目的:
- git commit:
- sPHENIX release:
- 入力:
- generator: particle=, pt/E=, eta=, phi=, vertex=
- reconstruction: tower node=, cluster node=, split=, threshold=, vertex mode=
- jobs: processID=, events/job=, number of jobs=, job_offset=
- 実行コマンド:
- 主な出力:
- 確認結果: entries=, matched events=, mass peak=, その他=
- 注意・未解決:
```

すべて埋める必要はありません。最低限 `git commit`、入力、変更した条件、主な出力の4項目があれば、
コードと生成物を結び付けられます。

## 2026-07-19: ドキュメント作成時点の設定スナップショット

これは完了したproductionの記録ではなく、現行コードの既定値を例として書いたものです。

- 目的: 現行の単一pi0生成とsplit / no-split比較設定を記録
- git commit: `d595650` を基準。ドキュメント作成時の未コミット変更あり
- sPHENIX release: wrapperでは `ana`。固定releaseではない
- 入力: generator生成のため外部入力なし
- generator: `particle=pi0`, `pt=5 GeV`, `eta=0.5`, `phi=[-pi,pi]`, `vertex=(0,0,0) cm`
- reconstruction: `TOWERINFO_CALIB_CEMC`, split=`CLUSTERINFO_CEMC`, no-split=`CLUSTERINFO_CEMC_NO_SPLIT`
- no-split cluster builder threshold: `0.070 GeV`
- pi0 reconstruction: vertex mode=`Origin`, min cluster energy=`0.0 GeV`, mass histogram=`100 bins, 0--1 GeV`
- ローカル確認例: `processID=0`, `events=5`
- simulation出力: `SinglePi0GunSimulation/output/DST_pi0_5GeV_eta05_towerinfo/DST_single_pi0_reconstructedInfo_000000.root`
- reconstruction出力: `Pi0Reconstruction/output/root/pi0_reconstruction_{SPLIT,NO_SPLIT}_000000.root`
- 注意: Pi0Reconstructionの既定入力は `eta0` なので、現行simulationの `eta05` を読む場合は入力ディレクトリ引数を明示する
- 注意: EventDisplayの現行既定入力は別プロジェクトのsingle-gamma sampleで、このpi0手順とは独立
