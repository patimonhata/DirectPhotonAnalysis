# Run 2 pp Good Run List for isolated prompt photon A_N

このディレクトリにある isolated prompt photon A_N 用の good run list は、
- sPHENIX Calorimeter Golden Run List v3 に含まれる
- かつ、Spin QA を pass する
という条件で作られました。


## sPHENIX Calorimeter Golden Run List v3

`Full_ppGoldenRunList_Version3.list` は sPHENIX Calorimeter Golden Run List v3 のコピーです。
元コミットは `cd5bf373b6c87b0160496a75adf8f87969e0bd83` で、1,568 run が含まれます。
[上流の golden run list](https://github.com/sPHENIX-Collaboration/RunList/blob/cd5bf373b6c87b0160496a75adf8f87969e0bd83/Physics/Calo/goodRuns/grls/run2pp/Full_ppGoldenRunList_Version3.list)

## Spin QA run list

`PhotonAnalysisTree/workflows/physics_data/select_spin_qa_runs.C` を使って、sPHENIX Calorimeter Golden Run List v3 に含まれる全 run に対して Spin QA のチェックが行われました。

```bash
[]$ root PhotonAnalysisTree/workflows/physics_data/select_spin_qa_runs.C
```

入力: `PhotonAnalysisTree/input/physics_data/Full_ppGoldenRunList_Version3.list`  
出力: `PhotonAnalysisTree/input/physics_data/run2pp_spin_qa.list`

マクロは golden list の各 run について Spin DB の default 行（`is_default = TRUE`）を参照し、`GetBadRunFlag() == 0` の run を `run2pp_spin_qa.list` に昇順で出力します。
実行結果に表示される採用・非ゼロフラグ・default 行欠損の件数と、Spin DB の照会日をここに記録します。
再実行が成功すると既存の出力リストを置き換えます。

## DST リストの作成(`Jet` と `JETCALO`)

上で作成された run number のリストをもとに、`DST_Jet` と `DST_JETCALO` のファイルリストを作成しました。
Production tag は `ana521_2025p007_v001` を指定しました。

```bash
mkdir -p PhotonAnalysisTree/input/physics_data/Jet_and_JETCALO_ana521_2025p007_v001
cd PhotonAnalysisTree/input/physics_data/Jet_and_JETCALO_ana521_2025p007_v001
CreateDstList.pl --tag ana521_2025p007_v001 --list ../run2pp_spin_qa.list DST_Jet DST_JETCALO
```

出力は `dst_jet-<8桁のrun番号>.list` と `dst_jetcalo-<8桁のrun番号>.list` です。
両方を同じジョブで使う場合は、二つのリストで共通する `(run, segment)` を使います。
