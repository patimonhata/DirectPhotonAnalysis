# ビルドと production の対応

通常の開発用ビルドと、正式 production が使う固定ソフトウェアを分けます。既存の `install/` 直下と旧 build は保管し、新しいビルドの既定先として使いません。

## 開発用ビルド

repository root から実行します。

```bash
bash PhotonAnalysisTree/src/build.sh
```

既定 release は `ana.565` です。rolling alias の `ana` は受け付けません。Pi0Reconstruction と PhotonAnalysisTree の両方を同じ release でビルドし、次に配置します。

```text
PhotonAnalysisTree/
  build/dev/ana.565/
    pi0reconstruction/
    photon/
  install/dev/ana.565/
    build-manifest.json
    include/
    lib64/
    share/PhotonAnalysisTree/models/
    dependencies/Pi0Reconstruction/
      include/pi0reconstruction/
      lib/
```

Pi0Reconstruction の新しい CMake 手順は、この管理ビルドで使用します。既存の autotools build/install は変更しません。
`build-manifest.json` は Git commit・dirty 状態・source hash・実際の OFFLINE_MAIN・ROOT/compiler version・ビルド設定・配置後の artifact SHA-256・解決したライブラリ依存を記録します。
開発中の dirty source はビルドできます。正式 software の発行には使えません。

固定 release と出力先は明示的に変更できます。

```bash
bash PhotonAnalysisTree/src/build.sh --release ana.565 --jobs 4
# 別の検証先（共通 install/ や software/<id>/ は指定しない）
bash PhotonAnalysisTree/src/build.sh --build-dir /tmp/photon-build --install-prefix /tmp/photon-install
```

環境・ソース・install 先が異なる CMake cache は再利用しません。release を変える場合は別の build ディレクトリを使います。
旧 `build/current`、`candidate`、`eta_topology`、`topology_v11` は当時の作業領域です。名前だけでソース版が固定されているとは扱いません。
`build/layout_migration` と `recovery_v11_schema5` は過去の検証・復元記録なので、通常ビルドの整理対象と区別します。今回これらは移動・削除していません。

## 正式 software の作成

実行コードをコミットして作業ツリーを clean にした後、ビルドし直します。clean な現在のソースとビルド記録が一致した場合にだけ freeze できます。

```bash
bash PhotonAnalysisTree/src/build.sh
python3 PhotonAnalysisTree/tools/software.py freeze candidate_v1_ana565
python3 PhotonAnalysisTree/tools/software.py verify PhotonAnalysisTree/software/candidate_v1_ana565
```

```text
software/<software-id>/
  manifest.json
  include/・lib64/・share/・dependencies/    # 配置後のライブラリ・ヘッダ・モデル
  source/
    PhotonAnalysisTree/                   # src・analysis・workflows・validation・tests・tools 等
    Pi0Reconstruction/src/
    TopologyEventDisplay/
    ...
```

input list はデータとして扱い、ソフトウェア snapshot には複製しません。map 実行時に入力 manifest を明示します。
同じ software ID は上書きしません。実行前に保存ファイルを hash 検証し、改変や欠落・未記録ファイルがあれば停止します。
異なる解析条件の production が同じ software ID を共有することはできます。

ROOT のヘッダと library は選択した prefix の `ROOT_INCLUDE_PATH` / `LD_LIBRARY_PATH` から読みます。
`libPhotonAnalysisTree.so` の install RPATH は同梱 Pi0Reconstruction への相対パスです。実行時の sPHENIX release は software manifest から選びます。
reduce / merge も保存した source のコードを使うため、作業ツリーを変更しても正式 production のコードは切り替わりません。

## 新しい production の作成

保存済み結果と異なる、新しい configuration ディレクトリを指定します。

```bash
python3 PhotonAnalysisTree/tools/software.py production \
  PhotonAnalysisTree/software/candidate_v1_ana565 \
  PhotonAnalysisTree/output/photon_candidate_selection/productions/new_production/new_configuration/manifest.json \
  --production new_production --configuration new_configuration
```

manifest は software の絶対パス・ID・manifest SHA-256・Git commit・release を結び付けます。
`--settings settings.json` で解析条件の説明を保存できます。実行で実際に使った引数は出力ごとの execution record と ROOT metadata に記録します。
正式 production の結び付けと stage 履歴は、まず主解析 `photon_candidate_selection` を対象とします。その他の workflow は固定環境の開発用 runtime を共用します。
既存結果を持つディレクトリや既存 manifest は上書き・再ラベルしません。過去結果のビルド版は確認できる範囲だけを履歴として扱い、現在の install から推測して埋めません。

主解析の Condor submit は `production_directory` を必須の運用設定とします。既定の `UNCONFIGURED` では有効な production を実行できません。
以下は投入例です。スクリプトは自動投入しません。

```bash
condor_submit -append "production_directory = /absolute/path/to/new_configuration" \
  PhotonAnalysisTree/workflows/photon_candidate_selection/submit_jet5.job
```

同じ設定を reduce の submit にも渡します。各 submit の物理 cut は引き続き確認してください。Jet と PhotonJet は別 configuration に分けます。

ローカルでは、既存の引数をそのまま渡せます。

```bash
python3 PhotonAnalysisTree/tools/software.py run \
  --production-manifest /absolute/path/to/new_configuration/manifest.json \
  merge jet /absolute/path/to/new_configuration/partials/jet /absolute/path/to/new_configuration/results/jet
```

`run --software SOFTWARE_DIRECTORY map|reduce|merge ...` は software だけを選ぶ QA 用です。
production を選ばない通常の `run_*.sh` は、既定の開発用 install を使います。
別の開発 prefix を使う場合は `PHOTON_INSTALL_PREFIX`、release を変更する場合は `PHOTON_RELEASE` を明示します。

## 出力の履歴と入力の混在防止

主解析の shell wrapper は、出力 ROOT ごとに `<filename>.root.execution.json` を作ります。
記録内容は stage・引数・出力 SHA-256・software manifest SHA-256・source 版・実際の runtime release・production manifest SHA-256 です。
map の一時名は正式出力名と分けて扱い、出力名の公開時に execution record も移動します。

正式 production の reduce / merge は、入力の execution record と ROOT hash を確認し、異なる software / production の入力、改変された入力、履歴がない入力を拒否します。
旧 map を使う過去結果の再現は `run_reference.sh` 等の専用入口を使います。旧結果を新しい production として自動登録しません。

## 検証

```bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
bash PhotonAnalysisTree/tests/run_candidate_selection_tests.sh
python3 PhotonAnalysisTree/tests/test_software.py
```

後者は rolling release の拒否、artifact 改変・欠落、clean-source gate、既存結果の保護、production のソフトウェア切替拒否、stage 間の入力混在を検証します。

## 今回の移行で確認した範囲

- ana.565 で両ローカル library を新規ビルドし、既存 candidate / tagging mass テストと provenance の回帰テストを確認した。
- 実 DST 1 segment・2イベントの map → validator → reduce が成功し、従来の install による結果と histogram 2,202個と tree 内容が一致した。出力先とランダムな一時ファイル名だけを正規化した。
- 一時 clone に clean な版を作って freeze と production 初期化を確認した。元の開発 install を隠し、作業ツリーの reducer / merger を読み込めなくした状態でも保存したコードで実行できた。本作業の未コミット状態を正式 software として発行してはいない。
- 保存版で4 sample を reduce → merge し、全26 ROOT に execution record が揃った。従来結果と histogram 12,654個の全 bin 内容・誤差・軸、44 tree の内容、CSV 15ファイルが一致した。PDF 246ファイルは配置とページ数を確認した。
- 全28 Condor submit の dry-run と、過去 study の12再現設定の入力確認が通った。ジョブは投入していない。
- 既存 install 直下の両 library の SHA-256 は変化していない。過去 build、保存成果物、外部 map は保持した。
