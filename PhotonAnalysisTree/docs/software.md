# ビルドと production の対応

通常の開発用ビルドは `build/dev/` と `install/dev/`、正式 production 用の固定ソフトウェアは `software/` に置きます。

## 開発用ビルド

repository root から実行します。引数なしで `bash PhotonAnalysisTree/src/build.sh` を実行したときの既定 release は `ana.565` です。rolling alias の `ana` は受け付けません。Pi0Reconstruction と PhotonAnalysisTree の両方を同じ release でビルドし、既定では次に配置します。

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

ローカルで使う版のビルドでは、固定 release と出力先の組を明示できます。

```bash
version=ana.565-ppg15v1
bash PhotonAnalysisTree/src/build.sh --release ana.565 \
  --build-dir "$PWD/PhotonAnalysisTree/build/dev/$version" \
  --install-prefix "$PWD/PhotonAnalysisTree/install/dev/$version"
```

環境・ソース・install 先が異なる CMake cache は再利用しません。release を変える場合は別の build ディレクトリを使います。
`Build cache belongs to another release` は、キャッシュの `PHOTON_CONFIGURED_PREFIX` と setup が返す `OFFLINE_MAIN` の文字列が異なると、コンパイル前に発生します。
同じ `ana.565` でも CVMFS のリンク先変更により表記が変わると、実体が同じパスでも停止します。Git の clean 状態とは別の判定です。
既存キャッシュを保持する場合は、新しい build ディレクトリを指定します。`--build-dir` だけを指定した場合、install 先は既定の `install/dev/ana.565/` です。既存の開発 install も保持する場合は、別の `--install-prefix` を指定します。freeze 時にも同じ install prefix を渡します。
ビルドはキャッシュ検査より前に install 先の `build-manifest.json` を無効化します。失敗後にファイルが残っていても、完成済みのビルドとして freeze することはできません。

`build/topology_v11` は当時の build、`build/layout_migration` と `build/recovery_v11_schema5` は検証・復元の記録です。これらは過去参照用に保持します。

### 開発 install の選択

日常の ROOT 実行では、検証済み install を `install/current` で手動選択できます。上の例を選ぶリンクを初めて作る場合は、次のようにします。

```bash
ln -sT "dev/$version" PhotonAnalysisTree/install/current
```

このリンクを workflow の管理ラッパーでも使う場合は、`PHOTON_INSTALL_PREFIX` に明示します。未指定時の既定は `install/dev/<release>/` です。

```bash
export PHOTON_INSTALL_PREFIX="$PWD/PhotonAnalysisTree/install/current"
python3 PhotonAnalysisTree/tools/software.py resolve
```

## 正式 software の作成

実行コードをコミットして作業ツリーを clean にした後、ビルドし直します。clean な現在のソースとビルド記録が一致した場合にだけ freeze できます。

```bash
version=ana.565-candidate-v1
prefix="$PWD/PhotonAnalysisTree/install/dev/$version"
bash PhotonAnalysisTree/src/build.sh --release ana.565 \
  --build-dir "$PWD/PhotonAnalysisTree/build/dev/$version" --install-prefix "$prefix"
python3 PhotonAnalysisTree/tools/software.py freeze candidate_v1_ana565 --install-prefix "$prefix"
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
ソース snapshot は `git archive HEAD` から作ります。map が読むモデルは `source/PhotonAnalysisTree/models/` に必要なので、モデル本体も Git に収録してください。
`.gitignore` の `*.root` に対して、production 用モデルには個別の例外を追加します。CMake の `install(FILES ...)` に追加するだけではソース snapshot に含まれません。
現在の map モデル・SHA-256・11入力の順序は [candidate definitions](photon_candidate_definitions.md#bdt-model-and-input-order) を参照してください。
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
各 submit は `production_manifest = $(production_directory)/manifest.json` を `PHOTON_PRODUCTION_MANIFEST` に渡します。
`software.py production` は manifest を作成するだけで、呼び出し元のシェル環境変数を設定しません。
実行時は `runtime.sh` の `photon_dispatch` が manifest に結び付いた固定ソースへ切り替え、`photon_setup_runtime` が `PHOTON_SOURCE_ROOT` などを設定します。
通常の production submit に `PHOTON_SOFTWARE_ROOT` の別指定は不要です。この変数は固定 software を直接選択する場合に使い、manifest の結び付けと異なる値は拒否されます。

ローカルでは、既存の引数をそのまま渡せます。

```bash
python3 PhotonAnalysisTree/tools/software.py run \
  --production-manifest /absolute/path/to/new_configuration/manifest.json \
  merge jet /absolute/path/to/new_configuration/partials/jet /absolute/path/to/new_configuration/results/jet
```

`run --software SOFTWARE_DIRECTORY map|reduce|merge ...` は software だけを選ぶ QA 用です。
production を選ばない通常の `run_*.sh` は、既定の開発用 install を使います。
別の開発 prefix を使う場合は `PHOTON_INSTALL_PREFIX`、release を変更する場合は `PHOTON_RELEASE` を明示します。

## map だけを tg に保存する

新しい production の作成時に `--map-root` を指定します。指定先は新規または空のディレクトリに限ります。
production manifest は `/sphenix/user` に置き、外部 map root の実パスを manifest の `output_roots.map` に保存します。
`production/maps` のシンボリックリンクは自動作成するため、手動の `ln -s` は不要です。
map ROOT とその `.execution.json` は tg、partials・results・標準出力/標準エラーのログは `/sphenix/user` に置かれます。
HTCondor の event log（`Log`）はサイト規則により `/tmp/` が必須です。submit ファイルの既定値を使います。

この機能を含むコードをコミットし、ビルド・新規 freeze を済ませてから実行します。
旧 frozen software は新しい設定を解釈できないため、`--map-root` での production 作成を拒否します。
既存 manifest や frozen software の直接編集、既存 map directory の再登録は行いません。
以下は新 software ID を `common_composition_anchor_v14_maproot1_ana565` とした例です。

```bash
module=/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree
production_name=common_composition_anchor_v14_maproot1
configuration=jet_cluster0p12_partner0p5
software="$module/software/${production_name}_ana565"
prod="$module/output/photon_candidate_selection/productions/$production_name/$configuration"
map_root="/sphenix/tg/tg01/coldqcd/ryotaro/DirectPhotonAnalysis/photon_candidate_selection/$production_name/$configuration/maps"

python3 "$module/tools/software.py" production "$software" "$prod/manifest.json" \
  --production "$production_name" --configuration "$configuration" --map-root "$map_root"

mkdir -p "$prod/logs"
for sample in jet3 jet5 jet8 jet12 jet20 jet30 jet40; do
  condor_submit \
    -append "production_directory = $prod" \
    -append "log_directory = $prod/logs" \
    "$module/workflows/photon_candidate_selection/submit_${sample}.job" || break
done
```

reduce は従来どおり `$prod/maps` を読み、`$prod/partials/jet` に保存します。
map が全件完了したら、reduce の submit にも同じ `production_directory` と `log_directory` を指定してください。
reduce でも `Log` は上書きせず、既定の `/tmp/` を使います。`log_directory` は標準出力/標準エラーの保存先だけを変更します。
merge は `run --production-manifest "$prod/manifest.json" merge jet "$prod/partials/jet" "$prod/results/jet"` です。

外部 root への書き込み許可は map stage のみに適用します。reduce・merge の出力や、未登録の外部リンク先は許可しません。
入力検証でも stage の出力先制限と従来の software/production/hash 照合を行うため、別 production の map は混在できません。
`--map-root` を使わない既存 production の配置規則は変わりません。

## 出力の履歴と入力の混在防止

主解析の shell wrapper は、出力 ROOT ごとに `<filename>.root.execution.json` を作ります。
記録内容は stage・引数・出力 SHA-256・software manifest SHA-256・source 版・実際の runtime release・production manifest SHA-256 です。
map の一時名は正式出力名と分けて扱い、出力名の公開時に execution record も移動します。

正式 production の reduce / merge は、入力の execution record と ROOT hash を確認し、異なる software / production の入力、改変された入力、履歴がない入力を拒否します。
旧 map を使う過去結果の再現は `run_reference.sh` 等の専用入口を使います。旧結果を新しい production として自動登録しません。

## 検証

```bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
PHOTON_INSTALL_PREFIX="$PWD/PhotonAnalysisTree/install/current" \
  bash PhotonAnalysisTree/tests/run_candidate_selection_tests.sh \
  "$PWD/PhotonAnalysisTree/build/dev/ana.565-ppg15v1/photon"
python3 PhotonAnalysisTree/tests/test_software.py
```

後者は rolling release の拒否、artifact 改変・欠落、clean-source gate、既存結果の保護、production のソフトウェア切替拒否、stage 間の入力混在を検証します。

## ビルド管理を導入したときの検証

- ana.565 で両ローカル library を新規ビルドし、既存 candidate / tagging mass テストと provenance の回帰テストを確認した。
- 実 DST 1 segment・2イベントの map → validator → reduce が成功し、従来の install による結果と histogram 2,202個と tree 内容が一致した。出力先とランダムな一時ファイル名だけを正規化した。
- 一時 clone に clean な版を作って freeze と production 初期化を確認した。元の開発 install を隠し、作業ツリーの reducer / merger を読み込めなくした状態でも保存したコードで実行できた。本作業の未コミット状態を正式 software として発行してはいない。
- 保存版で4 sample を reduce → merge し、全26 ROOT に execution record が揃った。従来結果と histogram 12,654個の全 bin 内容・誤差・軸、44 tree の内容、CSV 15ファイルが一致した。PDF 246ファイルは配置とページ数を確認した。
- 全28 Condor submit の dry-run と、過去 study の12再現設定の入力確認が通った。ジョブは投入していない。
- 既存 install 直下の両 library の SHA-256 は変化していない。過去 build、保存成果物、外部 map は保持した。

## Common composition の開発版

現在の開発版 map は schema 11 / common composition 1 / pi0 topology 14 / eta topology 4 です。現行 reducer は schema 10 map を入力として partial schema 13 を出力します。旧 schema 8・9 map や異なる partial の再ラベル、reduce のみの移行はできません。新しい DST map が必要です。固定 software・production manifest の作成は上記の clean-source 手順に従い、保存済み production と異なる出力先を指定してください。新しいコードを検証する開発 build に既存 cache を流用できない場合は、`--build-dir` と `--install-prefix` に独立した検証先を指定します。
