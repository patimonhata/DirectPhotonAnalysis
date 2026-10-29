# Calorimeter の z-vertex 依存 acceptance

ana.565 の `TOWERGEOM_CEMC / HCALIN / HCALOUT` を geometry-only initialization で生成し、
各 tower の共通中心半径上で η bin の外側境界を vertex から見直す独立 study です。
粒子の生成・Geant4 transport・DST 入力は不要です。

## 実行

sPHENIX の **ana.565** を設定した環境で、リポジトリのルートから:

```bash
PhotonAnalysisTree/workflows/studies/calorimeter_acceptance/run.sh
```

`OFFLINE_MAIN` は末尾が `ana.565` の directory、`CALIBRATIONROOT` はその `share/calibrations` を要求します。
runner はこの release の include と library directory を先頭に置きます。
正式出力は次の配置です。

```text
PhotonAnalysisTree/output/calorimeter_acceptance/productions/ana565/tower_surface/
  manifest.json              # 定義、条件、source / library / XML / output SHA-256、検証結果
  extract.log
  plot.log
  results/
    tower_geometry.root      # 実際に生成した3つの TOWERGEOM node
    geometry.csv             # 代表半径、端点、bin数、container の半径・厚さ
    eta_bins.csv             # 全 η bin の境界
    acceptance_scan.csv      # -100～100 cm、1 cm刻み、3検出器
    acceptance_selected_vertices.csv  # 0～100 cm、10 cm刻み、3検出器
    acceptance_vs_z.pdf      # z_vertex–η の灰色 acceptance、3パネル、1ページ
    detector_comparison.pdf  # 境界と Δη の比較、1ページ
    acceptance_vertex_views.pdf  # z–r と η–φ、3検出器、11ページ
    acceptance_vs_z.png      # 表示確認用
    vertex_000.png
    vertex_100.png
    validation.json
```

既存 production の上書きは拒否します。再実行時は新しい production directory を指定してください。

```bash
PhotonAnalysisTree/workflows/studies/calorimeter_acceptance/run.sh /tmp/calo-acceptance-new-production
```

geometry の再抽出なしで保存 CSV から再描画する場合（PDF/PNG/acceptance CSV は上書き）:

```bash
root -l -b -q 'PhotonAnalysisTree/workflows/studies/calorimeter_acceptance/PlotCalorimeterAcceptance.C("/absolute/path/to/results")'
```

この直接呼び出しは既存 manifest の hash を更新しません。正式な再生成には新しい production を使用してください。

## 定義と近似

位置は cm。vertex は `(0, 0, zv)`。各検出器の全 tower について中心半径 `R` が共通であること、
`asinh(z_center/R)` が bin 中心 η に一致することを検証します。
container の `get_radius()` は内半径であり、HCAL の tower 中心半径とは異なります。

原点から見た最外側の bin 境界を `eta0_min/max` とすると、代表面上の端点は
`z_min/max = R*sinh(eta0_min/max)`、vertex から見た境界は
`eta_min/max(zv) = asinh((z_min/max-zv)/R)` です。

- 全 φ、横方向 vertex はゼロ。内部の η / φ gap、dead channel は無視した外側 envelope です。
- η の固定 cut は加えません。η–φ 図の灰色は `[-pi, pi]` 全域を覆います。
- これは tower 代表面での幾何学的 acceptance です。物理的入口面、材料全体への交差、シャワー包含、再構成効率ではありません。
- EMCAL の bin 境界は SPACAL geometry から導出されますが、代表面への投影は RawTowerBuilder の座標定義に従います。
- z–r 図の灰色 ribbon は代表面を見せるために半径方向 ±2 cm の描画幅を持たせています。この幅は計算に使いません。
- 全 vertex ページで軸範囲を固定します。正の vertex では acceptance が負の η 側へ移動します。

## Geometry の取得

実行入口は `workflows/studies/calorimeter_acceptance/`、実装は `analysis/calorimeter_acceptance/`、
独立検証は `validation/calorimeter_acceptance/` に置いています。

1. インストール済み `G4_CEmc_Spacal.C`, `G4_HcalIn_ref.C`, `G4_HcalOut_ref.C` で subsystem を登録します。
2. EMCAL は登録された subsystem の XML 読み込み後、macro override を適用します。
   `PHG4CylinderGeom_Spacalv3` → `PHG4FullProjSpacalCellReco::InitRun` → `RawTowerBuilder::InitRun` で tower geometry を生成します。
   cell builder の初期化に必要な空の hit container を用意しますが、event 処理はしません。
3. HCAL は標準 GDML subsystem の parameter を `G4GEOPARAM_*` に保存し、
   `HcalRawTowerBuilder::InitRun` で tower geometry を生成します。標準の `phistart` も builder に渡します。
   この代表面 geometry は subsystem の半径・bin 設定から構築されるので、CDB / GDML の物理 volume 読み込みは不要です。
4. tower の全数・半径・中心 η・全 bin 境界の有限性を検証して ROOT/CSV に保存します。

ana.565 での抽出結果:

| Detector | R [cm] | 原点の η 範囲 | z端点 [cm] | η × φ bins |
|---|---:|---:|---:|---:|
| EMCAL | 93.5 | ±1.1338096125 | ±130.2299453 | 96 × 256 |
| iHCAL | 127.5025 | ±1.1 | ±170.2983916 | 24 × 64 |
| oHCAL | 225.87 | ±1.1 | ±301.6826941 | 24 × 64 |

この表を計算コードに hard-code してはいません。実行ごとに installed geometry から取得します。

## 検証

runner は最後に独立した `atan2` / `-log(tan(theta/2))` 計算、境界の逆変換、
vertex に対する単調性、正負 vertex の鏡映対称性、bin の外側境界、scan 点の完全性、PDF のページ数を確認します。
PDF は 1 + 1 + 11 ページ。数値表は 603 + 33 行です。

```bash
python3 PhotonAnalysisTree/validation/calorimeter_acceptance/check_outputs.py /absolute/path/to/results
```

内部 gap の最大値も `validation.json` に記録します。これは envelope 近似で無視する量です。
