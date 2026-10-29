# 集計・描画の内部実装

ここは実行用マクロから利用する ROOT の集計・描画部品です。通常の解析操作は [workflows](../workflows/README.md) から行います。
`src/` の Fun4All module、truth matching、topology evaluator は引き続き共有の物理処理を所有します。

## Photon candidate selection

| 部品 | 責務 |
|---|---|
| `CandidateSelectionCommon.h` | 選別・sample 定義、map metadata と整合性検証、正規化の入力、基本 histogram |
| `CandidateSelectionPlots.h` | 共通の描画・caption・軸・legend |
| `CandidateComposition.h` | 排他的な起源分類、composition histogram、fraction と整合性検証 |
| `CandidateSelection.h` | reducer が利用する集計部品の入口 |
| `CandidateMerge.h` | partial metadata、shard 検証、統合・survival の共通処理 |
| `*Histograms.h` / `*Plots.h` | pi0/eta topology、eta decay、Other/decay-origin、BDT の集計と描画 |
| `ReduceCandidateSelection.C` / `MergeCandidateSelection.C` | reduce / merge の実装 |

実行用 `ReducePythiaPhotonCandidateSelection.C` と `MergePythiaPhotonCandidateSelection.C` はそれぞれの実装へ委譲します。
Merge は Reduce の実行マクロや実装を読み込まず、両方が共通ヘッダを参照します。
ROOT interpreter による従来の実行方式を維持しており、このディレクトリを別ライブラリとして install する手順はありません。
共通ヘッダは include guard で保護し、依存する内部ヘッダを明示しています。
通常の shell wrapper は各 stage を別の ROOT process で実行します。synthetic integration test は、独立した両入口を
一つの一時マクロに include して読み込み、ROOT の macro autoload と従来の単一翻訳単位での実行を維持しています。

## 共通・補助解析

- `calorimeter_acceptance/`: ana.565 tower geometry の抽出と解析的な z-vertex acceptance の描画。入口は `workflows/studies/calorimeter_acceptance/`。
- `common/PhotonCandidateSettings.h`: map の条件読み取りと partial 間の条件検証。複数の診断で共有します。
- `common/plotting/`: sPHENIX 描画スタイル。
- `tagging_mass/` / `tagging_working_point_scan/`: 補助解析の histogram・metadata・描画部品。

主解析が補助 workflow のコードを読み込む依存は作りません。共有する条件・表現は `common/`、共有する物理判定は `src/` に置きます。
