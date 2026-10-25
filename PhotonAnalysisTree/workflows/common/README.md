# Shared workflow entry points

- `scoring/`: single_particle / pythia_tree 共通の score 追加マクロ。
- `make_pythia_input_manifest.sh`: 同期した4 stream list の suffix manifest 検証・生成。
- `replay_study.py`: `workflows/studies/<study>/reproduction.json` に保存されたパラメータによる再実行。

描画スタイルの内部実装は `analysis/common/plotting/` に置きます。
物理的な分類・matching は `src/` の共有実装を使います。

`runtime.sh` は固定 release と開発/保存済み software の選択、source snapshot への dispatch を担当します。[設計と手順](../../docs/software.md)を参照してください。
