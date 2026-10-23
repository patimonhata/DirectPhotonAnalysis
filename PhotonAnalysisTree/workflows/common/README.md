# Shared workflow helpers

- `plotting/`: sPHENIX 描画スタイル。
- `scoring/`: single_particle / pythia_tree 共通の score 追加マクロ。
- `make_pythia_input_manifest.sh`: 同期した4 stream list の suffix manifest 検証・生成。
- `replay_study.py`: 過去 study の保存パラメータによる再実行。

物理的な分類・matching は `src/` の共有実装を使います。
