# Eta decay flight diagnostic

`CheckEtaDecayFlight.C` は G4Hits DST を読み、eta の飛行距離・崩壊を独立に調べるマクロです。
主解析の composition / topology 選別を変更しません。引数・物理的な解釈・使用例は
[定義ノート](../../../docs/photon_candidate_definitions.md) の `CheckEtaDecayFlight` の説明を参照してください。
既定の出力先は引き続き `output/photon_candidate_selection/qa/eta_decay_flight/` です。
自動テストは `tests/test_eta_decay_flight.C` にあります。
