"""Write the weighted comparison and raw pi0-cluster counts from a v2 summary."""
import argparse
import csv
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('output', type=Path)
a = p.parse_args()
rows = list(csv.DictReader((a.output / 'summary.tsv').open(), delimiter='\t'))
order = ['e0p5_m0p10_0p20', 'e0p15_m0p10_0p20', 'e0p15_m0p08_0p18', 'e0p15_m0p08_0p20', 'e0p12_m0p08_0p18', 'e0p12_m0p08_0p20']
labels = ['0.50 / 0.10–0.20', '0.15 / 0.10–0.20', '0.15 / 0.08–0.18', '0.15 / 0.08–0.20', '0.12 / 0.08–0.18', '0.12 / 0.08–0.20']
lines = ['# π⁰ tagging 比較（基準: E > 0.5 GeV, 0.10 < m < 0.20 GeV）', '',
         '同じ topology-v11 map の Region A を全設定で使用。η veto は E > 0.5 GeV、0.45 < m < 0.65 GeV に固定。', '',
         'Purity と π⁰ 割合は sample weight を適用した値。π⁰ は既存 composition と同じ排他的 majority 分類（anchor valid、main fraction > 0.5）。π⁰ 数は生成 π⁰ 粒子数ではなく、選別後の π⁰ 由来 candidate cluster 数。生カウントは異なる Jet sample の生成量の影響を受けるため、設定間の物理的比較には重み付き収量 [pb] を使う。', '',
         '追加 prompt 損失（正味）= 1 − N_prompt(setting) / N_prompt(reference)。質量窓の変更により基準で落ちた prompt が戻ることもあるため、追加除去率と回復率も別表に示す。すべて基準で残った prompt の重み付き収量を分母とする。', '',
         'Purity と prompt 生存率の誤差は summary.tsv に保存。設定間には相関があり、ここでは有意差を主張しない。', '']
for lo, hi in [('5','6'),('6','8'),('8','10'),('10','18'),('18','30'),('5','30')]:
    group = {r['working_point']: r for r in rows if r['et_low'] == lo and r['et_high'] == hi}
    lines += [f'## Candidate ET {lo}–{hi} GeV', '', '| Partner 閾値 / mass window [GeV] | Purity | 追加 prompt 損失（正味） | π⁰ 生カウント | π⁰ 収量 [pb] | π⁰ 割合 |', '|---|---:|---:|---:|---:|---:|']
    for key, label in zip(order, labels):
        r = group[key]
        lines.append(f"| {label} | {100*float(r['purity']):.2f}% | {100*float(r['prompt_net_loss_fraction']):.2f}% | {float(r['pi0_after_count']):,.0f} | {float(r['pi0_after_pb']):,.2f} | {100*float(r['pi0_fraction']):.2f}% |")
    lines += ['', '| 設定 [GeV] | 追加 prompt 除去率 | prompt 回復率 |', '|---|---:|---:|']
    for key, label in zip(order, labels):
        r=group[key]
        lines.append(f"| {label} | {100*float(r['prompt_gross_loss_fraction']):.2f}% | {100*float(r['prompt_recovery_fraction']):.2f}% |")
    lines += ['']
(a.output / 'RESULTS.md').write_text('\n'.join(lines))
print(a.output / 'RESULTS.md')
