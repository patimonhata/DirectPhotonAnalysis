#!/usr/bin/env python3
"""Project final photon purity with perfect prompt TightBDT efficiency after preselection.

Run with the sPHENIX Python environment (PyROOT). The final background yields
are fixed; recovered prompt candidates have the nominal prompt isolation/veto
acceptance separately in each ET bin. This is a hypothetical central-value
projection, not a trained-model prediction or a statistical uncertainty estimate.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
BASE = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = BASE / 'output/plots/photon_candidate_selection/ClusterE0p12_Pi0Partner0p5/jet/candidate_composition.root'
CATEGORIES = ('prompt', 'pi0_separated', 'pi0_merged', 'pi0_single_contaminated', 'pi0_missing', 'pi0_other', 'eta', 'other')
SELECTIONS = ('preselection', 'preselection_tight', 'region_a', 'region_a_tagging_veto')


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9)


def histogram(source, selection, category):
    name = f'{selection}/h_candidate_{category}_et_pb'
    hist = source.Get(name)
    if not hist or not hist.InheritsFrom('TH1'):
        raise ValueError(f'Missing yield histogram: {name}')
    return hist


def calculate(source, low, high):
    histograms = {s: {c: histogram(source, s, c) for c in ('denominator',) + CATEGORIES} for s in SELECTIONS}
    reference = histograms['preselection']['prompt']
    edges = [reference.GetXaxis().GetBinLowEdge(i) for i in range(1, reference.GetNbinsX() + 2)]
    if not any(close(low, x) for x in edges) or not any(close(high, x) for x in edges):
        raise ValueError('ET limits must coincide with histogram edges; partial bins are not supported.')
    for selection in SELECTIONS:
        for h in histograms[selection].values():
            h_edges = [h.GetXaxis().GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 2)]
            if len(h_edges) != len(edges) or any(not close(a, b) for a, b in zip(edges, h_edges)):
                raise ValueError('Yield histograms have inconsistent ET binning.')
    rows = []
    for i in range(1, reference.GetNbinsX() + 1):
        left, right = edges[i - 1:i + 1]
        if left < low - 1e-9 or right > high + 1e-9:
            continue
        yields = {s: {c: h.GetBinContent(i) for c, h in histograms[s].items()} for s in SELECTIONS}
        for selection, values in yields.items():
            if any(not math.isfinite(x) or x < 0 for x in values.values()):
                raise ValueError(f'Invalid yield in {selection}, ET bin {i}')
            if not close(sum(values[c] for c in CATEGORIES), values['denominator']):
                raise ValueError(f'Category partition fails in {selection}, ET bin {i}')
        pre, tight, isolated, final = [yields[s]['prompt'] for s in SELECTIONS]
        if not pre >= tight >= isolated >= final or tight <= 0:
            raise ValueError(f'Prompt selection chain invalid or TightBDT yield zero in ET bin {i}')
        total = yields['region_a_tagging_veto']['denominator']
        background = total - final
        efficiency = tight / pre
        downstream = final / tight
        projected = pre * downstream
        nominal_purity = final / total
        stored_purity = source.Get('region_a_tagging_veto/h_photon_candidate_purity')
        stored_efficiency = source.Get('preselection_tight/h_candidate_prompt_survival_fraction_relative_to_previous_cut')
        if not stored_purity or not close(nominal_purity, stored_purity.GetBinContent(i)):
            raise ValueError(f'Nominal purity differs from stored purity in ET bin {i}')
        if stored_efficiency:
            if not close(efficiency, stored_efficiency.GetBinContent(i)):
                raise ValueError(f'TightBDT efficiency differs from stored previous-cut survival in ET bin {i}')
        else:
            # Schema 7 stores survival relative to kinematics, but not the previous cut.
            kinematic = histogram(source, 'kinematic', 'prompt').GetBinContent(i)
            stored_survival = source.Get('preselection_tight/h_candidate_prompt_survival_fraction_relative_to_kinematic')
            if kinematic <= 0 or not stored_survival or not close(tight / kinematic, stored_survival.GetBinContent(i)):
                raise ValueError(f'TightBDT yield differs from stored kinematic survival in ET bin {i}')
        new_purity = projected / (projected + background)
        row = dict(et_low_GeV=left, et_high_GeV=right, prompt_preselection_pb=pre, prompt_tight_pb=tight,
                   prompt_region_a_pb=isolated, prompt_final_pb=final, background_final_pb=background,
                   tight_efficiency=efficiency, isolation_acceptance=isolated / tight,
                   veto_acceptance=final / isolated if isolated else None, downstream_acceptance=downstream,
                   prompt_final_projected_pb=projected, prompt_yield_gain_percent=100 * (projected / final - 1) if final else None,
                   purity_nominal=nominal_purity, purity_projected=new_purity,
                   purity_gain_percentage_points=100 * (new_purity - nominal_purity))
        row.update({f'{c}_final_pb': yields['region_a_tagging_veto'][c] for c in CATEGORIES if c != 'prompt'})
        rows.append(row)
    if not rows:
        raise ValueError('No bins in the requested range.')
    sums = {key: sum(row[key] for row in rows) for key in rows[0] if key.endswith('_pb')}
    signal, projected, background = (sums[k] for k in ('prompt_final_pb', 'prompt_final_projected_pb', 'background_final_pb'))
    nominal = signal / (signal + background)
    improved = projected / (projected + background)
    summary = dict(et_low_GeV=low, et_high_GeV=high, **sums,
                   tight_efficiency=sums['prompt_tight_pb'] / sums['prompt_preselection_pb'],
                   prompt_yield_gain_percent=100 * (projected / signal - 1), purity_nominal=nominal,
                   purity_projected=improved, purity_gain_percentage_points=100 * (improved - nominal))
    return rows, summary


def plot(rows, output, settings):
    canvas = ROOT.TCanvas('purity_projection', 'Final photon purity projection', 900, 650)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.27)
    frame = canvas.DrawFrame(rows[0]['et_low_GeV'], 0, rows[-1]['et_high_GeV'], 100)
    frame.SetTitle(';Candidate cluster E_{T} [GeV];Prompt photon purity [%]')
    graphs = []
    legend = ROOT.TLegend(0.16, 0.86, 0.82, 0.94)
    legend.SetBorderSize(0)
    for field, color, style, label in [('purity_nominal', ROOT.kBlack, 20, 'Current TightBDT'),
                                       ('purity_projected', ROOT.kRed + 1, 24, '100% prompt TightBDT efficiency after preselection')]:
        graph = ROOT.TGraph(len(rows))
        for i, row in enumerate(rows):
            graph.SetPoint(i, (row['et_low_GeV'] + row['et_high_GeV']) / 2, 100 * row[field])
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(style)
        graph.SetMarkerSize(1.3)
        graph.SetLineWidth(2)
        graph.Draw('LP SAME')
        legend.AddEntry(graph, label, 'lp')
        graphs.append(graph)
    legend.Draw()
    caption = ROOT.TLatex()
    caption.SetNDC(True)
    caption.SetTextSize(0.027)
    caption.DrawLatex(0.16, 0.82, 'Fixed background; nominal prompt isolation/veto acceptance in each E_{T} bin')
    caption.DrawLatex(0.16, 0.77, f"#pi^{{0}}: E_{{partner}}>{settings['pi0_partner_min_energy']:g} GeV, 0.10<m<0.20 GeV; "
                      f"#eta: E_{{partner}}>{settings['eta_partner_min_energy']:g} GeV, 0.45<m<0.65 GeV")
    canvas.SaveAs(str(output / 'final_purity_projection.pdf'))
    canvas.SaveAs(str(output / 'final_purity_projection.png'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--pi0-partner-threshold', type=float, default=0.5, help='Expected input pi0 partner energy threshold [GeV]')
    parser.add_argument('--eta-partner-threshold', type=float, default=0.5, help='Expected input eta partner energy threshold [GeV]')
    parser.add_argument('--et-min', type=float, default=5)
    parser.add_argument('--et-max', type=float, default=10)
    args = parser.parse_args()
    if args.et_min >= args.et_max:
        parser.error('--et-min must be below --et-max')
    input_path = args.input.resolve()
    output = args.output_dir or input_path.parent / 'prompt_efficiency_projection'
    source = ROOT.TFile.Open(str(input_path), 'READ')
    if not source or source.IsZombie():
        raise ValueError(f'Cannot open {input_path}')
    metadata = source.Get('metadata')
    if not metadata or metadata.GetEntries() != 1:
        raise ValueError('Expected one metadata entry.')
    metadata.GetEntry(0)
    settings = {name: float(getattr(metadata, name)) for name in (
        'min_cluster_energy', 'pi0_partner_min_energy', 'eta_partner_min_energy',
        'pi0_mass_min', 'pi0_mass_max', 'eta_mass_min', 'eta_mass_max')}
    expected = dict(min_cluster_energy=0.12, pi0_partner_min_energy=args.pi0_partner_threshold, eta_partner_min_energy=args.eta_partner_threshold,
                    pi0_mass_min=0.10, pi0_mass_max=0.20, eta_mass_min=0.45, eta_mass_max=0.65)
    if str(metadata.family) != 'jet' or any(not close(settings[k], v) for k, v in expected.items()):
        raise ValueError('Input does not match the requested jet sample/tagging configuration.')
    rows, summary = calculate(source, args.et_min, args.et_max)
    output.mkdir(parents=True, exist_ok=True)
    result = dict(input=str(input_path), input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
                  settings=settings, family=str(metadata.family), weight_definition=str(metadata.weight_definition),
                  assumptions=['Prompt TightBDT efficiency relative to preselection is set to 1 in the selected ET range.',
                               'Final background yields in every category are fixed.',
                               'Recovered prompt isolation and veto acceptance equals nominal TightBDT-selected prompt acceptance in each ET bin.',
                               'Outside the selected ET range all yields are unchanged.',
                               'Central values only; statistical and scenario uncertainties are not evaluated.',
                               'Integrated purity uses sums of bin yields, not an average of bin purities.'],
                  formula='S_new[i] = S_final[i] * S_preselection[i] / S_tight[i]; P_new[i] = S_new[i] / (S_new[i] + B_final[i])',
                  bins=rows, integrated=summary)
    (output / 'projection.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    with (output / 'projection.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    lines = ['Final prompt photon purity: perfect TightBDT prompt efficiency relative to preselection',
             f'Input: {input_path}',
             f"Partner energy thresholds: pi0 = {settings['pi0_partner_min_energy']:g} GeV, eta = {settings['eta_partner_min_energy']:g} GeV",
             'Mass windows: pi0 = 0.10-0.20 GeV, eta = 0.45-0.65 GeV',
             'Yields are weighted cross sections in pb, not raw event counts.',
             result['formula'], '',
             'ET [GeV] | Tight efficiency [%] | Current purity [%] | Projected purity [%] | Gain [percentage points]']
    for row in rows + [summary]:
        lines.append(f"{row['et_low_GeV']:g}-{row['et_high_GeV']:g} | {100 * row['tight_efficiency']:.3f} | "
                     f"{100 * row['purity_nominal']:.3f} | {100 * row['purity_projected']:.3f} | {row['purity_gain_percentage_points']:.3f}")
    lines.extend(['', f"Integrated final prompt: {summary['prompt_final_pb']:.9g} -> {summary['prompt_final_projected_pb']:.9g} pb",
                  f"Integrated final background (fixed): {summary['background_final_pb']:.9g} pb",
                  f"Integrated prompt yield gain: {summary['prompt_yield_gain_percent']:.3f}%", '', 'Assumptions:'] + result['assumptions'])
    report = '\n'.join(lines) + '\n'
    (output / 'report.txt').write_text(report)
    plot(rows, output, settings)
    source.Close()
    print(report)
    print(f'Outputs: {output.resolve()}')


if __name__ == '__main__':
    main()
