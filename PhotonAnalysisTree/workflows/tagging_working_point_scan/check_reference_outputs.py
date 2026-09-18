"""Compare a v2 scan to existing photon compositions and the v1 scan."""
import argparse
import math
from pathlib import Path
import ROOT

p = argparse.ArgumentParser()
p.add_argument('scan', type=Path)
p.add_argument('composition_base', type=Path)
p.add_argument('previous_scan', type=Path)
a = p.parse_args()
ROOT.gROOT.SetBatch(True)

def open_root(path):
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f'Cannot open {path}')
    return f

scan = open_root(a.scan)
checks = 0

def compare(x, y, context):
    global checks
    if not math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-7):
        raise AssertionError((context, x, y))
    checks += 1

for config, bit in [('ClusterE0p12_Pi0Partner0p5_v11', 4), ('ClusterE0p12_Pi0Partner0p15_v11', 5), ('ClusterE0p12_Pi0Partner0p15_v11_Pi0Window0p08-0p18', 0)]:
    reference = open_root(a.composition_base / config / 'jet/candidate_composition.root')
    for kind in ['count', 'pb']:
        for selection, selected in [('region_a', False), ('region_a_tagging_veto', True)]:
            for category, origins, components in [
                ('prompt', ['prompt'], ['prompt']),
                ('all', ['prompt', 'background'], ['denominator']),
                ('pi0', ['pi0'], ['pi0_' + k for k in ['separated', 'merged', 'single_contaminated', 'missing', 'other']])]:
                actual = [scan.Get(origin + '_survival_mask_' + kind) for origin in origins]
                expected = [reference.Get(selection + '/h_candidate_' + c + '_et_' + kind) for c in components]
                for h in expected:
                    assert h and h.GetNbinsX() == actual[0].GetNbinsX()
                    for edge in range(1, h.GetNbinsX() + 2):
                        assert h.GetXaxis().GetBinLowEdge(edge) == actual[0].GetXaxis().GetBinLowEdge(edge)
                masks = [m + 1 for m in range(64) if not selected or m & (1 << bit)]
                for x in range(actual[0].GetNbinsX() + 2):
                    for variance in [False, True]:
                        av = sum(h.GetBinError(x, y)**2 if variance else h.GetBinContent(x, y) for h in actual for y in masks)
                        ev = sum(h.GetBinError(x)**2 if variance else h.GetBinContent(x) for h in expected)
                        compare(av, ev, (config, kind, selection, category, x, variance))
    reference.Close()
previous = open_root(a.previous_scan)
for origin in ['prompt', 'background']:
    for kind in ['count', 'pb']:
        name = origin + '_survival_mask_' + kind
        actual, expected = scan.Get(name), previous.Get(name)
        for x in range(actual.GetNbinsX() + 2):
            for mask in range(16):
                bins = [m + 1 for m in range(64) if (m & 15) == mask]
                compare(sum(actual.GetBinContent(x, y) for y in bins), expected.GetBinContent(x, mask + 1), (name, x, mask))
                compare(sum(actual.GetBinError(x, y)**2 for y in bins), expected.GetBinError(x, mask + 1)**2, (name, x, mask, 'variance'))
print(f'Existing three compositions and v1 scan agree: {checks} comparisons')
