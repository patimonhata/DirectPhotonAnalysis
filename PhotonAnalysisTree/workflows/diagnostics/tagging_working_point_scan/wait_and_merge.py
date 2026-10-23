"""Wait for this scan's own Condor jobs, then merge and validate the baseline once."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import time

p = argparse.ArgumentParser()
p.add_argument('partial_root', type=Path)
p.add_argument('output', type=Path)
p.add_argument('baseline', type=Path)
p.add_argument('--clusters', nargs='+', required=True)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
state_path = a.output / 'execution_status.json'

def state(phase, **extra):
    value = {'phase': phase, 'updated_unix': time.time(), 'clusters': a.clusters, **extra}
    temporary = state_path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n'); temporary.replace(state_path)
    print(json.dumps(value), flush=True)

expected = [a.partial_root / sample / f'shard_{shard}' / 'working_point_scan.root'
            for sample in ['jet3', 'jet5', 'jet8', 'jet12', 'jet20', 'jet30', 'jet40']
            for shard in range(10 if sample == 'jet12' else 1)]
try:
    start = time.monotonic()
    last_count = -1
    while True:
        count = sum(path.is_file() for path in expected)
        if count != last_count:
            state('reducing', completed_partials=count, expected_partials=len(expected)); last_count = count
        if count == len(expected):
            break
        query = subprocess.run(['condor_q', *a.clusters, '-json', '-attributes', 'ClusterId,ProcId,JobStatus,HoldReason'], check=True, capture_output=True, text=True)
        jobs = json.loads(query.stdout)
        failed = [j for j in jobs if j['JobStatus'] in [3, 5]]
        if failed:
            raise RuntimeError('Scan job held or removed: ' + json.dumps(failed))
        if not jobs:
            time.sleep(5)
            if not all(path.is_file() for path in expected):
                raise RuntimeError('Scan jobs left the queue without all expected outputs')
        if time.monotonic() - start > 24 * 3600:
            raise RuntimeError('24-hour monitor limit reached; jobs have not been cancelled')
        time.sleep(30)
    state('merging', completed_partials=len(expected))
    subprocess.run([str(Path(__file__).resolve().parent / 'run_merge.sh'), 'jet', str(a.partial_root), str(a.output)], check=True)
    state('validating_baseline')
    import ROOT
    ROOT.gROOT.SetBatch(True)
    actual = ROOT.TFile.Open(str(a.output / 'working_point_scan.root'))
    baseline = ROOT.TFile.Open(str(a.baseline))
    if not actual or actual.IsZombie() or not baseline or baseline.IsZombie():
        raise RuntimeError('Cannot open merged scan/baseline ROOT')
    comparisons = 0
    for kind in ['count', 'pb']:
        prompt = actual.Get('prompt_survival_mask_' + kind)
        background = actual.Get('background_survival_mask_' + kind)
        pi0 = actual.Get('pi0_survival_mask_' + kind)
        for selection, selected in [('region_a', False), ('region_a_tagging_veto', True)]:
            refp = baseline.Get(selection + '/h_candidate_prompt_et_' + kind)
            refall = baseline.Get(selection + '/h_candidate_denominator_et_' + kind)
            refpi0 = [baseline.Get(selection + '/h_candidate_pi0_' + category + '_et_' + kind) for category in ['separated', 'merged', 'single_contaminated', 'missing', 'other']]
            if refp.GetNbinsX() != prompt.GetNbinsX():
                raise RuntimeError('Baseline ET axes differ')
            for x in range(prompt.GetNbinsX() + 2):
                bins = [mask + 1 for mask in range(prompt.GetNbinsY()) if not selected or mask & (1 << 4)]
                pv = sum(prompt.GetBinContent(x, y) for y in bins)
                bv = sum(background.GetBinContent(x, y) for y in bins)
                for observed, expected_value in [(pv, refp.GetBinContent(x)), (pv + bv, refall.GetBinContent(x)), (sum(pi0.GetBinContent(x, y) for y in bins), sum(h.GetBinContent(x) for h in refpi0))]:
                    if not math.isclose(observed, expected_value, rel_tol=1e-9, abs_tol=1e-7):
                        raise RuntimeError(f'Baseline mismatch: {kind}, {selection}, ET bin {x}: {observed} vs {expected_value}')
                    comparisons += 1
    actual.Close(); baseline.Close()
    subprocess.run([os.sys.executable, str(Path(__file__).resolve().parent / 'report_scan.py'), str(a.output)], check=True)
    state('complete', baseline_comparisons= comparisons, summary=str(a.output / 'summary.tsv'))
except Exception as e:
    state('failed', error=str(e))
    raise
