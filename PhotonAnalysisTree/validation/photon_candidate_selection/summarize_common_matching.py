#!/usr/bin/env python3
"""Summarize bounded DST diagnostics and verify them against production map_000000.

The input diagnostics must use manifest row 0. Counterfactual gates only count
saved labels: they are NOT validated physics classifications or production fixes.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import statistics
import ROOT


def flag(x):
    return bool(ord(x)) if isinstance(x, (str, bytes)) else bool(x)


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def summarize(sample, directory, production):
    inputs = dict(line.split('=', 1) for line in (directory / 'input.txt').read_text().splitlines())
    if inputs['row'] != '0': raise ValueError('This checker requires manifest row 0')
    events = read(directory / 'events.csv')
    candidates = read(directory / 'candidates.csv')
    anchors = read(directory / 'anchors.csv')
    failures = read(directory / 'first_failures.csv')
    expected_events = int(inputs['requested_events'])
    if {int(r['event_in_file']) for r in events} != set(range(expected_events)): raise RuntimeError('Incomplete DST run')
    candidate_by_key = {(int(r['event_in_file']), int(r['pdg']), int(r['candidate_index'])): r for r in candidates}
    anchor_by_key = {(int(r['event_in_file']), int(r['pdg']), int(r['cluster_id'])): r for r in anchors}
    checked_candidates, checked_anchors = set(), set()
    path = production / f'maps/{sample}/map_000000.root'
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie(): raise RuntimeError(f'Cannot open {path}')
    meta = f.Get('metadata'); meta.GetEntry(0)
    if str(meta.first_input_suffix) != inputs['suffix']: raise RuntimeError('DST/map input mismatch')
    t = f.Get('event_tree'); t.SetBranchStatus('*', 0)
    for k in ('source_file_id', 'event_in_file', 'hepmc_event_number'):
        t.SetBranchStatus(k, 1)
    for species in ('pi0', 'eta'):
        for kind in ('anchor_cluster_id', 'anchor_candidate_index', 'anchor_topology', 'candidate_evaluated', 'candidate_matching_usable', 'candidate_parent_id'):
            t.SetBranchStatus(species + '_topology_' + kind, 1)
    for entry in range(t.GetEntries()):
        if t.GetEntry(entry) <= 0: raise RuntimeError('Unreadable event')
        event = int(t.event_in_file)
        if int(t.source_file_id) != 0 or event >= expected_events: break
        for species, pdg in [('pi0', 111), ('eta', 221)]:
            get = lambda k: list(getattr(t, species + '_topology_' + k))
            for ci, evaluated in enumerate(get('candidate_evaluated')):
                if not flag(evaluated): continue
                key = event, pdg, ci
                row = candidate_by_key[key]
                if int(row['hepmc_event_number']) != int(t.hepmc_event_number): raise RuntimeError('Event number mismatch')
                if flag(get('candidate_matching_usable')[ci]) != bool(int(row['global_usable'])) or int(get('candidate_parent_id')[ci]) != int(row['parent_id']):
                    raise RuntimeError(f'Candidate mismatch: {key}')
                checked_candidates.add(key)
            for cluster, ci, topology in zip(get('anchor_cluster_id'), get('anchor_candidate_index'), get('anchor_topology')):
                key = event, pdg, int(cluster)
                row = anchor_by_key[key]
                if int(ci) != int(row['candidate_index']) or int(topology) != int(row['topology']): raise RuntimeError(f'Anchor mismatch: {key}')
                checked_anchors.add(key)
    f.Close()
    if checked_candidates != set(candidate_by_key) or checked_anchors != set(anchor_by_key): raise RuntimeError('Incomplete map/DST agreement')
    gates = Counter()
    for r in anchors:
        c = candidate_by_key[int(r['event_in_file']), int(r['pdg']), int(r['candidate_index'])]
        prerequisites = int(c['projection_ok']) and int(c['interaction_ok'])
        gates['current_applicable'] += int(r['applicability']) == 1
        gates['omit_below_threshold_only'] += bool(prerequisites and int(c['above_failures']) == 0)
        gates['anchor_only_diagnostic'] += bool(prerequisites and int(r['anchor_match_usable']))
        gates['anchor_and_two_representatives_diagnostic'] += bool(prerequisites and int(r['anchor_match_usable']) and int(c['representatives_usable']))
    failure_groups = Counter((r['stage'], r['failure'], r['ancestry_complete'], r['is_anchor'], r['is_representative']) for r in failures)
    summary = dict(sample=sample, requested_events=expected_events, accepted_events=sum(r['pdg'] == '111' and r['status'] == '2' for r in events),
                   evaluated_candidates=len(candidates), anchors=len(anchors), production_map_candidates_and_topologies_identical=True,
                   usable_candidates=sum(int(r['global_usable']) for r in candidates),
                   usable_anchor_matches=sum(int(r['anchor_match_usable']) for r in anchors),
                   all_candidate_anchors_usable=sum(int(r['all_anchors_usable']) for r in candidates),
                   median_clusters={k: statistics.median(float(r[k]) for r in candidates) if candidates else None
                                    for k in ('above_count', 'above_failures', 'below_count', 'below_failures')},
                   first_failures=[dict(stage=k[0], failure=k[1], ancestry_complete=int(k[2]), is_anchor=int(k[3]), is_representative=int(k[4]), count=v)
                                   for k,v in sorted(failure_groups.items())],
                   first_failure_min_delta_r=min((float(r['min_photon_delta_r']) for r in failures), default=None),
                   counterfactual_gate_counts=dict(gates),
                   saved_topology_counts=dict(Counter(r['topology'] for r in anchors)))
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('production', type=Path); p.add_argument('output', type=Path)
    p.add_argument('diagnostics', nargs='+', help='SAMPLE=DIRECTORY')
    args = p.parse_args()
    if args.output.exists(): p.error('Output already exists')
    summaries = []
    for item in args.diagnostics:
        sample, directory = item.split('=', 1)
        summaries.append(summarize(sample, Path(directory), args.production))
    result = dict(scope='Bounded first-file DST events; no stitching/kinematic selection. Gate variants are diagnostic only.', samples=summaries)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
