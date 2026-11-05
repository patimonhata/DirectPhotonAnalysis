#!/usr/bin/env python3
"""Check bounded schema-9 -> schema-10 maps against the independent DST matching audit.

Compare every pre-existing event branch except applicability/reason, using source
file/event identity. The DST audit supplies the independent per-anchor expectation.
Neither input nor existing output is overwritten.
"""
import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import ROOT


def plain(value):
    if isinstance(value, (int, bool, str, bytes)):
        return value
    if isinstance(value, float):
        return 'NaN' if math.isnan(value) else value
    return tuple(plain(x) for x in value)


def flag(value):
    return bool(ord(value)) if isinstance(value, (str, bytes)) else bool(value)


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def compare(old_path, new_path, diagnostic):
    old_file, new_file = (ROOT.TFile.Open(str(p)) for p in (old_path, new_path))
    if any(not f or f.IsZombie() for f in (old_file, new_file)):
        raise RuntimeError('Cannot read input map')
    old_meta, new_meta = (f.Get('metadata') for f in (old_file, new_file))
    old_meta.GetEntry(0); new_meta.GetEntry(0)
    assert (old_meta.schema_version, new_meta.schema_version) == (9, 10)
    assert (old_meta.pi0_topology_algorithm_version, new_meta.pi0_topology_algorithm_version) == (13, 14)
    assert (old_meta.eta_topology_algorithm_version, new_meta.eta_topology_algorithm_version) == (3, 4)
    assert str(old_meta.first_input_suffix) == str(new_meta.first_input_suffix)
    for field in ('sample_name', 'common_composition_algorithm_version', 'model_sha256', 'min_cluster_energy', 'pi0_partner_min_energy',
                  'eta_partner_min_energy', 'pi0_mass_min', 'pi0_mass_max', 'eta_mass_min', 'eta_mass_max',
                  'missing_energy_min', 'missing_energy_max', 'min_photon_energy_recovery'):
        assert str(getattr(old_meta, field)) == str(getattr(new_meta, field)), field
    candidates = {(int(r['event_in_file']), int(r['pdg']), int(r['candidate_index'])): r for r in rows(diagnostic / 'candidates.csv')}
    anchors = {(int(r['event_in_file']), int(r['pdg']), int(r['cluster_id'])): r for r in rows(diagnostic / 'anchors.csv')}
    old, new = (f.Get('event_tree') for f in (old_file, new_file))
    prefix = 'split_cluster_common_'
    changed = {prefix + 'applicability', prefix + 'topology_reason'}
    branches = {b.GetName() for b in old.GetListOfBranches()}
    new_branches = {b.GetName() for b in new.GetListOfBranches()}
    assert new_branches - branches == {prefix + 'anchor_matching_' + k for k in ('status', 'coverage', 'usable')}
    assert not branches - new_branches
    stable = sorted(branches - changed)
    old.SetBranchStatus('*', 0)
    for b in ('source_file_id', 'event_in_file'):
        old.SetBranchStatus(b, 1)
    lookup = {}
    for i in range(old.GetEntries()):
        old.GetEntry(i)
        key = int(old.source_file_id), int(old.event_in_file)
        assert key not in lookup
        lookup[key] = i
    old.SetBranchStatus('*', 1)
    transitions, topology, totals, checked = Counter(), Counter(), Counter(), set()
    for i in range(new.GetEntries()):
        new.GetEntry(i)
        key = int(new.source_file_id), int(new.event_in_file)
        old.GetEntry(lookup[key])
        for b in stable:
            if plain(getattr(old, b)) != plain(getattr(new, b)):
                raise AssertionError(f'Unexpected change: event {key}, branch {b}')
        for c, cluster in enumerate(new.split_cluster_id):
            main = int(new.split_cluster_common_main[c])
            app = int(new.split_cluster_common_applicability[c])
            reason = int(new.split_cluster_common_topology_reason[c])
            old_app = int(old.split_cluster_common_applicability[c])
            local_status = int(new.split_cluster_common_anchor_matching_status[c])
            if local_status < 0:
                assert app == old_app and reason == int(old.split_cluster_common_topology_reason[c])
                continue
            pdg = 111 if main == 2 else 221 if main == 3 else None
            assert pdg is not None and key[0] == 0
            anchor_key = key[1], pdg, int(cluster)
            a = anchors[anchor_key]; candidate = candidates[key[1], pdg, int(a['candidate_index'])]
            usable = bool(int(a['anchor_match_usable']))
            projection, interaction = int(candidate['projection_ok']), int(candidate['interaction_ok'])
            expected_reason = 4 if not usable else 5 if not projection else 6 if not interaction else 0
            assert app == int(expected_reason == 0) and reason == expected_reason
            assert flag(new.split_cluster_common_anchor_matching_usable[c]) == usable
            assert local_status == int(a['anchor_match_status'])
            assert math.isclose(new.split_cluster_common_anchor_matching_coverage[c], float(a['anchor_match_coverage']), abs_tol=1e-9)
            checked.add(anchor_key)
            species = 'pi0' if pdg == 111 else 'eta'
            transitions[f'{species}:{old_app}->{app}'] += 1
            for selection, passed in (('all_anchors', True), ('kinematic', flag(new.split_cluster_pass_kinematics[c])),
                                      ('region_a', flag(new.split_cluster_pass_region_a[c]))):
                if not passed or (selection != 'all_anchors' and not flag(new.sample_stitching_pass)):
                    continue
                totals[f'{species}:{selection}:applicable'] += app == 1
                totals[f'{species}:{selection}:unknown'] += app == 0
                if app:
                    topology[f'{species}:{selection}:{getattr(new, "split_cluster_" + species + "_anchor_topology")[c]}'] += 1
    assert checked == set(anchors), 'Not all independent diagnostic anchors were checked'
    result = dict(sample=str(new_meta.sample_name), events=int(new.GetEntries()), stable_branches=len(stable),
                  stable_branch_comparison='exact equality including nested vectors and NaNs',
                  checked_anchors=len(checked), applicability_transitions=dict(transitions), populations=dict(totals), topology=dict(topology))
    old_file.Close(); new_file.Close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('old_map', 'new_map', 'diagnostic', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Output already exists')
    result = compare(args.old_map, args.new_map, args.diagnostic)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
