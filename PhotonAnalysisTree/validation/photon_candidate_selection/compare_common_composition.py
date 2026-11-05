#!/usr/bin/env python3
"""Compare a legacy map with a common-composition map made from identical input/cuts/model."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import ROOT


def equal(a, b):
    if isinstance(a, (float, int)) and isinstance(b, (float, int)):
        return a == b or (math.isnan(float(a)) and math.isnan(float(b)))
    return a == b


def flag(value):
    return bool(ord(value)) if isinstance(value, (str, bytes)) else bool(value)


def compare(old_path, new_path):
    old_file, new_file = (ROOT.TFile.Open(str(p)) for p in (old_path, new_path))
    assert old_file and new_file and not old_file.IsZombie() and not new_file.IsZombie()
    old, new = old_file.Get('event_tree'), new_file.Get('event_tree')
    assert old and new and old.GetEntries() == new.GetEntries()
    old_meta, new_meta = old_file.Get('metadata'), new_file.Get('metadata')
    assert old_meta.GetEntry(0) > 0 and new_meta.GetEntry(0) > 0
    for key in ('model_sha256', 'min_cluster_energy', 'pi0_partner_min_energy', 'eta_partner_min_energy',
                'pi0_mass_min', 'pi0_mass_max', 'eta_mass_min', 'eta_mass_max', 'missing_energy_min', 'missing_energy_max'):
        assert getattr(old_meta, key) == getattr(new_meta, key), key
    assert int(old_meta.schema_version) == 8 and int(new_meta.schema_version) in (9, 10)
    keys = [b.GetName() for b in old.GetListOfBranches() if b.GetName().startswith((
        'split_cluster_pass_', 'split_cluster_bdt_', 'split_cluster_iso_'))]
    keys += ['split_cluster_id', 'split_cluster_e', 'split_cluster_et', 'split_cluster_pi0_tag', 'split_cluster_eta_tag',
             'split_cluster_pi0_partner_cluster_id', 'split_cluster_eta_partner_cluster_id',
             'split_cluster_pi0_partner_mass', 'split_cluster_eta_partner_mass']
    clusters = 0
    transitions, reasons = Counter(), Counter()
    for entry in range(old.GetEntries()):
        assert old.GetEntry(entry) > 0 and new.GetEntry(entry) > 0
        assert old.event_uid == new.event_uid and old.split_ncluster == new.split_ncluster
        for key in keys:
            left, right = list(getattr(old, key)), list(getattr(new, key))
            assert len(left) == len(right) and all(equal(a, b) for a, b in zip(left, right)), (entry, key)
        for c in range(new.split_ncluster):
            prompt = flag(old.split_cluster_truth_prompt_cluster[c]) and old.split_cluster_truth_dominant_fraction[c] > .5
            pi0 = flag(old.split_cluster_pi0_anchor_valid[c]) and old.split_cluster_pi0_anchor_main_fraction[c] > .5
            eta = not pi0 and old.split_cluster_truth_eta_main_fraction[c] > .5
            origin = (0 if prompt else 1 if pi0 else 2) if sum((prompt, pi0, eta)) == 1 else 3
            assert new.split_cluster_common_legacy_origin[c] == origin, (entry, c, 'legacy origin')
            assert equal(new.original_cluster_truth_total_edep[c], new.split_cluster_common_total[c])
            total, attributed = new.split_cluster_common_total[c], new.split_cluster_common_attributed[c]
            assert math.isclose(total, attributed, rel_tol=1e-5, abs_tol=1e-12) or new.split_cluster_common_quality[c]
            main = int(new.split_cluster_common_main[c])
            transitions[f'{origin}->{main}'] += 1
            reasons[int(new.split_cluster_common_topology_reason[c])] += 1
        clusters += new.split_ncluster
    return {'events': old.GetEntries(), 'clusters': clusters, 'identical_reconstruction_branches': len(keys),
            'model_sha256': str(new_meta.model_sha256), 'all_origin_transitions': dict(transitions), 'topology_reasons': dict(reasons)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('old_map', type=Path)
    parser.add_argument('new_map', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)
    result = json.dumps(compare(args.old_map, args.new_map), indent=2)
    print(result)
    if args.output:
        args.output.write_text(result + '\n')


if __name__ == '__main__':
    main()
