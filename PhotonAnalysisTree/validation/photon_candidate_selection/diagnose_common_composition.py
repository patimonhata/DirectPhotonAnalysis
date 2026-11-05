#!/usr/bin/env python3
"""Read-only production audit: full partial migrations plus bounded, paired map examples.

Outputs are created in a new directory. Full-sample weights come from existing
partials; the bounded map scan is diagnostic and is never extrapolated.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import ROOT

SELECTIONS = ('kinematic', 'preselection', 'preselection_tight', 'preselection_isolation', 'region_a', 'region_a_tagging_veto')
ORIGINS = ('prompt', 'pi0', 'eta', 'other')
DETAILS = ('prompt', 'pi0_separated', 'pi0_merged', 'pi0_single_contaminated', 'pi0_missing', 'pi0_other',
           'pi0_not_applicable', 'pi0_unknown', 'eta_separated', 'eta_merged', 'eta_single_contaminated',
           'eta_missing', 'eta_other', 'eta_gg_unknown', 'eta_parent', 'eta_daughters_unresolved', 'decay_photon',
           'charged', 'other_parent', 'other_parent_unresolved', 'residual_particle', 'mixed', 'unknown',
           'non_signal', 'photon_origin_unknown', 'eta_gg_not_applicable')


def root_file(path):
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie():
        raise RuntimeError(f'Cannot read {path}')
    return f


def flag(x):
    return bool(ord(x)) if isinstance(x, (str, bytes)) else bool(x)


def scalar(x):
    return ord(x) if isinstance(x, (str, bytes)) and len(x) == 1 else x


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    with path.open('w') as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def production_summary(new, old, output):
    # key = sample, selection, bin, origin, detail. Include flow bins explicitly.
    sums = defaultdict(lambda: [0., 0., 0.])
    bin_edges = {}
    normalizations = {}
    partials = sorted((new / 'partials/jet').glob('*/shard_*/photon_candidate_selection.root'))
    if not partials:
        raise RuntimeError('No partials found')
    for path in partials:
        sample = path.parent.parent.name
        f = root_file(path)
        m = f.Get('metadata'); m.GetEntry(0)
        if int(m.schema_version) != 11:
            raise RuntimeError(f'Unexpected partial schema: {path}')
        denom = float(m.sample_sum_generator_weights[0])
        if sample in normalizations and normalizations[sample] != denom:
            raise RuntimeError(f'Inconsistent full-sample normalization: {sample}')
        normalizations[sample] = denom
        for selection in SELECTIONS:
            h = f.Get(f'common_composition/{selection}/migration_count')
            w = f.Get(f'common_composition/{selection}/migration_weighted')
            if h.GetNbinsY() != 4 or h.GetNbinsZ() != len(DETAILS):
                raise RuntimeError('Unexpected migration axes')
            for b in range(h.GetNbinsX() + 2):
                bin_edges[b] = (h.GetXaxis().GetBinLowEdge(b), h.GetXaxis().GetBinUpEdge(b))
                for origin in range(4):
                    for detail in range(len(DETAILS)):
                        count = h.GetBinContent(b, origin + 1, detail + 1)
                        if not count:
                            continue
                        weight = w.GetBinContent(b, origin + 1, detail + 1)
                        sumw2 = w.GetBinError(w.GetBin(b, origin + 1, detail + 1)) ** 2
                        for scope in (sample, 'all'):
                            dest = sums[scope, selection, b, origin, detail]
                            for i, value in enumerate((count, weight, sumw2)):
                                dest[i] += value
        f.Close()
    # Check partial sums against the independently written merged migration CSV.
    merged = defaultdict(lambda: [0., 0.])
    with (new / 'results/jet/common_composition_migration.csv').open() as source:
        for row in csv.DictReader(source):
            key = ('all', row['selection'], int(row['bin']), ORIGINS.index(row['legacy_origin']), DETAILS.index(row['common_detail']))
            merged[key] = [float(row['count']), float(row['weighted_pb'])]
    for key in set(merged) | {k for k in sums if k[0] == 'all'}:
        for a, b in zip(sums.get(key, [0., 0.])[:2], merged.get(key, [0., 0.])):
            if not math.isclose(a, b, rel_tol=2e-5, abs_tol=1e-9):
                raise RuntimeError(f'Merged/partial disagreement: {key}: {a}, {b}')
    write_csv(output / 'migration_by_sample_et.csv', (
        dict(sample=s, selection=sel, bin=b, et_low=bin_edges[b][0], et_high=bin_edges[b][1], legacy_origin=ORIGINS[o],
             common_detail=DETAILS[d], count=int(v[0]), weighted_pb=v[1], sumw2=v[2])
        for (s, sel, b, o, d), v in sorted(sums.items())))
    totals = defaultdict(lambda: [0., 0.])
    for (s, sel, b, o, d), values in sums.items():
        for name in ('all', 'old_' + ORIGINS[o]):
            for i in range(2): totals[s, sel, name][i] += values[i]
        if o == 3 and DETAILS[d].startswith('pi0_'):
            for name in ('other_to_pi0', 'other_to_' + DETAILS[d]):
                for i in range(2): totals[s, sel, name][i] += values[i]
    rates = []
    for s in sorted({key[0] for key in sums}):
        for sel in SELECTIONS:
            denominator, weighted_denominator = totals[s, sel, 'old_other']
            numerator, weighted_numerator = totals[s, sel, 'other_to_pi0']
            rates.append(dict(sample=s, selection=sel, old_other_count=int(denominator), other_to_pi0_count=int(numerator),
                              fraction_of_old_other=numerator / denominator if denominator else None,
                              old_other_weighted_pb=weighted_denominator, other_to_pi0_weighted_pb=weighted_numerator,
                              weighted_fraction_of_old_other=weighted_numerator / weighted_denominator if weighted_denominator else None,
                              other_to_pi0_unknown_count=int(totals[s, sel, 'other_to_pi0_unknown'][0]),
                              other_to_pi0_not_applicable_count=int(totals[s, sel, 'other_to_pi0_not_applicable'][0])))
    write_csv(output / 'migration_rates.csv', rates)
    files = [root_file(p / 'results/jet/candidate_composition.root') for p in (old, new)]
    denominators = {}
    for sel in SELECTIONS:
        hs = [f.Get(f'{sel}/h_candidate_denominator_et_count') for f in files]
        if hs[0].GetNcells() != hs[1].GetNcells(): raise RuntimeError('Different ET binning')
        equal = all(hs[0].GetBinContent(i) == hs[1].GetBinContent(i) for i in range(hs[0].GetNcells()))
        denominators[sel] = dict(old=hs[0].GetEntries(), new=hs[1].GetEntries(), all_bins_identical=equal)
    for f in files: f.Close()
    return dict(partials=len(partials), full_sample_sum_generator_weights=normalizations, denominators=denominators,
                aggregate_rates=[r for r in rates if r['sample'] == 'all'])


def old_origin(v, c):
    prompt = flag(v['split_cluster_truth_prompt_cluster'][c]) and v['split_cluster_truth_dominant_fraction'][c] > .5
    pi0 = flag(v['split_cluster_pi0_anchor_valid'][c]) and v['split_cluster_pi0_anchor_main_fraction'][c] > .5
    eta = not pi0 and v['split_cluster_truth_eta_main_fraction'][c] > .5
    return (0 if prompt else 1 if pi0 else 2) if sum((prompt, pi0, eta)) == 1 else 3


def map_audit(args, output):
    scan_summary = []
    reason_rows, migration_rows = Counter(), Counter()
    with (output / 'representative_clusters.jsonl').open('w') as examples:
        for sample in args.samples:
            for index in args.map_indices:
                paths = [p / f'maps/{sample}/map_{index:06d}.root' for p in (args.old, args.new)]
                files = [root_file(p) for p in paths]
                trees = [f.Get('event_tree') for f in files]
                metas = [f.Get('metadata') for f in files]
                for m in metas: m.GetEntry(0)
                settings = ('model_sha256', 'min_cluster_energy', 'pi0_partner_min_energy', 'eta_partner_min_energy',
                            'pi0_mass_min', 'pi0_mass_max', 'eta_mass_min', 'eta_mass_max', 'missing_energy_min', 'missing_energy_max',
                            'n_events_processed', 'sum_generator_weight_processed', 'first_input_suffix', 'last_input_suffix')
                if any(getattr(metas[0], k) != getattr(metas[1], k) for k in settings): raise RuntimeError(f'Metadata mismatch: {sample}')
                if trees[0].GetEntries() != trees[1].GetEntries(): raise RuntimeError('Entry count mismatch')
                compare = ['split_cluster_id', 'split_cluster_e', 'split_cluster_et', 'split_cluster_pi0_tag', 'split_cluster_eta_tag']
                compare += [b.GetName() for b in trees[0].GetListOfBranches() if b.GetName().startswith(('split_cluster_pass_', 'split_cluster_bdt_', 'split_cluster_iso_'))]
                old_keys = set(compare) | {'split_cluster_truth_prompt_cluster', 'split_cluster_truth_dominant_fraction',
                    'split_cluster_pi0_anchor_valid', 'split_cluster_pi0_anchor_main_fraction', 'split_cluster_truth_eta_main_fraction',
                    'split_cluster_truth_dominant_g4_pdg_id', 'split_cluster_truth_dominant_g4_track_id',
                    'split_cluster_truth_dominant_immediate_parent_pdg', 'split_cluster_truth_dominant_classification_parent_pdg'}
                common = ('main', 'selected_pdg', 'selected_source', 'selected_id', 'selected_embedding', 'fraction', 'quality',
                          'applicability', 'topology_reason', 'legacy_origin', 'offset', 'ancestor_offset', 'table', 'pdg', 'source', 'id',
                          'decay', 'deposit', 'ancestor_pdg', 'ancestor_id', 'ancestor_source', 'piece_quality', 'unresolved_piece')
                new_keys = set(compare) | {'split_cluster_common_' + k for k in common}
                new_keys |= {f'split_cluster_{p}_anchor_{k}' for p in ('pi0', 'eta') for k in ('valid', 'topology')}
                scalars = ('event_uid', 'source_file_id', 'event_in_file', 'hepmc_event_number', 'sample_stitching_pass')
                for t, keys in zip(trees, (old_keys, new_keys)):
                    t.SetBranchStatus('*', 0)
                    for k in keys | set(scalars):
                        if not t.GetBranch(k): raise RuntimeError(f'Missing branch {k}')
                        t.SetBranchStatus(k, 1)
                stats = Counter(); saved = Counter()
                for entry in range(min(args.events_per_map, trees[0].GetEntries())):
                    for t in trees:
                        if t.GetEntry(entry) <= 0: raise RuntimeError('Cannot read entry')
                    if any(getattr(trees[0], k) != getattr(trees[1], k) for k in scalars): raise RuntimeError(f'Event identity mismatch: {entry}')
                    vs = [{k: list(getattr(t, k)) for k in keys} for t, keys in zip(trees, (old_keys, new_keys))]
                    a, n = vs
                    for k in compare:
                        if len(a[k]) != len(n[k]) or any(x != y and not (isinstance(x, float) and isinstance(y, float) and math.isnan(x) and math.isnan(y)) for x, y in zip(a[k], n[k])):
                            raise RuntimeError(f'Reconstruction changed: {sample}/{index}/{entry}/{k}')
                    stats['events_compared'] += 1
                    if not flag(trees[1].sample_stitching_pass): continue
                    stats['events_stitch_pass'] += 1
                    for c, cid in enumerate(n['split_cluster_id']):
                        if not flag(n['split_cluster_pass_kinematics'][c]): continue
                        stats['kinematic_clusters'] += 1
                        get = lambda k: n['split_cluster_common_' + k][c]
                        origin = old_origin(a, c)
                        if origin != get('legacy_origin'): stats['legacy_origin_mismatches'] += 1
                        main, reason, app = int(get('main')), int(get('topology_reason')), int(get('applicability'))
                        species = 'pi0' if main == 2 else 'eta'
                        topology = int(n[f'split_cluster_{species}_anchor_topology'][c]) if main in (2, 3) else -999
                        reason_rows[sample, main, app, reason, topology] += 1
                        if main in (2, 3) and app == 0 and topology in range(5): stats['unknown_with_saved_topology'] += 1
                        group = 'other_to_pi0' if origin == 3 and main == 2 else 'retained_pi0' if origin == 1 and main == 2 else 'eta_gg_unknown' if main == 3 and app == 0 else None
                        if not group: continue
                        start, end = n['split_cluster_common_offset'][c:c+2]
                        selected_rows = [j for j in range(start, end) if n['split_cluster_common_table'][j] == 0 and
                                         n['split_cluster_common_id'][j] == get('selected_id') and n['split_cluster_common_source'][j] == get('selected_source')]
                        decays = tuple(sorted({n['split_cluster_common_decay'][j] for j in selected_rows}))
                        ancestors = tuple(sorted({n['split_cluster_common_ancestor_pdg'][k] for j in selected_rows
                            for k in range(n['split_cluster_common_ancestor_offset'][j], n['split_cluster_common_ancestor_offset'][j+1])}))
                        migration_rows[sample, group, int(get('selected_source')), str(decays), str(ancestors),
                                       int(a['split_cluster_truth_dominant_g4_pdg_id'][c])] += 1
                        signature = (group, int(get('selected_source')), decays, ancestors, int(a['split_cluster_truth_dominant_g4_pdg_id'][c]))
                        if saved[signature] >= args.examples_per_group: continue
                        saved[signature] += 1
                        pieces = []
                        for j in range(start, end):
                            piece = {k: scalar(n['split_cluster_common_' + k][j]) for k in ('table', 'pdg', 'source', 'id', 'decay', 'deposit', 'piece_quality', 'unresolved_piece')}
                            lo, hi = n['split_cluster_common_ancestor_offset'][j:j+2]
                            piece['ancestors'] = [{k: scalar(n['split_cluster_common_ancestor_' + k][q]) for k in ('pdg', 'id', 'source')} for q in range(lo, hi)]
                            pieces.append(piece)
                        record = dict(sample=sample, map_index=index, old_map=str(paths[0]), new_map=str(paths[1]), entry=entry,
                                      source_file_id=int(trees[1].source_file_id), event_in_file=int(trees[1].event_in_file),
                                      event_uid=int(trees[1].event_uid), hepmc_event_number=int(trees[1].hepmc_event_number),
                                      cluster_id=int(cid), et=n['split_cluster_et'][c], group=group, old_origin=ORIGINS[origin],
                                      saved_topology=topology, pieces=pieces)
                        record.update({k: scalar(get(k)) for k in ('main', 'selected_source', 'selected_id', 'selected_embedding', 'fraction', 'quality', 'applicability', 'topology_reason')})
                        record['old_dominant'] = {k: scalar(a['split_cluster_truth_dominant_' + k][c]) for k in ('g4_pdg_id', 'g4_track_id', 'fraction', 'immediate_parent_pdg', 'classification_parent_pdg')}
                        record['old_pi0_main_fraction'] = a['split_cluster_pi0_anchor_main_fraction'][c]
                        examples.write(json.dumps(record) + '\n')
                scan_summary.append(dict(sample=sample, map_index=index, compared_vector_branches=len(compare), **stats))
                print(sample, index, dict(stats), flush=True)
                for f in files: f.Close()
    write_csv(output / 'sampled_topology_reasons.csv', (dict(sample=k[0], main=k[1], applicability=k[2], reason=k[3], saved_topology=k[4], count=v) for k,v in sorted(reason_rows.items())))
    write_csv(output / 'sampled_origin_details.csv', (dict(sample=k[0], group=k[1], selected_source=k[2], decay_codes=k[3], ancestor_pdgs=k[4], old_dominant_pdg=k[5], count=v) for k,v in sorted(migration_rows.items())))
    return scan_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('old', type=Path); parser.add_argument('new', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--samples', nargs='+', default=['jet3', 'jet12', 'jet40'])
    parser.add_argument('--map-indices', nargs='+', type=int, default=[0])
    parser.add_argument('--events-per-map', type=int, default=1000)
    parser.add_argument('--examples-per-group', type=int, default=2)
    args = parser.parse_args()
    if args.events_per_map <= 0 or args.examples_per_group <= 0 or any(i < 0 for i in args.map_indices): parser.error('Counts must be positive and indices nonnegative')
    args.output.mkdir(parents=True, exist_ok=False)
    ROOT.gROOT.SetBatch(True)
    summary = dict(old=str(args.old.resolve()), new=str(args.new.resolve()), sampling='First N written events in the specified maps; not a random sample.',
                   events_per_map=args.events_per_map, map_indices=args.map_indices)
    summary['production'] = production_summary(args.new, args.old, args.output)
    summary['map_scan'] = map_audit(args, args.output)
    summary['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
