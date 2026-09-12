"""Synthetic integration checks. Run with ana.565 Python/PyROOT; writes only a temporary directory."""
import array
import math
from pathlib import Path
import tempfile
import ROOT

ROOT.gROOT.SetBatch(True)
workflow = Path(__file__).resolve().parent
ROOT.gROOT.ProcessLine(f'.L {workflow}/ReduceTaggingMassDiagnostic.C')
ROOT.gROOT.ProcessLine(f'.L {workflow}/MergeTaggingMassDiagnostic.C')


def make_map(path, pi0_cut=0.15, wrong_veto=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    f = ROOT.TFile(str(path), 'RECREATE')
    metadata = ROOT.TTree('metadata', '')
    keep = []
    def scalar(tree, name, value, code='d'):
        data = array.array(code, [value]); keep.append(data)
        tree.Branch(name, data, name + '/' + {'d': 'D', 'i': 'I', 'I': 'i', 'q': 'L', 'B': 'b'}[code])
        return data
    for name, value in [('schema_version', 5), ('pi0_topology_algorithm_version', 11)]: scalar(metadata, name, value, 'i')
    scalar(metadata, 'map_chunk_id', 0, 'I')
    scalar(metadata, 'manifest_begin', 0, 'q'); scalar(metadata, 'manifest_end', 10000, 'q')
    for name, value in {'sum_generator_weight_processed': 10., 'min_cluster_energy': .12, 'partner_diagnostic_min_cluster_energy': 0.,
                        'pi0_partner_min_energy': pi0_cut, 'eta_partner_min_energy': .5, 'pi0_mass_min': .1, 'pi0_mass_max': .2,
                        'eta_mass_min': .45, 'eta_mass_max': .65, 'missing_energy_min': .2, 'missing_energy_max': .5,
                        'min_photon_energy_recovery': 0., 'sample_cross_section_pb': 10., 'sample_window_min': 8.,
                        'sample_window_max': 12., 'sample_upper_unbounded': 0., 'shower_shape_min_tower_energy': 0.,
                        'candidate_et_min': 5., 'candidate_et_max': 100., 'candidate_abs_eta_max': .7, 'max_abs_vertex_z': 60.,
                        'isolation_radius': .3, 'isolation_scale': 1., 'isolation_offset': 0., 'nonisolation_gap': 1.}.items(): scalar(metadata, name, value)
    for name, value in [('sample_name', 'jet8'), ('analysis_release', 'synthetic'), ('model_sha256', 'synthetic')]:
        data = ROOT.std.string(value); keep.append(data); metadata.Branch(name, data)
    metadata.Fill(); metadata.Write()
    tree = ROOT.TTree('event_tree', '')
    for name in ['event_weight_valid', 'sample_stitching_valid', 'sample_stitching_pass']: scalar(tree, name, 1, 'B')
    numerator = scalar(tree, 'weight_numerator_pb', 20.)
    n = scalar(tree, 'split_ncluster', 0, 'I')
    vectors = {}
    types = {'unsigned int': ['id'], 'double': ['e', 'et', 'eta', 'phi'], 'unsigned char': ['pass_region_a', 'pi0_anchor_valid', 'truth_prompt_cluster', 'pi0_tag', 'eta_tag'],
             'float': ['pi0_anchor_main_fraction', 'pi0_anchor_truth_partner_mass', 'pi0_anchor_truth_partner_cluster_e', 'truth_dominant_fraction'],
             'int': ['pi0_anchor_topology', 'pi0_anchor_truth_partner_cluster_id', 'pi0_anchor_truth_partner_tag_status']}
    for typ, names in types.items():
        for name in names:
            v = ROOT.std.vector(typ)(); vectors[name] = v; tree.Branch('split_cluster_' + name, v)
    # Three eligible pi0 window partners, one eta window partner, plus an exact-threshold excluded pi0 partner.
    energies = [10., 1., 1., 1., .15, .5]
    target_masses = [0., .12, .16, .55, .13, .14]
    for event in range(5):
        prompt = event < 2
        numerator[0] = 20. if event != 1 else -5.
        n[0] = len(energies) if prompt else 1
        for v in vectors.values(): v.clear()
        for j in range(n[0]):
            values = {'id': j + 1, 'e': energies[j], 'et': energies[j], 'eta': 0.,
                      'phi': 0. if j == 0 else math.acos(1 - target_masses[j] ** 2 / (20 * energies[j])),
                      'pass_region_a': j == 0, 'pi0_anchor_valid': not prompt, 'pi0_anchor_main_fraction': 1. if not prompt else 0.,
                      'pi0_anchor_topology': 3 if event == 2 else 1 if not prompt else 0,
                      'pi0_anchor_truth_partner_tag_status': 4 if event == 2 else 1, 'pi0_anchor_truth_partner_mass': .08 if event == 2 else -1. if event == 4 else .135,
                      'pi0_anchor_truth_partner_cluster_e': .1 if event == 2 else 1., 'pi0_anchor_truth_partner_cluster_id': 99,
                      'truth_prompt_cluster': prompt and j == 0, 'truth_dominant_fraction': 1.,
                      'pi0_tag': prompt and j == 0 and not wrong_veto, 'eta_tag': prompt and j == 0}
            for key, v in vectors.items(): v.push_back(values[key])
        tree.Fill()
    tree.Write(); f.Close()


def hist(path, name):
    f = ROOT.TFile.Open(str(path)); h = f.Get(name).Clone(); h.SetDirectory(0); f.Close(); return h


with tempfile.TemporaryDirectory(prefix='tagging-mass-test-') as temporary:
    base = Path(temporary)
    make_map(base / 'maps/jet8/map_000000.root')
    args = ('jet', str(base / 'maps'))
    assert ROOT.ReduceTaggingMassDiagnostic(*args, str(base / 'serial/jet8/shard_0'), True, 'jet8', 0, 1, 0) == 0
    for shard in range(2):
        assert ROOT.ReduceTaggingMassDiagnostic(*args, str(base / f'split/jet8/shard_{shard}'), True, 'jet8', shard, 2, 0) == 0
    assert ROOT.MergeTaggingMassDiagnostic('jet', str(base / 'split'), str(base / 'merged'), 'jet8', 2, True, True) == 0
    serial = base / 'serial/jet8/shard_0/tagging_mass_diagnostic.root'
    merged = base / 'merged/tagging_mass_diagnostic.root'
    groups = ['separated_truth_pair', 'separated_truth_pair_below_threshold', 'prompt_pi0_all_pairs', 'prompt_pi0_window_pairs', 'prompt_eta_all_pairs', 'prompt_eta_window_pairs']
    expected = [1, 0, 8, 6, 6, 2]
    for group, count in zip(groups, expected):
        name = f'{group}/h_{group}_mass_count'
        h = hist(serial, name); assert h.Integral(0, h.GetNbinsX() + 1) == count, (group, h.Integral(), count)
        for suffix in ['mass_count', 'mass_pb', 'mass_vs_anchor_et_count', 'mass_vs_anchor_et_pb', 'window_count', 'window_pb', 'window_pb_fraction']:
            a = hist(serial, f'{group}/h_{group}_{suffix}'); b = hist(merged, f'{group}/h_{group}_{suffix}')
            for bin in range(a.GetNcells()):
                assert math.isclose(a.GetBinContent(bin), b.GetBinContent(bin), abs_tol=1e-12)
                assert math.isclose(a.GetBinError(bin), b.GetBinError(bin), abs_tol=1e-12)
    h = hist(serial, 'prompt_pi0_window_pairs/h_prompt_pi0_window_pairs_mass_pb')
    assert math.isclose(h.Integral(), 4.5)  # 3 partners * (2 - .5)
    flow = hist(serial, 'candidate_flow_count')
    assert flow.Integral(0, flow.GetNbinsX() + 1, 5, 5) == 2 and flow.Integral(0, flow.GetNbinsX() + 1, 8, 8) == 2
    assert flow.Integral(0, flow.GetNbinsX() + 1, 1, 1) == 2 and flow.Integral(0, flow.GetNbinsX() + 1, 4, 4) == 1
    assert len(list((base / 'merged').rglob('*.pdf'))) == 24
    et_axis = flow.GetXaxis()
    assert [et_axis.GetBinLowEdge(i) for i in range(1, et_axis.GetNbins() + 2)] == [0, 5, 6, 8, 10, 18, 30, 40]
    # Strict mass boundaries, no nominal-mass best-pair selection.
    ROOT.gInterpreter.Declare('bool check_tagging_mass_boundaries() { tagging_mass::Histograms h("boundary"); h.fill(.1, 10, 1, .1, .2); h.fill(.2, 10, 1, .1, .2); h.fill(.15, 10, 1, .1, .2); return h.window_count->Integral(0, h.window_count->GetNbinsX() + 1, 2, 2) == 1; }')
    assert ROOT.check_tagging_mass_boundaries()
    # Refuse overwrite, incomplete coverage and mixed map settings.
    assert ROOT.ReduceTaggingMassDiagnostic(*args, str(base / 'serial/jet8/shard_0'), True, 'jet8', 0, 1, 0) != 0
    assert ROOT.MergeTaggingMassDiagnostic('jet', str(base / 'split'), str(base / 'bad_coverage'), 'jet8', 1, True, False) != 0
    make_map(base / 'bad_veto_maps/jet8/map_000000.root', wrong_veto=True)
    assert ROOT.ReduceTaggingMassDiagnostic('jet', str(base / 'bad_veto_maps'), str(base / 'bad_veto'), True, 'jet8') == 6
    assert not (base / 'bad_veto/tagging_mass_diagnostic.root').exists()
    make_map(base / 'bad_maps/jet8/map_000000.root', .1)
    assert ROOT.ReduceTaggingMassDiagnostic('jet', str(base / 'bad_maps'), str(base / 'bad_threshold'), True, 'jet8') != 0
    bad = ROOT.TFile.Open(str(base / 'split/jet8/shard_1/tagging_mass_diagnostic.root'), 'UPDATE')
    meta = ROOT.tagging_mass.Metadata(); assert meta.read(bad)
    meta.selection_settings[1] = .2
    bad.Delete('metadata;*'); meta.write(); bad.Close()
    assert ROOT.MergeTaggingMassDiagnostic('jet', str(base / 'split'), str(base / 'bad_settings'), 'jet8', 2, True, False) != 0
print('Synthetic reduce/merge, pair counts, negative weights, boundaries, veto flags, plots and rejection checks passed')
