"""Independently check a small scan against input maps; not intended for full production."""
import math
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
pattern, result = sys.argv[1:]
f = ROOT.TFile.Open(result)
meta = f.Get('sources/partial_0/metadata'); meta.GetEntry(0)
sumw = float(meta.sumw)
cols = ['split_cluster_id', 'split_cluster_e', 'split_cluster_et', 'split_cluster_eta', 'split_cluster_phi',
        'split_cluster_pass_region_a', 'split_cluster_truth_prompt_cluster', 'split_cluster_truth_dominant_fraction', 'weight_numerator_pb', 'split_cluster_pi0_anchor_valid', 'split_cluster_pi0_anchor_main_fraction',
        'split_cluster_truth_contributor_offset', 'split_cluster_truth_contributor_g4_pdg_id', 'split_cluster_truth_contributor_embedding_id',
        'split_cluster_truth_contributor_fraction', 'split_cluster_truth_contributor_photon_source', 'split_cluster_truth_contributor_classification_parent_pdg']
data = ROOT.RDataFrame('event_tree', pattern).Filter('event_weight_valid && sample_stitching_valid && sample_stitching_pass && region_a_count > 0').AsNumpy(cols)
expected = {}
for origin in ['prompt', 'background', 'pi0']:
    for kind in ['count', 'pb']:
        name = origin + '_survival_mask_' + kind
        h = ROOT.TH2D('expected_' + name, '', 40, 0, 40, 64, -.5, 63.5); h.SetDirectory(0); h.Sumw2(); expected[name] = h

def flag(x):
    return ord(x) != 0 if isinstance(x, str) else bool(x)

candidates = 0
for row in range(len(data[cols[0]])):
    values = {name: data['split_cluster_' + name][row] for name in ['id', 'e', 'et', 'eta', 'phi', 'pass_region_a', 'truth_prompt_cluster', 'truth_dominant_fraction', 'pi0_anchor_valid', 'pi0_anchor_main_fraction', 'truth_contributor_offset', 'truth_contributor_g4_pdg_id', 'truth_contributor_embedding_id', 'truth_contributor_fraction', 'truth_contributor_photon_source', 'truth_contributor_classification_parent_pdg']}
    for i in range(len(values['id'])):
        if not flag(values['pass_region_a'][i]):
            continue
        candidates += 1
        pairs = []
        for j in range(len(values['id'])):
            if values['id'][i] == values['id'][j]:
                continue
            mass2 = 2 * values['e'][i] * values['e'][j] / (math.cosh(values['eta'][i]) * math.cosh(values['eta'][j])) * (math.cosh(values['eta'][i] - values['eta'][j]) - math.cos(values['phi'][i] - values['phi'][j]))
            pairs.append((values['e'][j], math.sqrt(max(0, mass2))))
        eta_tag = any(e > .5 and .45 < mass < .65 for e, mass in pairs)
        mask = sum(1 << w for w, (cut, low, high) in enumerate([(.15, .08, .18), (.15, .08, .20), (.12, .08, .18), (.12, .08, .20), (.5, .10, .20), (.15, .10, .20)])
                   if not eta_tag and not any(e > cut and low < mass < high for e, mass in pairs))
        is_prompt = flag(values['truth_prompt_cluster'][i]) and values['truth_dominant_fraction'][i] > .5
        is_pi0 = flag(values['pi0_anchor_valid'][i]) and values['pi0_anchor_main_fraction'][i] > .5
        ef = 0.
        for j in range(values['truth_contributor_offset'][i], values['truth_contributor_offset'][i+1]):
            pdg = values['truth_contributor_g4_pdg_id'][j]
            if values['truth_contributor_embedding_id'][j] == 1 and (abs(pdg) == 221 or (pdg == 22 and (values['truth_contributor_photon_source'][j] == 3 or values['truth_contributor_classification_parent_pdg'][j] == 221))):
                ef += values['truth_contributor_fraction'][j]
        exclusive = int(is_prompt) + int(is_pi0) + int(ef > .5) <= 1
        origins = ['prompt' if is_prompt and exclusive else 'background']
        if is_pi0 and exclusive:
            origins.append('pi0')
        for origin in origins:
            expected[origin + '_survival_mask_count'].Fill(values['et'][i], mask)
            expected[origin + '_survival_mask_pb'].Fill(values['et'][i], mask, data['weight_numerator_pb'][row] / sumw)

for name, h in expected.items():
    actual = f.Get(name)
    for b in range(h.GetNcells()):
        assert math.isclose(h.GetBinContent(b), actual.GetBinContent(b), rel_tol=1e-10, abs_tol=1e-9), (name, b)
        assert math.isclose(h.GetBinError(b), actual.GetBinError(b), rel_tol=1e-10, abs_tol=1e-9), (name, b, 'error')
print(f'Independent six-working-point and pi0 comparison passed: {candidates} candidates')
