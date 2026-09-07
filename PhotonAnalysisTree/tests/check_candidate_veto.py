"""Recompute all-pair vetoes from a QA map whose partner cuts exceed its stored-cluster cut."""
import math
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
source = ROOT.TFile.Open(sys.argv[1])
metadata = source.Get('metadata')
metadata.GetEntry(0)
settings = [(meson, float(getattr(metadata, meson + '_partner_min_energy')),
             float(getattr(metadata, meson + '_mass_min')), float(getattr(metadata, meson + '_mass_max')),
             nominal) for meson, nominal in [('pi0', 0.134977), ('eta', 0.547862)]]
assert all(cut >= metadata.min_cluster_energy for _, cut, _, _, _ in settings), 'QA input must store every eligible partner'
def flag(value):
    return ord(value) != 0 if isinstance(value, str) else bool(value)

checked = tags = 0
for event in source.Get('event_tree'):
    clusters = list(zip(event.split_cluster_id, event.split_cluster_e, event.split_cluster_eta, event.split_cluster_phi))
    for i, (cid, energy, eta, phi) in enumerate(clusters):
        for meson, cut, low, high, nominal in settings:
            pairs = []
            for pid, pe, peta, pphi in clusters:
                if pid == cid or pe <= cut:
                    continue
                mass2 = 2 * energy * pe / (math.cosh(eta) * math.cosh(peta)) * (math.cosh(eta - peta) - math.cos(phi - pphi))
                mass = math.sqrt(max(0, mass2))
                if low < mass < high:
                    pairs.append((abs(mass - nominal), pid, mass))
            assert flag(getattr(event, 'split_cluster_' + meson + '_tag')[i]) == bool(pairs)
            if pairs:
                best = min(pairs)
                assert getattr(event, 'split_cluster_' + meson + '_partner_cluster_id')[i] == best[1]
                assert math.isclose(getattr(event, 'split_cluster_' + meson + '_partner_mass')[i], best[2], abs_tol=1e-6)
                tags += 1
            checked += 1
        expected = flag(event.split_cluster_pass_region_a[i]) and not flag(event.split_cluster_pi0_tag[i]) and not flag(event.split_cluster_eta_tag[i])
        assert flag(event.split_cluster_pass_final_photon[i]) == expected
print(f'all-pair veto: {checked} decisions and {tags} tags passed')
