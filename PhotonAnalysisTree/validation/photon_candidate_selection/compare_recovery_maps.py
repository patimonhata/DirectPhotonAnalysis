"""Check that changing only recovery preserves reconstructed selections on identical inputs."""
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
a = ROOT.TFile.Open(sys.argv[1])
b = ROOT.TFile.Open(sys.argv[2])
left, right = a.Get('event_tree'), b.Get('event_tree')
assert left.GetEntries() == right.GetEntries()
fixed = [branch.GetName() for branch in left.GetListOfBranches()
         if branch.GetName().startswith('split_cluster_pass_')]
fixed += ['split_cluster_id', 'split_cluster_e', 'split_cluster_et', 'split_cluster_eta', 'split_cluster_phi',
          'split_cluster_pi0_tag', 'split_cluster_eta_tag', 'split_cluster_pi0_partner_cluster_id',
          'split_cluster_eta_partner_cluster_id', 'split_cluster_pi0_partner_mass', 'split_cluster_eta_partner_mass']
changed = total = 0
for entry in range(left.GetEntries()):
    left.GetEntry(entry)
    right.GetEntry(entry)
    assert left.event_uid == right.event_uid
    for name in fixed:
        assert list(getattr(left, name)) == list(getattr(right, name)), (entry, name)
    old = list(left.split_cluster_pi0_anchor_topology)
    new = list(right.split_cluster_pi0_anchor_topology)
    assert len(old) == len(new)
    total += sum(x != -999 for x in old)
    changed += sum(x != y for x, y in zip(old, new))
assert changed > 0, 'Use a QA sample with a measurable recovery effect'
print(f'recovery comparison: {left.GetEntries()} events, all {len(fixed)} selection/kinematic branches unchanged; {changed}/{total} anchor topologies changed')
