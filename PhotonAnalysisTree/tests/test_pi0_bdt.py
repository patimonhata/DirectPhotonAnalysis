import ROOT
import os
import tempfile
from pathlib import Path

repo = Path(__file__).resolve().parents[2]
os.chdir(repo)
ROOT.gROOT.SetBatch(True)
assert ROOT.gROOT.LoadMacro('PhotonAnalysisTree/workflows/photon_candidate_selection/MergePythiaPhotonCandidateSelection.C') == 0
with tempfile.TemporaryDirectory(prefix="pi0_bdt_test_") as temporary:
    os.chdir(temporary)
    ROOT.gInterpreter.Declare(r'''
int test_pi0_bdt_roundoff() {
  using namespace candidate_composition;
  Histograms composition(40, 40); Spectra topology(40, 40); pi0_bdt::Histograms bdt;
  // Jet20's merged count and per-cluster weight reproduced the production hold.
  constexpr int count = 442703;
  constexpr double weight = 62623.0 / 10000000.0;
  for (int i = 0; i < count; ++i) {
    const double et = 5.5 + i % 29;
    composition.fill(pi0_merged, et, weight);
    topology.fill(3, et, weight);
    if (!bdt.fill(1, .9, true, false, false, et, weight)) return 1;
  }
  const double sequential = bdt.cells[0][1]->summary.GetBinContent(3);
  const double binned = composition.weighted[pi0_merged]->Integral(0, 41);
  if (same_double(sequential, binned)) return 2; // The old check must reject this regression fixture.
  if (!consistent_pi0_bdt(composition, topology, bdt)) return 3;
  // A real weight mismatch must still fail, even with a consistent isolation partition.
  for (int group : {0, 3}) bdt.cells[group][1]->summary.SetBinContent(3, sequential + 0.01);
  if (consistent_pi0_bdt(composition, topology, bdt)) return 4;
  for (int group : {0, 3}) {
    bdt.cells[group][1]->summary.SetBinContent(3, sequential);
    bdt.cells[group][1]->summary.SetBinContent(1, count + 1);
  }
  if (consistent_pi0_bdt(composition, topology, bdt)) return 5;
  if (same_weighted_integral(NAN, binned) || same_weighted_integral(INFINITY, binned)) return 6;
  return 0;
}

int test_pi0_bdt() {
  using namespace candidate_composition;
  if (origin_category(false, .5 > kMajorityThreshold, false, 1) != other) return 20;
  if (origin_category(false, .6 > kMajorityThreshold, false, 1) != pi0_separated) return 21;
  if (origin_category(true, true, false, 1) != other) return 22;
  if (origin_category(false, true, false, 0) != pi0_other) return 23;
  pi0_bdt::Histograms h;
  if (!h.fill(0, .3, true, true, false, 10, 2)) return 1;
  if (!h.fill(0, 1, true, false, true, 20, 3)) return 2;
  if (!h.fill(0, NAN, false, false, false, 15, 4)) return 3;
  if (h.fill(0, .5, true, true, true, 10, 1)) return 4;
  if (!h.valid() || h.cells[0][0]->count.Integral() != 2 || h.cells[0][0]->weighted.Integral() != 5 || h.cells[0][0]->summary.GetBinContent(2) != 1) return 5;
  { TFile f("pi0_bdt_test.root", "RECREATE"); if (!h.write(f)) return 6; }
  pi0_bdt::Histograms merged;
  { TFile f("pi0_bdt_test.root"); if (!merged.add(f) || !merged.add(f)) return 7; }
  if (!merged.valid() || merged.cells[0][0]->weighted.Integral() != 10) return 8;
  candidate_composition::Histograms composition(40,40); Spectra topology(40,40);
  composition.fill(candidate_composition::pi0_separated,10,2); composition.fill(candidate_composition::pi0_separated,20,3); composition.fill(candidate_composition::pi0_separated,15,4);
  topology.fill(2,10,2); topology.fill(2,20,3); topology.fill(2,15,4);
  if (!consistent_pi0_bdt(composition,topology,h)) return 9;
  topology.fill(2,10,1); if (consistent_pi0_bdt(composition,topology,h)) return 10;
  gSystem->mkdir("pi0_bdt_qa",true);
  if (!h.export_files("pi0_bdt_qa","Kinematics",false)) return 11;
  pi0_bdt::Histograms isolated; isolated.fill(1,.9,true,true,false,20,1);
  gSystem->mkdir("pi0_bdt_qa_isolated",true);
  if (!isolated.export_files("pi0_bdt_qa_isolated","Region A",true)) return 12;
  for (int i = 0; i < 6; ++i) {
    const std::string directory = std::string("selections/") + kSelectionKeys[i];
    gSystem->mkdir(directory.c_str(), true);
    if (!isolated.export_files(directory, kSelectionLabels[i], i == 3 || i >= 4)) return 30 + i;
    for (const char* name : {"pi0_bdt_score_by_topology.pdf", "pi0_bdt_score_by_topology_isolation.pdf", "pi0_bdt_score.root", "pi0_bdt_score_summary.csv"})
      if (gSystem->AccessPathName((directory + "/" + name).c_str())) return 40 + i;
  }
  return 0;
}
''')
    roundoff_result = ROOT.test_pi0_bdt_roundoff()
    print("roundoff regression result:", roundoff_result)
    assert roundoff_result == 0
    result = ROOT.test_pi0_bdt()
    print('test result:', result)
    assert result == 0
