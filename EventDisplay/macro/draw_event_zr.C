#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
const double kZrEmcalInnerRadius = 93.0;
const double kZrEmcalOuterRadius = 113.0;
double marker_size(double energy)
{
  if (!std::isfinite(energy) || energy < 0.0)
  {
    return 0.8;
  }
  return std::min(2.4, 0.8 + 0.18 * energy);
}

bool has_track(const std::vector<int>& track_ids, int track_id)
{
  for (int stored_track_id : track_ids)
  {
    if (stored_track_id == track_id)
    {
      return true;
    }
  }
  return false;
}

int truth_segment_color(bool draw_pi0, bool draw_gamma, bool draw_electron)
{
  if (draw_pi0)
  {
    return kRed + 1;
  }
  if (draw_gamma)
  {
    return kBlue + 1;
  }
  if (draw_electron)
  {
    return kMagenta + 1;
  }
  return kGray + 1;
}
}

void draw_event_zr(const char* filename = "event_display.root", int event_id = 0, bool draw_other_truth = false, bool draw_hits = false)
{
  const std::string output_filename = Form("/sphenix/u/ryotaro/DirectPhotonAnalysis/EventDisplay/output/image/event_%d_zr.pdf",event_id);
  gStyle->SetOptStat(0);

  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie())
  {
    std::cout << "draw_event_zr - cannot open " << filename << std::endl;
    return;
  }

  TTree* segments = static_cast<TTree*>(file->Get("truth_segments"));
  TTree* clusters = static_cast<TTree*>(file->Get("cemc_clusters"));
  TTree* hits = static_cast<TTree*>(file->Get("cemc_hits"));
  if (!segments || !clusters)
  {
    std::cout << "draw_event_zr - missing truth_segments or cemc_clusters tree" << std::endl;
    file->Close();
    return;
  }

  gROOT->cd();
  TCanvas* canvas = new TCanvas("canvas_zr", "event display z-r", 1100, 800);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(0.04);
  canvas->SetTopMargin(0.07);
  canvas->SetBottomMargin(0.12);
  TH2F* frame = new TH2F("frame_zr", Form("Event %d;z [cm];r [cm]", event_id), 100, -180.0, 180.0, 100, 0.0, 140.0);
  frame->SetDirectory(nullptr);
  frame->Draw();

  TBox* emcal_band = new TBox(-180.0, kZrEmcalInnerRadius, 180.0, kZrEmcalOuterRadius);
  emcal_band->SetFillColorAlpha(kGray + 1, 0.18);
  emcal_band->SetLineColor(kGray + 2);
  emcal_band->SetLineStyle(2);
  emcal_band->SetLineWidth(2);
  emcal_band->Draw();

  int event = 0;
  int track_id = 0;
  int pid = 0;
  int parent_id = 0;
  int segment_type = 0;
  double z0 = 0.0;
  double z1 = 0.0;
  double r0 = 0.0;
  double r1 = 0.0;
  segments->SetBranchAddress("event", &event);
  segments->SetBranchAddress("track_id", &track_id);
  segments->SetBranchAddress("pid", &pid);
  segments->SetBranchAddress("parent_id", &parent_id);
  segments->SetBranchAddress("segment_type", &segment_type);
  segments->SetBranchAddress("z0", &z0);
  segments->SetBranchAddress("z1", &z1);
  segments->SetBranchAddress("r0", &r0);
  segments->SetBranchAddress("r1", &r1);

  std::vector<int> pi0_track_ids;
  for (Long64_t i = 0; i < segments->GetEntries(); ++i)
  {
    segments->GetEntry(i);
    if (event == event_id && segment_type == 1 && pid == 111)
    {
      pi0_track_ids.push_back(track_id);
    }
  }

  for (Long64_t i = 0; i < segments->GetEntries(); ++i)
  {
    segments->GetEntry(i);
    if (event != event_id)
    {
      continue;
    }

    const bool draw_pi0 = segment_type == 1;
    const bool draw_gamma = segment_type == 2;
    const bool draw_electron = has_track(pi0_track_ids, parent_id) && std::abs(pid) == 11;
    if (!draw_other_truth && !draw_pi0 && !draw_gamma && !draw_electron)
    {
      continue;
    }

    TLine* line = new TLine(z0, r0, z1, r1);
    line->SetLineWidth(draw_pi0 ? 3 : 2);
    line->SetLineColor(truth_segment_color(draw_pi0, draw_gamma, draw_electron));
    line->Draw();
  }

  double energy = 0.0;
  double z = 0.0;
  double r = 0.0;
  clusters->SetBranchAddress("event", &event);
  clusters->SetBranchAddress("energy", &energy);
  clusters->SetBranchAddress("z", &z);
  clusters->SetBranchAddress("r", &r);
  TLatex latex;
  latex.SetTextSize(0.026);
  latex.SetTextAlign(12);
  for (Long64_t i = 0; i < clusters->GetEntries(); ++i)
  {
    clusters->GetEntry(i);
    if (event != event_id)
    {
      continue;
    }
    TMarker* marker = new TMarker(z, r, 20);
    marker->SetMarkerColor(kGreen + 2);
    marker->SetMarkerSize(marker_size(energy));
    marker->Draw();
    latex.DrawLatex(z + 2.0, r + 2.0, Form("%.2f GeV", energy));
  }

  if (draw_hits && hits)
  {
    double hz0 = 0.0;
    double hz1 = 0.0;
    double hr0 = 0.0;
    double hr1 = 0.0;
    hits->SetBranchAddress("event", &event);
    hits->SetBranchAddress("z0", &hz0);
    hits->SetBranchAddress("z1", &hz1);
    hits->SetBranchAddress("r0", &hr0);
    hits->SetBranchAddress("r1", &hr1);
    for (Long64_t i = 0; i < hits->GetEntries(); ++i)
    {
      hits->GetEntry(i);
      if (event != event_id)
      {
        continue;
      }
      TMarker* marker = new TMarker(0.5 * (hz0 + hz1), 0.5 * (hr0 + hr1), 6);
      marker->SetMarkerColor(kGray + 1);
      marker->SetMarkerSize(0.35);
      marker->Draw();
    }
  }

  TLegend* legend = new TLegend(0.14, 0.70, 0.42, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(emcal_band, "EMCal 93-113 cm", "f");
  legend->AddEntry((TObject*) 0, "red: #pi^{0} truth", "");
  legend->AddEntry((TObject*) 0, "blue: #gamma truth", "");
  legend->AddEntry((TObject*) 0, "magenta: e^{#pm} truth", "");
  legend->Draw();

  canvas->SaveAs(output_filename.c_str());
  // canvas->SaveAs(Form("event_%d_zr.pdf", event_id));
  // canvas->SaveAs(Form("event_%d_zr.png", event_id));
  segments->ResetBranchAddresses();
  clusters->ResetBranchAddresses();
  if (draw_hits && hits)
  {
    hits->ResetBranchAddresses();
  }
  file->Close();
}
