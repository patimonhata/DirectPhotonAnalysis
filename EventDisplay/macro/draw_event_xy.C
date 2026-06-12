#include <TCanvas.h>
#include <TEllipse.h>
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

double xy_marker_size(double energy);
const double kXyEmcalInnerRadius = 93.0;
const double kXyEmcalOuterRadius = 113.0;

bool xy_has_track(const std::vector<int>& track_ids, int track_id)
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

int xy_truth_segment_color(bool draw_pi0, bool draw_gamma, bool draw_electron)
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

void draw_event_xy(const char* filename = "event_display.root", int event_id = 0, bool draw_other_truth = false, bool draw_hits = false)
{
  const std::string output_filename = Form("/sphenix/u/ryotaro/DirectPhotonAnalysis/EventDisplay/output/image/event_%d_xy.pdf",event_id);
  gStyle->SetOptStat(0);

  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie())
  {
    std::cout << "draw_event_xy - cannot open " << filename << std::endl;
    return;
  }

  TTree* segments = static_cast<TTree*>(file->Get("truth_segments"));
  TTree* clusters = static_cast<TTree*>(file->Get("cemc_clusters"));
  TTree* hits = static_cast<TTree*>(file->Get("cemc_hits"));
  if (!segments || !clusters)
  {
    std::cout << "draw_event_xy - missing truth_segments or cemc_clusters tree" << std::endl;
    file->Close();
    return;
  }

  const double frame_r = std::max(130.0, kXyEmcalOuterRadius + 30.0);
  gROOT->cd();
  TCanvas* canvas = new TCanvas("canvas_xy", "event display x-y", 900, 900);
  canvas->SetLeftMargin(0.14);
  canvas->SetRightMargin(0.04);
  canvas->SetTopMargin(0.06);
  canvas->SetBottomMargin(0.12);
  TH2F* frame = new TH2F("frame_xy", Form("Event %d;x [cm];y [cm]", event_id), 100, -frame_r, frame_r, 100, -frame_r, frame_r);
  frame->SetDirectory(nullptr);
  frame->Draw();

  TEllipse* emcal_band = new TEllipse(0.0, 0.0, kXyEmcalOuterRadius, kXyEmcalOuterRadius);
  emcal_band->SetFillColorAlpha(kGray + 1, 0.18);
  emcal_band->SetLineColor(kGray + 2);
  emcal_band->SetLineStyle(2);
  emcal_band->SetLineWidth(2);
  emcal_band->Draw();

  TEllipse* emcal_inner = new TEllipse(0.0, 0.0, kXyEmcalInnerRadius, kXyEmcalInnerRadius);
  emcal_inner->SetFillColor(kWhite);
  emcal_inner->SetLineColor(kGray + 2);
  emcal_inner->SetLineStyle(2);
  emcal_inner->SetLineWidth(2);
  emcal_inner->Draw();

  int event = 0;
  int track_id = 0;
  int pid = 0;
  int parent_id = 0;
  int segment_type = 0;
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
  segments->SetBranchAddress("event", &event);
  segments->SetBranchAddress("track_id", &track_id);
  segments->SetBranchAddress("pid", &pid);
  segments->SetBranchAddress("parent_id", &parent_id);
  segments->SetBranchAddress("segment_type", &segment_type);
  segments->SetBranchAddress("x0", &x0);
  segments->SetBranchAddress("y0", &y0);
  segments->SetBranchAddress("x1", &x1);
  segments->SetBranchAddress("y1", &y1);

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
    const bool draw_electron = xy_has_track(pi0_track_ids, parent_id) && std::abs(pid) == 11;
    if (!draw_other_truth && !draw_pi0 && !draw_gamma && !draw_electron)
    {
      continue;
    }

    TLine* line = new TLine(x0, y0, x1, y1);
    line->SetLineWidth(draw_pi0 ? 3 : 2);
    line->SetLineColor(xy_truth_segment_color(draw_pi0, draw_gamma, draw_electron));
    line->Draw();
  }

  double energy = 0.0;
  double x = 0.0;
  double y = 0.0;
  clusters->SetBranchAddress("event", &event);
  clusters->SetBranchAddress("energy", &energy);
  clusters->SetBranchAddress("x", &x);
  clusters->SetBranchAddress("y", &y);
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
    TMarker* marker = new TMarker(x, y, 20);
    marker->SetMarkerColor(kGreen + 2);
    marker->SetMarkerSize(xy_marker_size(energy));
    marker->Draw();
    latex.DrawLatex(x + 2.0, y + 2.0, Form("%.2f GeV", energy));
  }

  if (draw_hits && hits)
  {
    double hx0 = 0.0;
    double hy0 = 0.0;
    double hx1 = 0.0;
    double hy1 = 0.0;
    hits->SetBranchAddress("event", &event);
    hits->SetBranchAddress("x0", &hx0);
    hits->SetBranchAddress("y0", &hy0);
    hits->SetBranchAddress("x1", &hx1);
    hits->SetBranchAddress("y1", &hy1);
    for (Long64_t i = 0; i < hits->GetEntries(); ++i)
    {
      hits->GetEntry(i);
      if (event != event_id)
      {
        continue;
      }
      TMarker* marker = new TMarker(0.5 * (hx0 + hx1), 0.5 * (hy0 + hy1), 6);
      marker->SetMarkerColor(kGray + 1);
      marker->SetMarkerSize(0.35);
      marker->Draw();
    }
  }

  TLegend* legend = new TLegend(0.14, 0.74, 0.42, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(emcal_band, "EMCal 93-113 cm", "f");
  legend->AddEntry((TObject*) 0, "red: #pi^{0} truth", "");
  legend->AddEntry((TObject*) 0, "blue: #gamma truth", "");
  legend->AddEntry((TObject*) 0, "magenta: e^{#pm} truth", "");
  legend->Draw();

  canvas->SaveAs(output_filename.c_str());
  // canvas->SaveAs(Form("event_%d_xy.pdf", event_id));
  // canvas->SaveAs(Form("event_%d_xy.png", event_id));
  segments->ResetBranchAddresses();
  clusters->ResetBranchAddresses();
  if (draw_hits && hits)
  {
    hits->ResetBranchAddresses();
  }
  file->Close();
}

double xy_marker_size(double energy)
{
  if (!std::isfinite(energy) || energy < 0.0)
  {
    return 0.8;
  }
  return std::min(2.4, 0.8 + 0.18 * energy);
}


