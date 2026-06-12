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

namespace
{
double marker_size(double energy)
{
  if (!std::isfinite(energy) || energy < 0.0)
  {
    return 0.8;
  }
  return std::min(2.4, 0.8 + 0.18 * energy);
}

double estimate_emcal_radius(TTree* clusters)
{
  if (!clusters)
  {
    return 95.0;
  }

  int event = 0;
  double r = 0.0;
  double radius = 95.0;
  clusters->ResetBranchAddresses();
  clusters->SetBranchAddress("event", &event);
  clusters->SetBranchAddress("r", &r);

  for (Long64_t i = 0; i < clusters->GetEntries(); ++i)
  {
    clusters->GetEntry(i);
    if (std::isfinite(r) && r > 0.0)
    {
      radius = r;
      break;
    }
  }
  clusters->ResetBranchAddresses();
  return radius;
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

  const double emcal_r = estimate_emcal_radius(clusters);
  gROOT->cd();
  TCanvas* canvas = new TCanvas("canvas_zr", "event display z-r", 1100, 800);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(0.04);
  canvas->SetTopMargin(0.07);
  canvas->SetBottomMargin(0.12);
  TH2F* frame = new TH2F("frame_zr", Form("Event %d;z [cm];r [cm]", event_id), 100, -180.0, 180.0, 100, 0.0, 140.0);
  frame->SetDirectory(nullptr);
  frame->Draw();

  TLine* emcal_line = new TLine(-180.0, emcal_r, 180.0, emcal_r);
  emcal_line->SetLineColor(kGray + 2);
  emcal_line->SetLineStyle(2);
  emcal_line->SetLineWidth(2);
  emcal_line->Draw();

  int event = 0;
  int segment_type = 0;
  double z0 = 0.0;
  double z1 = 0.0;
  double r0 = 0.0;
  double r1 = 0.0;
  segments->SetBranchAddress("event", &event);
  segments->SetBranchAddress("segment_type", &segment_type);
  segments->SetBranchAddress("z0", &z0);
  segments->SetBranchAddress("z1", &z1);
  segments->SetBranchAddress("r0", &r0);
  segments->SetBranchAddress("r1", &r1);
  for (Long64_t i = 0; i < segments->GetEntries(); ++i)
  {
    segments->GetEntry(i);
    if (event != event_id)
    {
      continue;
    }
    if (!draw_other_truth && segment_type != 1 && segment_type != 2)
    {
      continue;
    }
    TLine* line = new TLine(z0, r0, z1, r1);
    line->SetLineWidth(segment_type == 1 ? 3 : 2);
    line->SetLineColor(segment_type == 1 ? kRed + 1 : (segment_type == 2 ? kBlue + 1 : kGray + 1));
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

  TLegend* legend = new TLegend(0.14, 0.75, 0.38, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(emcal_line, "EMCal radius", "l");
  legend->AddEntry((TObject*) 0, "red: #pi^{0} truth", "");
  legend->AddEntry((TObject*) 0, "blue: #gamma truth", "");
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
