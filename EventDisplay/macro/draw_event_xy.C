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


double xy_marker_size(double energy);
double xy_estimate_emcal_radius(TTree* clusters);

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

  const double emcal_r = xy_estimate_emcal_radius(clusters);
  const double frame_r = std::max(130.0, emcal_r + 30.0);
  gROOT->cd();
  TCanvas* canvas = new TCanvas("canvas_xy", "event display x-y", 900, 900);
  canvas->SetLeftMargin(0.14);
  canvas->SetRightMargin(0.04);
  canvas->SetTopMargin(0.06);
  canvas->SetBottomMargin(0.12);
  TH2F* frame = new TH2F("frame_xy", Form("Event %d;x [cm];y [cm]", event_id), 100, -frame_r, frame_r, 100, -frame_r, frame_r);
  frame->SetDirectory(nullptr);
  frame->Draw();

  TEllipse* emcal = new TEllipse(0.0, 0.0, emcal_r, emcal_r);
  emcal->SetFillStyle(0);
  emcal->SetLineColor(kGray + 2);
  emcal->SetLineStyle(2);
  emcal->SetLineWidth(2);
  emcal->Draw();

  int event = 0;
  int segment_type = 0;
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
  segments->SetBranchAddress("event", &event);
  segments->SetBranchAddress("segment_type", &segment_type);
  segments->SetBranchAddress("x0", &x0);
  segments->SetBranchAddress("y0", &y0);
  segments->SetBranchAddress("x1", &x1);
  segments->SetBranchAddress("y1", &y1);
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
    TLine* line = new TLine(x0, y0, x1, y1);
    line->SetLineWidth(segment_type == 1 ? 3 : 2);
    line->SetLineColor(segment_type == 1 ? kRed + 1 : (segment_type == 2 ? kBlue + 1 : kGray + 1));
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

  TLegend* legend = new TLegend(0.14, 0.78, 0.38, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(emcal, "EMCal radius", "l");
  legend->AddEntry((TObject*) 0, "red: #pi^{0} truth", "");
  legend->AddEntry((TObject*) 0, "blue: #gamma truth", "");
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


double xy_estimate_emcal_radius(TTree* clusters)
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

