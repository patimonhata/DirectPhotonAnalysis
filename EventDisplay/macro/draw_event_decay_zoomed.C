#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

struct DecayZoomParticle
{
  int track_id = 0;
  int pid = 0;
  int parent_id = 0;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double e = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
  double pt = 0.0;
  double p = 0.0;
  double eta = 0.0;
  double phi = 0.0;
};

struct DecayZoomSegment
{
  int track_id = 0;
  int pid = 0;
  int parent_id = 0;
  int segment_type = 0;
  double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
  double z1 = 0.0;
};

std::string decay_zoom_pid_label(int pid)
{
  switch (pid)
  {
  case 22:
    return "#gamma";
  case 111:
    return "#pi^{0}";
  case 11:
    return "e^{-}";
  case -11:
    return "e^{+}";
  case 211:
    return "#pi^{+}";
  case -211:
    return "#pi^{-}";
  default:
    return Form("pid %d", pid);
  }
}

bool decay_zoom_same_track(const std::vector<DecayZoomParticle>& particles, int track_id)
{
  for (const auto& particle : particles)
  {
    if (particle.track_id == track_id)
    {
      return true;
    }
  }
  return false;
}

std::vector<DecayZoomParticle> decay_zoom_daughters(const std::vector<DecayZoomParticle>& particles, int parent_id)
{
  std::vector<DecayZoomParticle> daughters;
  for (const auto& particle : particles)
  {
    if (particle.parent_id == parent_id)
    {
      daughters.push_back(particle);
    }
  }
  std::sort(daughters.begin(), daughters.end(), [](const DecayZoomParticle& lhs, const DecayZoomParticle& rhs) {
    return lhs.track_id < rhs.track_id;
  });
  return daughters;
}

std::string decay_zoom_pattern(const std::vector<DecayZoomParticle>& daughters)
{
  std::string pattern = "#pi^{0} #rightarrow";
  if (daughters.empty())
  {
    return pattern + " unknown";
  }
  for (const auto& daughter : daughters)
  {
    pattern += " " + decay_zoom_pid_label(daughter.pid);
  }
  return pattern;
}

bool decay_zoom_is_drawn_daughter(int pid)
{
  return pid == 22 || std::abs(pid) == 11;
}

int decay_zoom_daughter_color(int pid)
{
  if (pid == 22)
  {
    return kBlue + 1;
  }
  if (std::abs(pid) == 11)
  {
    return kMagenta + 1;
  }
  return kGray + 1;
}

bool decay_zoom_has_undrawn_daughters(const std::vector<DecayZoomParticle>& daughters)
{
  for (const auto& daughter : daughters)
  {
    if (!decay_zoom_is_drawn_daughter(daughter.pid))
    {
      return true;
    }
  }
  return false;
}

void draw_event_decay_zoomed(const char* filename = "event_display.root", int event_id = 0, double zoom_cm = 5.0)
{
  const std::string output_filename = Form("/sphenix/u/ryotaro/DirectPhotonAnalysis/EventDisplay/output/image/event_%d_decay_zoomed.pdf", event_id);
  gStyle->SetOptStat(0);

  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie())
  {
    std::cout << "draw_event_decay_zoomed - cannot open " << filename << std::endl;
    return;
  }

  TTree* particles_tree = static_cast<TTree*>(file->Get("truth_particles"));
  TTree* segments_tree = static_cast<TTree*>(file->Get("truth_segments"));
  if (!particles_tree || !segments_tree)
  {
    std::cout << "draw_event_decay_zoomed - missing truth_particles or truth_segments tree" << std::endl;
    file->Close();
    return;
  }

  std::vector<DecayZoomParticle> event_particles;
  std::vector<DecayZoomParticle> pi0s;
  int event = 0;
  DecayZoomParticle particle;
  particles_tree->SetBranchAddress("event", &event);
  particles_tree->SetBranchAddress("track_id", &particle.track_id);
  particles_tree->SetBranchAddress("pid", &particle.pid);
  particles_tree->SetBranchAddress("parent_id", &particle.parent_id);
  particles_tree->SetBranchAddress("px", &particle.px);
  particles_tree->SetBranchAddress("py", &particle.py);
  particles_tree->SetBranchAddress("pz", &particle.pz);
  particles_tree->SetBranchAddress("e", &particle.e);
  particles_tree->SetBranchAddress("vx", &particle.vx);
  particles_tree->SetBranchAddress("vy", &particle.vy);
  particles_tree->SetBranchAddress("vz", &particle.vz);
  particles_tree->SetBranchAddress("pt", &particle.pt);
  particles_tree->SetBranchAddress("p", &particle.p);
  particles_tree->SetBranchAddress("eta", &particle.eta);
  particles_tree->SetBranchAddress("phi", &particle.phi);
  for (Long64_t i = 0; i < particles_tree->GetEntries(); ++i)
  {
    particles_tree->GetEntry(i);
    if (event != event_id)
    {
      continue;
    }
    event_particles.push_back(particle);
    if (particle.pid == 111)
    {
      pi0s.push_back(particle);
    }
  }
  particles_tree->ResetBranchAddresses();

  std::vector<DecayZoomSegment> segments;
  DecayZoomSegment segment;
  segments_tree->SetBranchAddress("event", &event);
  segments_tree->SetBranchAddress("track_id", &segment.track_id);
  segments_tree->SetBranchAddress("pid", &segment.pid);
  segments_tree->SetBranchAddress("parent_id", &segment.parent_id);
  segments_tree->SetBranchAddress("x0", &segment.x0);
  segments_tree->SetBranchAddress("y0", &segment.y0);
  segments_tree->SetBranchAddress("z0", &segment.z0);
  segments_tree->SetBranchAddress("x1", &segment.x1);
  segments_tree->SetBranchAddress("y1", &segment.y1);
  segments_tree->SetBranchAddress("z1", &segment.z1);
  segments_tree->SetBranchAddress("segment_type", &segment.segment_type);
  for (Long64_t i = 0; i < segments_tree->GetEntries(); ++i)
  {
    segments_tree->GetEntry(i);
    if (event == event_id && (segment.segment_type == 1 || decay_zoom_same_track(pi0s, segment.parent_id)))
    {
      segments.push_back(segment);
    }
  }
  segments_tree->ResetBranchAddresses();

  gROOT->cd();
  TCanvas* canvas = new TCanvas("canvas_decay_zoomed", "pi0 decay zoom", 1200, 800);
  TPad* xy_pad = new TPad("decay_zoom_xy", "decay zoom x-y", 0.0, 0.5, 0.62, 1.0);
  TPad* zr_pad = new TPad("decay_zoom_zr", "decay zoom z-r", 0.0, 0.0, 0.62, 0.5);
  TPad* text_pad = new TPad("decay_zoom_text", "decay zoom text", 0.62, 0.0, 1.0, 1.0);
  xy_pad->Draw();
  zr_pad->Draw();
  text_pad->Draw();

  xy_pad->cd();
  xy_pad->SetLeftMargin(0.15);
  xy_pad->SetRightMargin(0.04);
  xy_pad->SetTopMargin(0.08);
  xy_pad->SetBottomMargin(0.12);
  const double range = std::max(0.1, zoom_cm);
  TH2F* frame = new TH2F("frame_decay_zoomed_xy", Form("Event %d #pi^{0} decay zoom x-y;x [cm];y [cm]", event_id), 100, -range, range, 100, -range, range);
  frame->SetDirectory(nullptr);
  frame->Draw();

  TMarker* origin = new TMarker(0.0, 0.0, 2);
  origin->SetMarkerColor(kGray + 2);
  origin->SetMarkerSize(1.2);
  origin->Draw();

  for (const auto& segment_to_draw : segments)
  {
    const bool draw_pi0 = segment_to_draw.segment_type == 1 && decay_zoom_same_track(pi0s, segment_to_draw.track_id);
    const bool draw_daughter = decay_zoom_same_track(pi0s, segment_to_draw.parent_id) && decay_zoom_is_drawn_daughter(segment_to_draw.pid);
    if (!draw_pi0 && !draw_daughter)
    {
      continue;
    }
    TLine* line = new TLine(segment_to_draw.x0, segment_to_draw.y0, segment_to_draw.x1, segment_to_draw.y1);
    line->SetLineColor(draw_pi0 ? kRed + 1 : decay_zoom_daughter_color(segment_to_draw.pid));
    line->SetLineWidth(draw_pi0 ? 3 : 2);
    line->Draw();

    if (draw_pi0)
    {
      TMarker* decay_vertex = new TMarker(segment_to_draw.x1, segment_to_draw.y1, 29);
      decay_vertex->SetMarkerColor(kBlack);
      decay_vertex->SetMarkerSize(1.3);
      decay_vertex->Draw();
    }
  }

  TLatex plot_latex;
  plot_latex.SetTextSize(0.030);
  plot_latex.SetTextAlign(12);
  plot_latex.DrawLatex(-0.92 * range, 0.88 * range, "red: #pi^{0}, blue: #gamma, magenta: e^{#pm}");

  zr_pad->cd();
  zr_pad->SetLeftMargin(0.15);
  zr_pad->SetRightMargin(0.04);
  zr_pad->SetTopMargin(0.10);
  zr_pad->SetBottomMargin(0.16);
  TH2F* frame_zr = new TH2F("frame_decay_zoomed_zr", Form("Event %d #pi^{0} decay zoom z-r;z [cm];r [cm]", event_id), 100, -range, range, 100, 0.0, range);
  frame_zr->SetDirectory(nullptr);
  frame_zr->Draw();

  TMarker* origin_zr = new TMarker(0.0, 0.0, 2);
  origin_zr->SetMarkerColor(kGray + 2);
  origin_zr->SetMarkerSize(1.2);
  origin_zr->Draw();

  for (const auto& segment_to_draw : segments)
  {
    const bool draw_pi0 = segment_to_draw.segment_type == 1 && decay_zoom_same_track(pi0s, segment_to_draw.track_id);
    const bool draw_daughter = decay_zoom_same_track(pi0s, segment_to_draw.parent_id) && decay_zoom_is_drawn_daughter(segment_to_draw.pid);
    if (!draw_pi0 && !draw_daughter)
    {
      continue;
    }
    const double r0 = std::sqrt(segment_to_draw.x0 * segment_to_draw.x0 + segment_to_draw.y0 * segment_to_draw.y0);
    const double r1 = std::sqrt(segment_to_draw.x1 * segment_to_draw.x1 + segment_to_draw.y1 * segment_to_draw.y1);
    TLine* line = new TLine(segment_to_draw.z0, r0, segment_to_draw.z1, r1);
    line->SetLineColor(draw_pi0 ? kRed + 1 : decay_zoom_daughter_color(segment_to_draw.pid));
    line->SetLineWidth(draw_pi0 ? 3 : 2);
    line->Draw();

    if (draw_pi0)
    {
      TMarker* decay_vertex = new TMarker(segment_to_draw.z1, r1, 29);
      decay_vertex->SetMarkerColor(kBlack);
      decay_vertex->SetMarkerSize(1.3);
      decay_vertex->Draw();
    }
  }

  TLatex zr_latex;
  zr_latex.SetTextSize(0.030);
  zr_latex.SetTextAlign(12);
  zr_latex.DrawLatex(-0.92 * range, 0.88 * range, "red: #pi^{0}, blue: #gamma, magenta: e^{#pm}");

  text_pad->cd();
  TLatex text;
  text.SetNDC(true);
  text.SetTextAlign(12);
  text.SetTextSize(0.035);
  double y = 0.94;
  text.DrawLatex(0.06, y, Form("Event %d decay summary", event_id));
  y -= 0.055;

  if (pi0s.empty())
  {
    text.SetTextSize(0.030);
    text.DrawLatex(0.06, y, "No #pi^{0} found in truth_particles.");
  }
  else
  {
    text.SetTextSize(0.028);
    for (const auto& pi0 : pi0s)
    {
      const std::vector<DecayZoomParticle> daughters = decay_zoom_daughters(event_particles, pi0.track_id);
      text.DrawLatex(0.06, y, decay_zoom_pattern(daughters).c_str());
      y -= 0.045;
      if (decay_zoom_has_undrawn_daughters(daughters))
      {
        text.DrawLatex(0.06, y, "Note: not all daughters are drawn; only #gamma/e^{#pm} are shown.");
        y -= 0.040;
      }
      text.DrawLatex(0.06, y, Form("#pi^{0} track %d: p = %.3f GeV, #eta = %.3f, #phi = %.3f", pi0.track_id, pi0.p, pi0.eta, pi0.phi));
      y -= 0.040;
      text.DrawLatex(0.06, y, Form("production: (x,y,z) = (%.4g, %.4g, %.4g) cm", pi0.vx, pi0.vy, pi0.vz));
      y -= 0.040;

      if (!daughters.empty())
      {
        const auto& first_daughter = daughters.front();
        text.DrawLatex(0.06, y, Form("decay:      (x,y,z) = (%.4g, %.4g, %.4g) cm", first_daughter.vx, first_daughter.vy, first_daughter.vz));
      }
      else
      {
        text.DrawLatex(0.06, y, "decay:      no direct daughter vertex found");
      }
      y -= 0.050;
      text.DrawLatex(0.06, y, "daughters:");
      y -= 0.038;

      for (const auto& daughter : daughters)
      {
        if (y < 0.08)
        {
          text.DrawLatex(0.06, y, "...");
          break;
        }
        text.DrawLatex(0.09, y, Form("%s track %d: E = %.3f GeV, #eta = %.3f, #phi = %.3f",
                                      decay_zoom_pid_label(daughter.pid).c_str(), daughter.track_id, daughter.e, daughter.eta, daughter.phi));
        y -= 0.038;
      }
      y -= 0.030;
      if (y < 0.08)
      {
        break;
      }
    }
  }

  canvas->SaveAs(output_filename.c_str());
  file->Close();
}
