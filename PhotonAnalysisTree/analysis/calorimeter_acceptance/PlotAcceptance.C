#include <TBox.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <iostream>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace calo_acceptance_plot
{
struct Surface
{
  std::string detector;
  double radius, etaLow, etaHigh, zLow, zHigh;
  const char* label() const { return detector == "CEMC" ? "EMCAL" : (detector == "HCALIN" ? "iHCAL" : "oHCAL"); }
  double low(double vertex) const { return std::asinh((zLow - vertex) / radius); }
  double high(double vertex) const { return std::asinh((zHigh - vertex) / radius); }
};

std::vector<Surface> read_surfaces(const std::string& path)
{
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot read " + path);
  std::string row;
  std::getline(input, row);
  std::vector<Surface> surfaces;
  while (std::getline(input, row))
  {
    std::replace(row.begin(), row.end(), ',', ' ');
    std::istringstream values(row);
    Surface s;
    if (!(values >> s.detector >> s.radius >> s.etaLow >> s.etaHigh >> s.zLow >> s.zHigh)) throw std::runtime_error("Invalid geometry CSV row");
    if (!std::isfinite(s.radius) || s.radius <= 0. || !std::isfinite(s.zLow) || !std::isfinite(s.zHigh) || s.zLow >= s.zHigh)
      throw std::runtime_error("Invalid representative surface");
    if (std::abs(s.low(0) - s.etaLow) > 1e-10 || std::abs(s.high(0) - s.etaHigh) > 1e-10)
      throw std::runtime_error("Geometry eta / z boundaries disagree");
    surfaces.push_back(s);
  }
  if (surfaces.size() != 3 || surfaces[0].detector != "CEMC" || surfaces[1].detector != "HCALIN" || surfaces[2].detector != "HCALOUT")
    throw std::runtime_error("Expected CEMC, HCALIN, HCALOUT surfaces");
  return surfaces;
}

void text(double x, double y, const std::string& caption, double size = 0.04)
{
  TLatex label; label.SetNDC(); label.SetTextFont(42); label.SetTextSize(size); label.DrawLatex(x, y, caption.c_str());
}

void frame(double x1, double y1, double x2, double y2, const char* title)
{
  gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.16);
  auto* axes = gPad->DrawFrame(x1, y1, x2, y2, title);
  axes->GetXaxis()->SetTitleSize(0.05); axes->GetYaxis()->SetTitleSize(0.05);
  axes->GetXaxis()->SetLabelSize(0.043); axes->GetYaxis()->SetLabelSize(0.043);
  axes->GetYaxis()->SetTitleOffset(1.1);
}

TGraph boundary(const Surface& s, bool upper)
{
  TGraph graph(201);
  for (int i = 0; i <= 200; ++i) graph.SetPoint(i, i - 100., upper ? s.high(i - 100.) : s.low(i - 100.));
  return graph;
}

void write_values(const std::string& path, const std::vector<Surface>& surfaces, bool selected)
{
  std::ofstream csv(path);
  if (!csv) throw std::runtime_error("Cannot write " + path);
  csv << std::setprecision(17) << "detector,z_vertex_cm,eta_min,eta_max,delta_eta\n";
  for (const auto& s : surfaces)
    for (int v = selected ? 0 : -100; v <= 100; v += selected ? 10 : 1)
      csv << s.detector << ',' << v << ',' << s.low(v) << ',' << s.high(v) << ',' << s.high(v) - s.low(v) << '\n';
}

void plot(const std::string& output)
{
  const auto surfaces = read_surfaces(output + "/geometry.csv");
  gStyle->SetOptStat(0); gStyle->SetTitleFont(42, "");
  double etaExtent = 0., zExtent = 100., rExtent = 0.;
  for (const auto& s : surfaces)
  {
    etaExtent = std::max({etaExtent, std::abs(s.low(100)), std::abs(s.high(-100))});
    zExtent = std::max({zExtent, std::abs(s.zLow), std::abs(s.zHigh)});
    rExtent = std::max(rExtent, s.radius);
  }
  etaExtent = std::ceil((etaExtent + 0.10) * 10.) / 10.;
  zExtent = std::ceil((zExtent + 20.) / 50.) * 50.;
  rExtent = std::ceil((rExtent + 20.) / 50.) * 50.;
  const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2};
  TCanvas scan("acceptance_scan", "Calorimeter acceptance", 1650, 650);
  scan.Divide(3, 1);
  for (std::size_t d = 0; d < surfaces.size(); ++d)
  {
    const auto& s = surfaces[d]; scan.cd(d + 1);
    frame(-100, -etaExtent, 100, etaExtent, ";z_{vertex} [cm];#eta (relative to vertex)");
    TGraph fill(402);
    for (int i = 0; i <= 200; ++i)
    {
      fill.SetPoint(i, i - 100., s.low(i - 100.));
      fill.SetPoint(201 + i, 100. - i, s.high(100. - i));
    }
    fill.SetFillColor(kGray); fill.SetLineColor(kGray + 2); fill.DrawClone("F");
    auto low = boundary(s, false), high = boundary(s, true); low.DrawClone("L"); high.DrawClone("L");
    text(0.14, 0.93, std::string(s.label()) + Form("   R = %.3f cm", s.radius), 0.045);
    text(0.14, 0.87, "ana.565 | tower representative surface", 0.032);
    gPad->RedrawAxis();
  }
  scan.Print((output + "/acceptance_vs_z.pdf").c_str());
  scan.Print((output + "/acceptance_vs_z.png").c_str());

  TCanvas compare("acceptance_comparison", "Calorimeter comparison", 1350, 650);
  compare.Divide(2, 1);
  compare.cd(1); frame(-100, -etaExtent, 100, etaExtent, ";z_{vertex} [cm];#eta boundary");
  TLegend legend(0.15, 0.68, 0.45, 0.83); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.035);
  for (std::size_t d = 0; d < surfaces.size(); ++d)
  {
    auto low = boundary(surfaces[d], false), high = boundary(surfaces[d], true);
    low.SetLineColor(colors[d]); high.SetLineColor(colors[d]); low.SetLineWidth(2); high.SetLineWidth(2); low.SetLineStyle(2);
    auto* drawn = high.DrawClone("L"); low.DrawClone("L"); legend.AddEntry(drawn, surfaces[d].label(), "l");
  }
  legend.DrawClone(); text(0.14, 0.93, "Boundary comparison", 0.045);
  text(0.14, 0.87, "solid: upper / dashed: lower", 0.035);
  compare.cd(2); frame(-100, 0, 100, 2.5, ";z_{vertex} [cm];#Delta#eta = #eta_{max} - #eta_{min}");
  for (std::size_t d = 0; d < surfaces.size(); ++d)
  {
    TGraph width(201);
    for (int i = 0; i <= 200; ++i) width.SetPoint(i, i - 100., surfaces[d].high(i - 100.) - surfaces[d].low(i - 100.));
    width.SetLineColor(colors[d]); width.SetLineWidth(2); width.DrawClone("L");
  }
  legend.DrawClone(); text(0.14, 0.93, "Accepted pseudorapidity width", 0.045);
  text(0.14, 0.87, "ana.565 | full #phi, x_{v} = y_{v} = 0", 0.035);
  compare.Print((output + "/detector_comparison.pdf").c_str());

  TCanvas display("acceptance_vertices", "Vertex acceptance", 1600, 1350);
  const std::string book = output + "/acceptance_vertex_views.pdf";
  display.Print((book + "[").c_str());
  for (int vertex = 0; vertex <= 100; vertex += 10)
  {
    display.Clear(); display.Divide(2, 3);
    for (std::size_t d = 0; d < surfaces.size(); ++d)
    {
      const auto& s = surfaces[d]; display.cd(2 * d + 1);
      frame(-zExtent, 0, zExtent, rExtent, ";z [cm];r [cm]");
      TBox surface(s.zLow, s.radius - 2., s.zHigh, s.radius + 2.); surface.SetFillColor(kGray + 1); surface.DrawClone();
      TLine rayLow(vertex, 0, s.zLow, s.radius), rayHigh(vertex, 0, s.zHigh, s.radius);
      rayLow.SetLineColor(colors[d]); rayHigh.SetLineColor(colors[d]); rayLow.SetLineWidth(2); rayHigh.SetLineWidth(2);
      rayLow.DrawClone(); rayHigh.DrawClone();
      TMarker point(vertex, 0, 20); point.SetMarkerColor(kRed + 1); point.SetMarkerSize(1.3); point.DrawClone();
      text(0.14, 0.93, std::string(s.label()) + Form("   z_{vertex} = %d cm", vertex), 0.052);
      text(0.14, 0.87, "Gray: tower surface (ribbon width for display only)", 0.035);
      display.cd(2 * d + 2);
      frame(-etaExtent, -TMath::Pi(), etaExtent, TMath::Pi(), ";#eta (relative to vertex);#phi [rad]");
      TBox acceptance(s.low(vertex), -TMath::Pi(), s.high(vertex), TMath::Pi());
      acceptance.SetFillColor(kGray); acceptance.SetLineColor(kGray + 2); acceptance.DrawClone();
      text(0.14, 0.93, Form("#eta #in [%.3f, %.3f]   #Delta#eta = %.3f", s.low(vertex), s.high(vertex), s.high(vertex) - s.low(vertex)), 0.047);
      text(0.14, 0.87, "ana.565 | full #phi | x_{v} = y_{v} = 0 | gaps ignored", 0.035);
      gPad->RedrawAxis();
    }
    display.Print(book.c_str(), Form("Title:zv %d cm", vertex));
    if (vertex == 0 || vertex == 100) display.Print((output + Form("/vertex_%03d.png", vertex)).c_str());
  }
  display.Print((book + "]").c_str());
  write_values(output + "/acceptance_scan.csv", surfaces, false);
  write_values(output + "/acceptance_selected_vertices.csv", surfaces, true);
}
}  // namespace calo_acceptance_plot
