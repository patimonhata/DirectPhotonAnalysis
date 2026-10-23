#ifndef CHECK_ETA_DECAY_FLIGHT_C
#define CHECK_ETA_DECAY_FLIGHT_C

#include "../../../analysis/common/plotting/sPhenixStyle.C"
#include <frog/FROG.h>
#include <g4main/PHG4MCProcessDefs.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>
#include <HepMC/WeightContainer.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libphhepmc.so)
R__LOAD_LIBRARY(libFROG.so)

namespace eta_decay_flight
{
constexpr double nan = std::numeric_limits<double>::quiet_NaN();
constexpr std::array<double, 4> thresholds_cm = {1e-4, 1e-3, 1e-2, 1e-1}; // 1 um, 10 um, 0.1 mm, 1 mm; diagnostics, never origin cuts.
constexpr std::array<const char*, 11> status_names = {
    "confirmed_decay", "inferred_decay_no_process", "missing_production_vertex", "no_saved_daughters", "missing_daughter_vertex",
    "ambiguous_decay_vertices", "non_decay_interaction_only", "unknown_process_or_decay", "invalid_coordinates", "negative_delta_t",
    "decay_with_other_vertices"};
constexpr std::array<const char*, 6> mode_names = {"unresolved", "gamma_gamma", "three_pi0", "pi_plus_pi_minus_pi0", "pi_plus_pi_minus_gamma", "other"};
constexpr std::array<double, 6> pt_edges = {0, 5, 10, 20, 40, 80};

struct Record
{
  int input_file_index = -1, source_file_id = -1, eta_track_id = 0, eta_barcode = 0, embedding_id = 0;
  Long64_t event_in_file = -1;
  int hepmc_event_number = -1, production_vertex_id = 0, decay_vertex_id = 0, decay_process = -1;
  int status = 3, mode = 0, n_direct_daughters = 0, n_vertex_groups = 0, n_missing_vertices = 0;
  bool distance_valid = false, usable = false, time_valid = false, same_vertex_id = false, zero_displacement = false;
  bool weight_valid = true, hepmc_weight_available = false;
  double weight = 1, energy = nan, pt = nan, eta = nan, px = nan, py = nan, pz = nan;
  double production_x = nan, production_y = nan, production_z = nan, production_t = nan;
  double decay_x = nan, decay_y = nan, decay_z = nan, decay_t = nan;
  double dx = nan, dy = nan, dz = nan, delta_t = nan, length_cm = nan, transverse_cm = nan, parallel_cm = nan, perpendicular_cm = nan, alignment = nan;
  std::vector<int> daughter_track_ids, daughter_pdgs, daughter_vertex_ids, daughter_processes;
};

inline bool finite_position(const PHG4VtxPoint* v)
{
  return v && std::isfinite(v->get_x()) && std::isfinite(v->get_y()) && std::isfinite(v->get_z());
}
inline int decay_mode(std::vector<int> pdgs)
{
  std::sort(pdgs.begin(), pdgs.end());
  if (pdgs == std::vector<int>{22, 22}) return 1;
  if (pdgs == std::vector<int>{111, 111, 111}) return 2;
  if (pdgs == std::vector<int>{-211, 111, 211}) return 3;
  if (pdgs == std::vector<int>{-211, 22, 211}) return 4;
  return 5;
}
inline bool unknown_process(int process) { return process < 0 || process == PHG4MCProcess::kPNoProcess; }

struct VertexGroup
{
  const PHG4VtxPoint* vertex = nullptr;
  std::vector<int> ids, pdgs;
};
// Direct daughter lookup is built once per event, including saved tracks with unusual sign conventions.
using Children = std::map<int, std::vector<const PHG4Particle*>>;
inline Children children_by_parent(const PHG4TruthInfoContainer& truth)
{
  Children result;
  for (const auto& [id, particle] : truth.GetMap())
    if (particle && particle->get_parent_id() != 0) result[particle->get_parent_id()].push_back(particle);
  return result;
}

inline Record analyze(const PHG4Particle& parent, PHG4TruthInfoContainer& truth, const Children& children)
{
  Record r;
  r.eta_track_id = parent.get_track_id(); r.eta_barcode = parent.get_barcode(); r.embedding_id = truth.isEmbeded(parent.get_track_id());
  r.energy = parent.get_e(); r.px = parent.get_px(); r.py = parent.get_py(); r.pz = parent.get_pz();
  r.pt = std::hypot(r.px, r.py);
  if (std::isfinite(r.pt) && r.pt > 0 && std::isfinite(r.pz)) r.eta = std::asinh(r.pz / r.pt);
  r.production_vertex_id = parent.get_vtx_id();
  const auto* production = truth.GetVtx(r.production_vertex_id);
  if (production)
  {
    r.production_x = production->get_x(); r.production_y = production->get_y(); r.production_z = production->get_z(); r.production_t = production->get_t();
  }
  std::vector<VertexGroup> groups;
  const auto found = children.find(parent.get_track_id());
  if (found != children.end()) for (const auto* child : found->second)
  {
    ++r.n_direct_daughters;
    r.daughter_track_ids.push_back(child->get_track_id()); r.daughter_pdgs.push_back(child->get_pid()); r.daughter_vertex_ids.push_back(child->get_vtx_id());
    const auto* vertex = truth.GetVtx(child->get_vtx_id());
    r.daughter_processes.push_back(vertex ? vertex->get_process() : -1);
    if (!vertex) { ++r.n_missing_vertices; continue; }
    auto group = std::find_if(groups.begin(), groups.end(), [&](const auto& g) {
      // Some productions assign separate IDs to exactly equal daughter vertices. Do not merge merely nearby vertices.
      const auto* v = g.vertex;
      return v == vertex || (finite_position(v) && finite_position(vertex) && v->get_x() == vertex->get_x() &&
          v->get_y() == vertex->get_y() && v->get_z() == vertex->get_z() && v->get_t() == vertex->get_t() && v->get_process() == vertex->get_process());
    });
    if (group == groups.end()) groups.push_back({vertex, {child->get_vtx_id()}, {child->get_pid()}});
    else { group->ids.push_back(child->get_vtx_id()); group->pdgs.push_back(child->get_pid()); }
  }
  r.n_vertex_groups = groups.size();
  if (!production) { r.status = 2; return r; }
  if (r.n_direct_daughters == 0) { r.status = 3; return r; }
  if (r.n_missing_vertices) { r.status = 4; return r; }
  std::vector<std::size_t> decays;
  for (std::size_t i = 0; i < groups.size(); ++i)
    if (groups[i].vertex->get_process() == PHG4MCProcess::kPDecay) decays.push_back(i);
  const VertexGroup* selected = nullptr;
  if (decays.size() > 1) { r.status = 5; return r; }
  if (decays.size() == 1)
  {
    selected = &groups[decays.front()];
    r.status = groups.size() == 1 ? 0 : 10;
  }
  else if (groups.size() == 1 && unknown_process(groups.front().vertex->get_process()) && decay_mode(groups.front().pdgs) != 5)
  {
    selected = &groups.front(); r.status = 1;
  }
  else
  {
    const bool explicit_non_decay = std::any_of(groups.begin(), groups.end(), [](const auto& g) { return !unknown_process(g.vertex->get_process()); });
    r.status = explicit_non_decay ? 6 : (groups.size() > 1 ? 5 : 7);
    return r;
  }
  const auto* decay = selected->vertex;
  r.mode = decay_mode(selected->pdgs); r.decay_vertex_id = selected->ids.front(); r.decay_process = decay->get_process();
  r.decay_x = decay->get_x(); r.decay_y = decay->get_y(); r.decay_z = decay->get_z(); r.decay_t = decay->get_t();
  r.same_vertex_id = std::find(selected->ids.begin(), selected->ids.end(), r.production_vertex_id) != selected->ids.end();
  if (!finite_position(production) || !finite_position(decay)) { r.status = 8; return r; }
  r.dx = r.decay_x - r.production_x; r.dy = r.decay_y - r.production_y; r.dz = r.decay_z - r.production_z;
  r.transverse_cm = std::hypot(r.dx, r.dy); r.length_cm = std::hypot(r.transverse_cm, r.dz);
  r.distance_valid = std::isfinite(r.length_cm); r.zero_displacement = r.length_cm == 0;
  if (!r.distance_valid) { r.status = 8; return r; }
  const double momentum = std::hypot(std::hypot(r.px, r.py), r.pz);
  if (std::isfinite(momentum) && momentum > 0)
  {
    r.parallel_cm = (r.dx * r.px + r.dy * r.py + r.dz * r.pz) / momentum;
    r.perpendicular_cm = std::sqrt(std::max(0.0, r.length_cm * r.length_cm - r.parallel_cm * r.parallel_cm));
    if (r.length_cm > 0) r.alignment = std::clamp(r.parallel_cm / r.length_cm, -1.0, 1.0);
  }
  if (std::isfinite(r.production_t) && std::isfinite(r.decay_t))
  {
    r.delta_t = r.decay_t - r.production_t; r.time_valid = true;
    if (r.delta_t < 0) { r.status = 9; return r; }
  }
  r.usable = r.status == 0 || r.status == 1;
  return r;
}

inline std::vector<std::string> csv_fields(const std::string& line)
{
  std::vector<std::string> fields; std::string value; bool quoted = false;
  for (std::size_t i = 0; i < line.size(); ++i)
  {
    if (line[i] == '"')
    {
      if (quoted && i + 1 < line.size() && line[i + 1] == '"') { value += '"'; ++i; }
      else quoted = !quoted;
    }
    else if (line[i] == ',' && !quoted) { fields.push_back(value); value.clear(); }
    else if (line[i] != '\r') value += line[i];
  }
  if (quoted) throw std::runtime_error("Unterminated CSV quote");
  fields.push_back(value); return fields;
}
inline std::string lower(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return std::tolower(c); }); return value;
}
inline std::string sample_from_path(const std::string& path)
{
  const auto begin = path.find("pythia8_"); if (begin == std::string::npos) return "";
  const auto end = path.find('-', begin); return lower(path.substr(begin + 8, end - begin - 8));
}
inline int segment_from_path(const std::string& path)
{
  const auto dot = path.rfind(".root"), dash = path.rfind('-');
  if (dot == std::string::npos || dash == std::string::npos || dot <= dash) return -1;
  try
  {
    const std::string token = path.substr(dash + 1, dot - dash - 1); std::size_t consumed = 0;
    const int value = std::stoi(token, &consumed); return consumed == token.size() ? value : -1;
  }
  catch (...) { return -1; }
}
struct Request
{
  std::string sample;
  int file = -1, track = 0, barcode = std::numeric_limits<int>::min(), hepmc_event = -1;
  Long64_t event = -1;
};
inline std::vector<Request> read_requests(const std::string& path)
{
  if (path.empty()) return {};
  std::vector<Request> requests;
  std::ifstream input(path); if (!input) throw std::runtime_error("Cannot read requests: " + path);
  std::string line; if (!std::getline(input, line)) throw std::runtime_error("Empty request CSV");
  const auto header = csv_fields(line);
  std::map<std::string, std::size_t> columns;
  for (std::size_t i = 0; i < header.size(); ++i) columns[header[i]] = i;
  const std::string track_column = columns.count("eta_track_id") ? "eta_track_id" : "g4_eta_track";
  for (const auto& name : {std::string("source_file_id"), std::string("event_in_file"), track_column})
    if (!columns.count(name)) throw std::runtime_error("Missing request column: " + name);
  while (std::getline(input, line))
  {
    if (line.empty()) continue;
    const auto fields = csv_fields(line); if (fields.size() != header.size()) throw std::runtime_error("Malformed request CSV row");
    Request r;
    r.file = std::stoi(fields[columns.at("source_file_id")]); r.event = std::stoll(fields[columns.at("event_in_file")]);
    r.track = std::stoi(fields[columns.at(track_column)]);
    if (columns.count("sample")) r.sample = lower(fields[columns.at("sample")]);
    if (columns.count("eta_barcode")) r.barcode = std::stoi(fields[columns.at("eta_barcode")]);
    if (columns.count("hepmc_event_number")) r.hepmc_event = std::stoi(fields[columns.at("hepmc_event_number")]);
    if (r.file < 0 || r.event < 0 || r.track <= 0) throw std::runtime_error("Invalid request file/event/primary track ID");
    requests.push_back(r);
  }
  if (requests.empty()) throw std::runtime_error("No requested eta tracks");
  return requests;
}
inline std::string quoted(const std::string& value)
{
  std::string result = "\""; for (char c : value) { result += c; if (c == '"') result += c; } return result + '"';
}

struct Observation { double length, transverse, weight; };
struct Summary
{
  Long64_t total = 0, usable = 0, confirmed = 0, inferred = 0, zero = 0, same_id = 0, missing_time = 0, invalid_weight = 0;
  std::array<Long64_t, status_names.size()> statuses = {};
  std::vector<Observation> distances;
  void add(const Record& r)
  {
    ++total; ++statuses.at(r.status);
    if (!r.weight_valid) ++invalid_weight;
    if (!r.usable) return;
    ++usable; if (r.status == 0) ++confirmed; else ++inferred;
    zero += r.zero_displacement; same_id += r.same_vertex_id; missing_time += !r.time_valid;
    distances.push_back({r.length_cm, r.transverse_cm, r.weight_valid ? r.weight : nan});
  }
};
inline std::vector<std::string> group_names(const Record& r)
{
  std::vector<std::string> names = {"all", std::string("mode_") + mode_names.at(r.mode)};
  if (std::isfinite(r.pt) && r.pt >= 0)
  {
    const auto upper = std::upper_bound(pt_edges.begin(), pt_edges.end(), r.pt); const int index = upper - pt_edges.begin() - 1;
    const std::string suffix = index + 1 == int(pt_edges.size()) ? "plus" : "to" + std::to_string(int(pt_edges[index + 1]));
    names.push_back("pt_" + std::to_string(int(pt_edges[index])) + suffix);
  }
  if (std::isfinite(r.energy) && r.energy >= 0)
  {
    const auto upper = std::upper_bound(pt_edges.begin(), pt_edges.end(), r.energy); const int index = upper - pt_edges.begin() - 1;
    const std::string suffix = index + 1 == int(pt_edges.size()) ? "plus" : "to" + std::to_string(int(pt_edges[index + 1]));
    names.push_back("energy_" + std::to_string(int(pt_edges[index])) + suffix);
  }
  return names;
}
struct Fraction { double value = nan, error = nan; };
inline Fraction fraction_statistics(const Summary& summary, double threshold, bool weighted = false)
{
  double numerator = 0, denominator = 0, numerator_w2 = 0, complement_w2 = 0;
  for (const auto& item : summary.distances)
  {
    const double w = weighted ? item.weight : 1.; if (!std::isfinite(w)) continue;
    denominator += w;
    if (item.length < threshold) { numerator += w; numerator_w2 += w * w; }
    else complement_w2 += w * w;
  }
  if (!(denominator > 0)) return {};
  const double fraction = numerator / denominator;
  const double variance = (std::pow(1 - fraction, 2) * numerator_w2 + fraction * fraction * complement_w2) / (denominator * denominator);
  return {fraction, std::sqrt(std::max(0., variance))};
}
inline double fraction_below(const Summary& summary, double threshold, bool weighted = false)
{
  return fraction_statistics(summary, threshold, weighted).value;
}
inline double quantile(const Summary& summary, double q)
{
  if (summary.distances.empty()) return nan;
  std::vector<double> values; for (const auto& item : summary.distances) values.push_back(item.length);
  std::sort(values.begin(), values.end()); return values[std::min(values.size() - 1, std::size_t(std::ceil(q * values.size()) - 1))];
}

inline void create_branches(TTree& tree, Record& r)
{
#define B(member) tree.Branch(#member, &r.member)
  B(input_file_index); B(source_file_id); B(event_in_file); B(hepmc_event_number); B(eta_track_id); B(eta_barcode); B(embedding_id);
  B(production_vertex_id); B(decay_vertex_id); B(decay_process); B(status); B(mode); B(n_direct_daughters); B(n_vertex_groups); B(n_missing_vertices);
  B(distance_valid); B(usable); B(time_valid); B(same_vertex_id); B(zero_displacement); B(weight_valid); B(hepmc_weight_available); B(weight);
  B(energy); B(pt); B(eta); B(px); B(py); B(pz); B(production_x); B(production_y); B(production_z); B(production_t);
  B(decay_x); B(decay_y); B(decay_z); B(decay_t); B(dx); B(dy); B(dz); B(delta_t); B(length_cm); B(transverse_cm);
  B(parallel_cm); B(perpendicular_cm); B(alignment); B(daughter_track_ids); B(daughter_pdgs); B(daughter_vertex_ids); B(daughter_processes);
#undef B
}

inline void write_summaries(TFile& output, const std::string& base, const std::map<std::string, Summary>& summaries)
{
  std::ofstream csv(base + "_summary.csv");
  csv << std::setprecision(17) << "group,total_eta,usable_eta,confirmed_decay,inferred_decay,zero_displacement,same_vertex_id,missing_time,invalid_weight,median_cm,p95_cm,p99_cm,max_cm";
  for (double cut : thresholds_cm) csv << ",fraction_L_lt_" << cut << "cm,weighted_fraction_L_lt_" << cut << "cm";
  for (const auto* name : status_names) csv << ",status_" << name;
  csv << '\n';
  std::ofstream cumulative(base + "_cumulative.csv"); cumulative << std::setprecision(17) << "group,threshold_cm,raw_fraction,raw_error,weighted_fraction,weighted_error\n";
  for (const auto& [name, s] : summaries)
  {
    csv << name << ',' << s.total << ',' << s.usable << ',' << s.confirmed << ',' << s.inferred << ',' << s.zero << ',' << s.same_id << ',' << s.missing_time << ',' << s.invalid_weight
        << ',' << quantile(s, .5) << ',' << quantile(s, .95) << ',' << quantile(s, .99) << ',' << quantile(s, 1);
    for (double cut : thresholds_cm) csv << ',' << fraction_below(s, cut) << ',' << fraction_below(s, cut, true);
    for (auto count : s.statuses) csv << ',' << count;
    csv << '\n';
    auto* directory = output.mkdir(name.c_str()); if (!directory) throw std::runtime_error("Cannot create ROOT summary group"); directory->cd();
    TH1D length("h_log10_length_cm", ";log_{10}(L / cm);#eta", 160, -14, 2);
    TH1D transverse("h_log10_transverse_cm", ";log_{10}(L_{T} / cm);#eta", 160, -14, 2);
    TH1D weighted("h_log10_length_weighted", ";log_{10}(L / cm);Weighted #eta", 160, -14, 2);
    TH1D cdf("h_cumulative_length", ";log_{10}(threshold / cm);Fraction with L < threshold", 160, -14, 2);
    TH1D weighted_cdf("h_cumulative_length_weighted", ";log_{10}(threshold / cm);Weighted fraction with L < threshold", 160, -14, 2);
    TH1D status("h_status", ";Diagnostic status;#eta", status_names.size(), -.5, status_names.size() - .5);
    for (auto* h : {&length, &transverse, &weighted, &cdf, &weighted_cdf, &status}) { h->SetDirectory(nullptr); h->Sumw2(); }
    for (const auto& item : s.distances)
    {
      // Exact zeros stay in the TTree and summary; they are never replaced with an arbitrary positive distance.
      if (item.length > 0) { length.Fill(std::log10(item.length)); if (std::isfinite(item.weight)) weighted.Fill(std::log10(item.length), item.weight); }
      if (item.transverse > 0) transverse.Fill(std::log10(item.transverse));
    }
    for (std::size_t i = 0; i < status_names.size(); ++i) { status.SetBinContent(i + 1, s.statuses[i]); status.SetBinError(i + 1, std::sqrt(double(s.statuses[i]))); status.GetXaxis()->SetBinLabel(i + 1, status_names[i]); }
    for (int bin = 1; bin <= cdf.GetNbinsX(); ++bin)
    {
      const double cut = std::pow(10., cdf.GetBinCenter(bin));
      const auto raw = fraction_statistics(s, cut), weight = fraction_statistics(s, cut, true);
      if (std::isfinite(raw.value)) { cdf.SetBinContent(bin, raw.value); cdf.SetBinError(bin, raw.error); }
      if (std::isfinite(weight.value)) { weighted_cdf.SetBinContent(bin, weight.value); weighted_cdf.SetBinError(bin, weight.error); }
      cumulative << name << ',' << cut << ',' << raw.value << ',' << raw.error << ',' << weight.value << ',' << weight.error << '\n';
    }
    length.Write(); transverse.Write(); weighted.Write(); cdf.Write(); weighted_cdf.Write(); status.Write();
  }
  if (!csv.good() || !cumulative.good()) throw std::runtime_error("Cannot write summary CSV");
  output.cd();
}

inline void draw_report(TFile& output, const std::string& base, const std::map<std::string, Summary>& summaries)
{
  SetsPhenixStyle(); TCanvas canvas("eta_flight_report", "Eta decay displacement", 1500, 1000);
  const std::string pdf = base + ".pdf"; canvas.Print((pdf + "[").c_str());
  for (const auto& [name, s] : summaries)
  {
    canvas.Clear(); canvas.Divide(2, 2);
    const auto draw = [&](int pad, const char* key, bool log_y) {
      canvas.cd(pad); gPad->SetLogy(log_y);
      auto* h = output.Get<TH1D>((name + "/" + key).c_str()); h->SetLineColor(kBlue + 1); h->SetLineWidth(3);
      if (std::string(key).find("cumulative") != std::string::npos) { h->SetMinimum(0); h->SetMaximum(1.1); }
      h->DrawCopy("HIST"); gPad->Modified(); gPad->Update();
    };
    draw(1, "h_log10_length_cm", s.usable > s.zero); draw(2, "h_log10_transverse_cm", false); draw(3, "h_cumulative_length", false);
    if (s.usable == s.zero)
    {
      canvas.cd(1); TLatex note; note.SetNDC(); note.SetTextSize(.05);
      note.DrawLatex(.20, .65, s.usable ? "All usable displacements are exactly zero" : "No usable decay vertices");
      note.DrawLatex(.20, .57, "Exact zeros are counted separately");
    }
    canvas.cd(4); TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(.045);
    double y = .92; const auto line = [&](const std::string& value) { text.DrawLatex(.04, y, value.c_str()); y -= .075; };
    line("#it{#bf{sPHENIX}} Internal: #eta decay flight"); line(name); line("All #eta: " + std::to_string(s.total) + "; usable: " + std::to_string(s.usable));
    line("Confirmed: " + std::to_string(s.confirmed) + "; inferred: " + std::to_string(s.inferred));
    line("Exact zero: " + std::to_string(s.zero) + "; same vertex ID: " + std::to_string(s.same_id));
    line("Unresolved / excluded: " + std::to_string(s.total - s.usable));
    for (double cut : thresholds_cm)
    {
      std::ostringstream message; message << std::setprecision(4) << "L < " << cut << " cm: ";
      const double fraction = fraction_below(s, cut); if (std::isfinite(fraction)) message << 100 * fraction << "%"; else message << "n/a";
      line(message.str());
    }
    line("Fractions use usable #eta, including exact zeros."); line("Zeros excluded from log-distance histograms.");
    canvas.Print(pdf.c_str());
    if (name == "all") canvas.SaveAs((base + "_overview.png").c_str());
  }
  canvas.Print((pdf + "]").c_str());
}

inline void enable_branch(TBranch& branch)
{
  branch.SetStatus(true);
  const auto* nested = branch.GetListOfBranches();
  for (int i = 0; i < nested->GetEntriesFast(); ++i) enable_branch(*static_cast<TBranch*>(nested->UncheckedAt(i)));
}

inline int run(const std::string& input_list, const std::string& base, Long64_t first_file, Long64_t n_files,
               Long64_t max_events_per_file, const std::string& requests_csv, int embedding_id, double sample_weight)
{
  if (input_list.empty() || base.empty() || first_file < 0 || n_files < 0 || max_events_per_file < 0 || !std::isfinite(sample_weight) || sample_weight <= 0)
    throw std::runtime_error("Invalid arguments");
  std::ifstream list(input_list); if (!list) throw std::runtime_error("Cannot read G4Hits list: " + input_list);
  std::vector<std::string> inputs; std::set<std::string> unique; std::string line;
  while (std::getline(list, line))
  {
    const auto begin = line.find_first_not_of(" \t\r"); if (begin == std::string::npos || line[begin] == '#') continue;
    line = line.substr(begin, line.find_last_not_of(" \t\r") - begin + 1);
    if (!unique.insert(line).second) throw std::runtime_error("Duplicate input: " + line);
    inputs.push_back(line);
  }
  if (first_file >= Long64_t(inputs.size())) throw std::runtime_error("Input file range is empty");
  const Long64_t end_file = n_files == 0 ? inputs.size() : std::min(Long64_t(inputs.size()), first_file + n_files);
  const auto requests = read_requests(requests_csv);
  for (const char* suffix : {".root", ".pdf", "_summary.csv", "_cumulative.csv", "_overview.png"})
    if (!gSystem->AccessPathName((base + suffix).c_str())) throw std::runtime_error("Refusing to overwrite " + base + suffix);
  const auto slash = base.find_last_of('/');
  if (slash != std::string::npos && gSystem->mkdir(base.substr(0, slash).c_str(), true) != 0 && gSystem->AccessPathName(base.substr(0, slash).c_str()))
    throw std::runtime_error("Cannot create output directory");
  TFile output((base + ".root").c_str(), "CREATE"); if (output.IsZombie()) throw std::runtime_error("Cannot create output ROOT file");
  Record record; TTree eta_tree("eta_tree", "One row per G4-primary eta; displacement is not an origin selection cut"); create_branches(eta_tree, record);
  int event_file_index = 0; Long64_t event_entry = 0; int event_eta = 0, event_usable = 0, event_unresolved = 0;
  std::array<int, 4> event_long_eta = {}; double event_weight = 1;
  TTree event_tree("event_tree", "Event-level counts retain events with unmeasurable eta decays");
  event_tree.Branch("input_file_index", &event_file_index); event_tree.Branch("event_in_file", &event_entry);
  event_tree.Branch("n_eta", &event_eta); event_tree.Branch("n_usable_eta", &event_usable); event_tree.Branch("n_unresolved_eta", &event_unresolved);
  event_tree.Branch("n_eta_ge_threshold", event_long_eta.data(), "n_eta_ge_threshold[4]/I"); event_tree.Branch("weight", &event_weight);
  TH2D length_pt("h_log10_length_vs_pt", ";#eta p_{T} [GeV];log_{10}(L / cm)", 100, 0, 200, 160, -14, 2);
  TH2D length_energy("h_log10_length_vs_energy", ";#eta energy [GeV];log_{10}(L / cm)", 100, 0, 200, 160, -14, 2);
  TH1D alignment("h_alignment", ";cos(angle between displacement and initial momentum);#eta", 100, -1, 1);
  length_pt.SetDirectory(nullptr); length_energy.SetDirectory(nullptr); alignment.SetDirectory(nullptr);
  std::map<std::string, Summary> summaries; summaries["all"];
  for (const auto* name : mode_names) summaries[std::string("mode_") + name];
  std::vector<std::string> resolved_inputs;
  Long64_t events_read = 0, events_with_eta = 0, requests_matched = 0;
  std::set<std::tuple<std::string, int, Long64_t, int>> matched;
  for (Long64_t index = first_file; index < end_file; ++index)
  {
    const auto& name = inputs[index]; const std::string sample = sample_from_path(name); const int segment = segment_from_path(name);
    std::map<Long64_t, std::map<int, Request>> selected;
    if (!requests_csv.empty())
    {
      for (const auto& r : requests) if (r.file == segment && (r.sample.empty() || r.sample == sample))
      {
        auto& targets = selected[r.event]; const auto found = targets.find(r.track);
        if (found != targets.end() && (found->second.barcode != r.barcode || found->second.hepmc_event != r.hepmc_event)) throw std::runtime_error("Conflicting duplicate requests");
        targets[r.track] = r;
      }
      if (selected.empty()) continue;
    }
    std::string path = name;
    if (gSystem->AccessPathName(path.c_str()))
    {
      const auto slash = input_list.find_last_of('/');
      const std::string relative = slash == std::string::npos ? name : input_list.substr(0, slash + 1) + name;
      if (!gSystem->AccessPathName(relative.c_str())) path = relative;
      else { FROG frog; path = frog.location(name.c_str()); }
    }
    if (std::find(resolved_inputs.begin(), resolved_inputs.end(), path) != resolved_inputs.end()) throw std::runtime_error("Same resolved DST listed more than once");
    resolved_inputs.push_back(path);
    std::unique_ptr<TFile> input(TFile::Open(path.c_str())); if (!input || input->IsZombie()) throw std::runtime_error("Cannot open DST: " + path);
    auto* tree = input->Get<TTree>("T"); if (!tree) throw std::runtime_error("Missing DST tree T: " + path);
    tree->SetBranchStatus("*", false);
    auto* truth_branch = tree->GetBranch("DST#G4TruthInfo");
    if (!truth_branch) throw std::runtime_error("Missing G4TruthInfo: " + path);
    enable_branch(*truth_branch);
    PHG4TruthInfoContainer* truth = nullptr; PHHepMCGenEventMap* hepmc = nullptr;
    if (tree->SetBranchAddress("DST#G4TruthInfo", &truth) < 0) throw std::runtime_error("Missing G4TruthInfo: " + path);
    const bool has_hepmc = tree->GetBranch("DST#PHHepMCGenEventMap");
    if (has_hepmc) { enable_branch(*tree->GetBranch("DST#PHHepMCGenEventMap")); if (tree->SetBranchAddress("DST#PHHepMCGenEventMap", &hepmc) < 0) throw std::runtime_error("Invalid HepMC branch"); }
    Long64_t scanned = 0;
    const Long64_t entries = tree->GetEntries();
    std::vector<Long64_t> entry_indices;
    if (requests_csv.empty()) for (Long64_t entry = 0; entry < entries && (max_events_per_file == 0 || entry < max_events_per_file); ++entry) entry_indices.push_back(entry);
    else for (const auto& [entry, targets] : selected)
    {
      if (entry >= entries) throw std::runtime_error("Requested event outside DST: " + path);
      if (max_events_per_file == 0 || Long64_t(entry_indices.size()) < max_events_per_file) entry_indices.push_back(entry);
    }
    for (Long64_t entry : entry_indices)
    {
      if (tree->GetEntry(entry) <= 0 || !truth) throw std::runtime_error("Cannot read DST event: " + path);
      ++events_read; ++scanned;
      const auto* subevent = hepmc ? hepmc->get(embedding_id) : nullptr;
      const auto* gen = subevent ? subevent->getEvent() : nullptr;
      const bool gen_weight_available = gen && !gen->weights().empty();
      event_weight = sample_weight * (gen_weight_available ? gen->weights().front() : 1.);
      event_file_index = index; event_entry = entry; event_eta = event_usable = event_unresolved = 0; event_long_eta.fill(0);
      const auto children = children_by_parent(*truth); const auto primaries = truth->GetPrimaryParticleRange();
      for (auto it = primaries.first; it != primaries.second; ++it)
      {
        const auto* parent = it->second; if (!parent || parent->get_pid() != 221 || truth->isEmbeded(parent->get_track_id()) != embedding_id) continue;
        const Request* target = nullptr;
        if (!requests_csv.empty())
        {
          const auto found = selected.at(entry).find(parent->get_track_id()); if (found == selected.at(entry).end()) continue;
          target = &found->second;
          if ((target->barcode != std::numeric_limits<int>::min() && target->barcode != parent->get_barcode()) ||
              (target->hepmc_event >= 0 && (!gen || gen->event_number() != target->hepmc_event))) throw std::runtime_error("Requested eta provenance mismatch");
          matched.insert({target->sample, target->file, target->event, target->track});
        }
        record = analyze(*parent, *truth, children);
        record.input_file_index = index; record.source_file_id = segment; record.event_in_file = entry;
        record.hepmc_event_number = gen ? gen->event_number() : -1;
        record.weight = event_weight; record.hepmc_weight_available = gen_weight_available;
        // Positive-weight CDFs are probabilities; signed or invalid weights remain in the raw census and are explicitly counted.
        record.weight_valid = std::isfinite(event_weight) && event_weight >= 0;
        if (record.usable && record.length_cm > 0)
        {
          if (std::isfinite(record.pt)) length_pt.Fill(record.pt, std::log10(record.length_cm));
          if (std::isfinite(record.energy)) length_energy.Fill(record.energy, std::log10(record.length_cm));
          if (std::isfinite(record.alignment)) alignment.Fill(record.alignment == 1 ? std::nextafter(1., 0.) : record.alignment);
        }
        eta_tree.Fill(); ++event_eta; event_usable += record.usable; event_unresolved += !record.usable;
        for (std::size_t cut = 0; cut < thresholds_cm.size(); ++cut) event_long_eta[cut] += record.usable && record.length_cm >= thresholds_cm[cut];
        for (const auto& group : group_names(record)) summaries[group].add(record);
      }
      events_with_eta += event_eta > 0; event_tree.Fill();
    }
    tree->ResetBranchAddresses(); delete truth; delete hepmc;
    std::cout << "CheckEtaDecayFlight: file " << index << ", events " << scanned << ", eta so far " << summaries.at("all").total << std::endl;
  }
  requests_matched = matched.size();
  if (resolved_inputs.empty()) throw std::runtime_error("No DST files matched the selected range/requests");
  if (!requests_csv.empty())
  {
    // Requests outside first_file/n_files are allowed; every selected request must be found unless an event limit was requested.
    for (Long64_t index = first_file; index < end_file; ++index) for (const auto& r : requests)
      if (r.file == segment_from_path(inputs[index]) && (r.sample.empty() || r.sample == sample_from_path(inputs[index])) && max_events_per_file == 0 && !matched.count({r.sample, r.file, r.event, r.track}))
        throw std::runtime_error("Requested eta missing from G4-primary truth");
  }
  output.cd(); eta_tree.Write(); event_tree.Write(); length_pt.Write(); length_energy.Write(); alignment.Write();
  int schema_version = 1; bool completed = true;
  std::string scope = requests_csv.empty() ? "inclusive_G4_primary_eta" : "requested_G4_primary_eta_only";
  std::string length_definition = "Euclidean_displacement_decay_vertex_minus_eta_production_vertex_cm_not_accumulated_track_length";
  std::string time_definition = "raw_PHG4VtxPoint_get_t_difference_no_assumed_unit_conversion";
  std::string usable_definition = "confirmed_or_common_mode_inferred_decay_single_vertex_no_missing_vertices_no_negative_delta_t";
  std::string weight_definition = "sample_weight_times_first_HepMC_weight_else_one_nonnegative_finite_weights_for_weighted_CDF";
  std::vector<double> cuts(thresholds_cm.begin(), thresholds_cm.end());
  std::vector<std::string> statuses(status_names.begin(), status_names.end()), modes(mode_names.begin(), mode_names.end());
  std::string metadata_input = input_list, metadata_requests = requests_csv;
  TTree metadata("metadata", "Eta decay flight diagnostic settings");
#define M(member) metadata.Branch(#member, &member)
  M(schema_version); M(completed); M(scope); M(length_definition); M(time_definition); M(usable_definition); M(weight_definition);
  M(metadata_input); M(metadata_requests); M(first_file); M(n_files); M(max_events_per_file); M(embedding_id); M(sample_weight);
  M(resolved_inputs); M(events_read); M(events_with_eta); M(requests_matched); M(cuts); M(statuses); M(modes);
#undef M
  write_summaries(output, base, summaries); draw_report(output, base, summaries);
  output.cd(); metadata.Fill(); metadata.Write(); output.Close();
  const auto& all = summaries.at("all");
  std::cout << "CheckEtaDecayFlight: complete; events/eta/usable/unresolved = " << events_read << '/' << all.total << '/' << all.usable << '/' << all.total - all.usable
            << "; median/p99/max displacement [cm] = " << quantile(all, .5) << '/' << quantile(all, .99) << '/' << quantile(all, 1) << std::endl;
  return 0;
}
}

int CheckEtaDecayFlight(const std::string input_list, const std::string output_base = "PhotonAnalysisTree/output/photon_candidate_selection/qa/eta_decay_flight/check",
                       Long64_t first_file = 0, Long64_t n_files = 1, Long64_t max_events_per_file = 0,
                       const std::string requests_csv = "", int embedding_id = 1, double sample_weight = 1.)
{
  try { return eta_decay_flight::run(input_list, output_base, first_file, n_files, max_events_per_file, requests_csv, embedding_id, sample_weight); }
  catch (const std::exception& error) { std::cerr << "CheckEtaDecayFlight: " << error.what() << std::endl; return 1; }
}
#endif
