#include "ReducePythiaRegionAAnchorTopology.C"

#include <numeric>

namespace candidate_composition
{
enum Category : std::size_t
{
  denominator,
  prompt,
  pi0_separated,
  pi0_merged,
  pi0_single_contaminated,
  pi0_missing,
  pi0_other,
  eta,
  other,
  category_count
};

constexpr double kMajorityThreshold = 0.5;
constexpr int kSignalEmbeddingId = 1;
constexpr std::array<const char*, category_count> kKeys = {
    "denominator", "prompt", "pi0_separated", "pi0_merged", "pi0_single_contaminated", "pi0_missing", "pi0_other", "eta", "other"};
constexpr std::array<const char*, category_count> kLabels = {
    "Selected candidates", "Prompt #gamma", "#pi^{0}: separated", "#pi^{0}: merged", "#pi^{0}: single contaminated",
    "#pi^{0}: missing", "#pi^{0}: other", "#eta", "Other"};
constexpr std::array<int, category_count> kColors = {
    kBlack, kRed + 1, kAzure + 7, kMagenta + 1, kCyan + 2, kGreen + 2, kGray + 1, kOrange + 7, kGray + 2};

struct Histograms
{
  std::array<std::unique_ptr<TH1D>, category_count> counts;
  std::array<std::unique_ptr<TH1D>, category_count> weighted;

  Histograms(int n_bins, double et_max)
  {
    for (std::size_t index = 0; index < category_count; ++index)
    {
      counts[index] = std::make_unique<TH1D>((std::string("h_candidate_") + kKeys[index] + "_et_count").c_str(), "", n_bins, 0.0, et_max);
      weighted[index] = std::make_unique<TH1D>((std::string("h_candidate_") + kKeys[index] + "_et_pb").c_str(), "", n_bins, 0.0, et_max);
      counts[index]->Sumw2();
      weighted[index]->Sumw2();
    }
  }

  void fill(std::size_t index, double et, double weight)
  {
    counts[index]->Fill(et);
    weighted[index]->Fill(et, weight);
  }
};

bool valid_partition(const std::array<std::unique_ptr<TH1D>, category_count>& histograms)
{
  for (int bin = 0; bin <= histograms.front()->GetNbinsX() + 1; ++bin)
  {
    double sum = 0.0;
    for (std::size_t index = 1; index < category_count; ++index) sum += histograms[index]->GetBinContent(bin);
    const double scale = std::max({1.0, std::abs(histograms[denominator]->GetBinContent(bin)), std::abs(sum)});
    if (std::abs(histograms[denominator]->GetBinContent(bin) - sum) > 1e-9 * scale) return false;
  }
  return true;
}

std::unique_ptr<TH1D> fraction_histogram(const TH1D& numerator, const TH1D& total, const std::string& name)
{
  auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(numerator.Clone(name.c_str())));
  result->Reset("ICES");
  result->SetDirectory(nullptr);
  for (int bin = 0; bin <= result->GetNbinsX() + 1; ++bin)
  {
    const double den = total.GetBinContent(bin);
    if (den == 0.0) continue;
    const double num = numerator.GetBinContent(bin);
    const double value = num / den;
    const double num_w2 = std::pow(numerator.GetBinError(bin), 2);
    const double den_w2 = std::pow(total.GetBinError(bin), 2);
    const double complement_w2 = std::max(0.0, den_w2 - num_w2);
    const double variance = (std::pow(1.0 - value, 2) * num_w2 + std::pow(value, 2) * complement_w2) / std::pow(den, 2);
    result->SetBinContent(bin, value);
    result->SetBinError(bin, std::sqrt(std::max(0.0, variance)));
  }
  return result;
}

bool valid_contributors(unsigned int ncluster, const std::vector<unsigned int>& offset, const std::vector<int>& g4_pdg,
                        const std::vector<int>& embedding, const std::vector<float>& fraction,
                        const std::vector<int>& source, const std::vector<int>& parent)
{
  const std::size_t size = fraction.size();
  return offset.size() == static_cast<std::size_t>(ncluster) + 1U && !offset.empty() && offset.front() == 0U && offset.back() == size &&
      g4_pdg.size() == size && embedding.size() == size && source.size() == size && parent.size() == size &&
      std::adjacent_find(offset.begin(), offset.end(), std::greater<unsigned int>()) == offset.end();
}

double eta_fraction(std::size_t cluster, const std::vector<unsigned int>& offset, const std::vector<int>& g4_pdg,
                    const std::vector<int>& embedding, const std::vector<float>& fraction,
                    const std::vector<int>& source, const std::vector<int>& parent)
{
  double result = 0.0;
  for (std::size_t contributor = offset[cluster]; contributor < offset[cluster + 1U]; ++contributor)
  {
    if (embedding[contributor] != kSignalEmbeddingId) continue;
    const bool g4_eta = std::abs(g4_pdg[contributor]) == 221;
    const bool generator_eta_photon = g4_pdg[contributor] == 22 && (parent[contributor] == 221 || source[contributor] == 3);
    if (g4_eta || generator_eta_photon) result += fraction[contributor];
  }
  return result;
}

std::size_t pi0_category(int topology)
{
  if (topology == 1) return pi0_separated;
  if (topology == 2) return pi0_merged;
  if (topology == 4) return pi0_single_contaminated;
  if (topology == 3) return pi0_missing;
  return topology == 0 ? pi0_other : category_count;
}

void draw_stack(std::array<std::unique_ptr<TH1D>, category_count>& fractions, const std::string& output,
                const std::string& family, const std::string& selection, double min_cluster_energy, bool detailed)
{
  SetsPhenixStyle();
  const std::string detail = detailed ? "detailed" : "summary";
  TCanvas canvas(("c_" + detail + "_photon_candidate_composition").c_str(), "", kCanvasWidth, kCanvasHeight);
  auto plot_pad = make_plot_pad(detail + "_photon_candidate_composition_pad");
  THStack stack(("stack_" + detail + "_photon_candidate_composition").c_str(), "");
  std::unique_ptr<TH1D> pi0_fraction;
  if (detailed)
  {
    for (std::size_t index = 1; index < category_count; ++index)
    {
      fractions[index]->SetFillColor(kColors[index]);
      fractions[index]->SetLineColor(kBlack);
      stack.Add(fractions[index].get());
    }
  }
  else
  {
    pi0_fraction.reset(static_cast<TH1D*>(fractions[pi0_separated]->Clone("h_candidate_pi0_fraction")));
    pi0_fraction->Reset("ICES");
    pi0_fraction->SetDirectory(nullptr);
    for (std::size_t index = pi0_separated; index <= pi0_other; ++index) pi0_fraction->Add(fractions[index].get());
    pi0_fraction->SetFillColor(kColors[pi0_separated]);
    pi0_fraction->SetLineColor(kBlack);
    for (std::size_t index : {prompt, eta, other})
    {
      fractions[index]->SetFillColor(kColors[index]);
      fractions[index]->SetLineColor(kBlack);
    }
    stack.Add(fractions[prompt].get());
    stack.Add(pi0_fraction.get());
    stack.Add(fractions[eta].get());
    stack.Add(fractions[other].get());
  }
  stack.SetMinimum(0.0);
  stack.SetMaximum(1.05);
  stack.Draw("HIST");
  stack.GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
  stack.GetYaxis()->SetTitle("Fraction of selected candidates");
  for (TAxis* axis : {stack.GetXaxis(), stack.GetYaxis()})
  {
    axis->SetLabelSize(0.05);
    axis->SetTitleSize(0.05);
    axis->CenterTitle();
  }
  stack.GetXaxis()->SetTitleOffset(1.4);
  stack.GetYaxis()->SetTitleOffset(1.15);
  canvas.cd();
  TLegend legend(detailed ? 0.48 : 0.56, detailed ? 0.60 : 0.73, 0.95, 0.97);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(detailed ? 0.018 : 0.024);
  if (detailed)
  {
    for (std::size_t index = 1; index < category_count; ++index) legend.AddEntry(fractions[index].get(), kLabels[index], "f");
  }
  else
  {
    legend.AddEntry(fractions[prompt].get(), kLabels[prompt], "f");
    legend.AddEntry(pi0_fraction.get(), "#pi^{0}", "f");
    legend.AddEntry(fractions[eta].get(), kLabels[eta], "f");
    legend.AddEntry(fractions[other].get(), kLabels[other], "f");
  }
  legend.Draw();
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(0.026);
  label.DrawLatex(0.06, 0.96, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.06, 0.91, (family == "jet" ? "Pythia8 p+p Jet samples" : "Pythia8 p+p PhotonJet samples"));
  label.DrawLatex(0.06, 0.86, selection == "region_a" ? "Region A: isolated and tight" : "Region A after #pi^{0}/#eta tag veto");
  std::ostringstream threshold;
  threshold << "E_{cluster} > " << min_cluster_energy << " GeV; truth contribution > 50%";
  label.DrawLatex(0.06, 0.81, threshold.str().c_str());
  plot_pad->cd();
  plot_pad->RedrawAxis();
  canvas.cd();
  canvas.SaveAs(output.c_str());
}
}

int ReducePythiaPhotonCandidateComposition(
    const std::string family = "jet",
    const std::string map_root = "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p1",
    const std::string output_base = "",
    const std::string selection = "region_a",
    const bool require_complete = true,
    const int n_bins = 200,
    const double et_max = 40.0)
{
  using namespace candidate_composition;
  const std::vector<SampleDefinition> samples = sample_definitions(family);
  if (samples.empty() || map_root.empty() || (selection != "region_a" && selection != "final_photon") ||
      n_bins <= 0 || !std::isfinite(et_max) || et_max <= 0.0) return 1;
  std::string normalized_map_root = map_root;
  while (normalized_map_root.size() > 1U && normalized_map_root.back() == '/') normalized_map_root.pop_back();
  const std::string configuration = normalized_map_root.substr(normalized_map_root.find_last_of('/') + 1U);
  const std::string resolved_output_base = output_base.empty()
      ? "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/plots/photon_candidate_selection/candidate_composition/" +
            configuration + "/" + selection + "/" + family + "/photon_candidate_composition"
      : output_base;
  if (!make_output_directory(resolved_output_base)) return 2;

  Histograms histograms(n_bins, et_max);
  std::vector<std::string> sample_names;
  std::vector<unsigned long long> sample_map_counts;
  std::vector<double> sample_sum_generator_weights;
  std::string analysis_release, model_sha256;
  double min_cluster_energy = -1.0;
  double partner_diagnostic_min_cluster_energy = -1.0;
  double meson_partner_min_energy = -1.0;
  int pi0_topology_algorithm_version = -1;
  unsigned long long selected_count = 0, prompt_count = 0, pi0_count = 0, eta_count = 0, other_count = 0;
  unsigned long long overlap_count = 0, half_boundary_count = 0, invalid_truth_count = 0;

  for (const SampleDefinition& sample : samples)
  {
    std::vector<MapMetadata> maps;
    if (!collect_maps(map_root + "/" + sample.name + "/map_*.root", sample, require_complete, maps)) return 3;
    if (maps.empty()) continue;
    if (analysis_release.empty())
    {
      analysis_release = maps.front().analysis_release;
      model_sha256 = maps.front().model_sha256;
      min_cluster_energy = maps.front().min_cluster_energy;
      partner_diagnostic_min_cluster_energy = maps.front().partner_diagnostic_min_cluster_energy;
      meson_partner_min_energy = maps.front().meson_partner_min_energy;
      pi0_topology_algorithm_version = maps.front().pi0_topology_algorithm_version;
    }
    else if (analysis_release != maps.front().analysis_release || model_sha256 != maps.front().model_sha256 ||
        !same_double(min_cluster_energy, maps.front().min_cluster_energy) ||
        !same_double(partner_diagnostic_min_cluster_energy, maps.front().partner_diagnostic_min_cluster_energy) ||
        !same_double(meson_partner_min_energy, maps.front().meson_partner_min_energy) ||
        pi0_topology_algorithm_version != maps.front().pi0_topology_algorithm_version) return 3;
    double sample_sumw = 0.0;
    for (const auto& map : maps) sample_sumw += map.sum_generator_weight_processed;
    if (!std::isfinite(sample_sumw) || sample_sumw <= 0.0) return 4;
    sample_names.push_back(sample.name);
    sample_map_counts.push_back(maps.size());
    sample_sum_generator_weights.push_back(sample_sumw);

    for (const MapMetadata& map : maps)
    {
      TFile file(map.path.c_str(), "READ");
      auto* tree = file.Get<TTree>("event_tree");
      if (file.IsZombie() || !tree) return 5;
      unsigned int ncluster = 0;
      unsigned char weight_valid = 0, stitch_valid = 0, stitch_pass = 0;
      double weight_numerator = 0.0;
      std::vector<double>* et = nullptr;
      std::vector<unsigned char>* selected = nullptr;
      std::vector<unsigned char>* truth_valid = nullptr;
      std::vector<unsigned char>* prompt_flag = nullptr;
      std::vector<float>* dominant_fraction = nullptr;
      std::vector<unsigned char>* pi0_valid = nullptr;
      std::vector<float>* pi0_fraction = nullptr;
      std::vector<int>* topology = nullptr;
      std::vector<unsigned int>* offset = nullptr;
      std::vector<int>* g4_pdg = nullptr;
      std::vector<int>* embedding = nullptr;
      std::vector<float>* fraction = nullptr;
      std::vector<int>* source = nullptr;
      std::vector<int>* parent = nullptr;
      tree->SetBranchStatus("*", false);
      tree->SetCacheSize(kEventTreeCacheSize);
      bool ok = bind_active(tree, "split_ncluster", &ncluster) &&
          bind_active(tree, "event_weight_valid", &weight_valid) &&
          bind_active(tree, "sample_stitching_valid", &stitch_valid) &&
          bind_active(tree, "sample_stitching_pass", &stitch_pass) &&
          bind_active(tree, "weight_numerator_pb", &weight_numerator) &&
          bind_active(tree, "split_cluster_et", &et) &&
          bind_active(tree, selection == "region_a" ? "split_cluster_pass_region_a" : "split_cluster_pass_final_photon", &selected) &&
          bind_active(tree, "split_cluster_truth_valid", &truth_valid) &&
          bind_active(tree, "split_cluster_truth_prompt_cluster", &prompt_flag) &&
          bind_active(tree, "split_cluster_truth_dominant_fraction", &dominant_fraction) &&
          bind_active(tree, "split_cluster_pi0_anchor_valid", &pi0_valid) &&
          bind_active(tree, "split_cluster_pi0_anchor_main_fraction", &pi0_fraction) &&
          bind_active(tree, "split_cluster_pi0_anchor_topology", &topology) &&
          bind_active(tree, "split_cluster_truth_contributor_offset", &offset) &&
          bind_active(tree, "split_cluster_truth_contributor_g4_pdg_id", &g4_pdg) &&
          bind_active(tree, "split_cluster_truth_contributor_embedding_id", &embedding) &&
          bind_active(tree, "split_cluster_truth_contributor_fraction", &fraction) &&
          bind_active(tree, "split_cluster_truth_contributor_photon_source", &source) &&
          bind_active(tree, "split_cluster_truth_contributor_classification_parent_pdg", &parent);
      if (!ok) return 5;
      tree->StopCacheLearningPhase();

      for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
      {
        tree->GetEntry(entry);
        if (!et || !selected || !truth_valid || !prompt_flag || !dominant_fraction || !pi0_valid || !pi0_fraction || !topology ||
            !offset || !g4_pdg || !embedding || !fraction || !source || !parent || et->size() != ncluster || selected->size() != ncluster ||
            truth_valid->size() != ncluster || prompt_flag->size() != ncluster || dominant_fraction->size() != ncluster ||
            pi0_valid->size() != ncluster || pi0_fraction->size() != ncluster || topology->size() != ncluster ||
            !valid_contributors(ncluster, *offset, *g4_pdg, *embedding, *fraction, *source, *parent)) return 6;
        if (!weight_valid || !stitch_valid || !stitch_pass) continue;
        const double weight = weight_numerator / sample_sumw;
        if (!std::isfinite(weight)) return 6;
        for (std::size_t cluster = 0; cluster < ncluster; ++cluster)
        {
          if (!(*selected)[cluster]) continue;
          const double cluster_et = (*et)[cluster];
          const double eta_contribution = eta_fraction(cluster, *offset, *g4_pdg, *embedding, *fraction, *source, *parent);
          if (!std::isfinite(cluster_et) || !std::isfinite(eta_contribution) || eta_contribution < -1e-6 || eta_contribution > 1.0 + 1e-5) return 6;
          ++selected_count;
          histograms.fill(denominator, cluster_et, weight);
          if (!(*truth_valid)[cluster]) ++invalid_truth_count;
          const bool is_prompt = (*prompt_flag)[cluster] && (*dominant_fraction)[cluster] > kMajorityThreshold;
          const bool is_pi0 = (*pi0_valid)[cluster] && (*pi0_fraction)[cluster] > kMajorityThreshold;
          const bool is_eta = eta_contribution > kMajorityThreshold;
          const int majority_count = static_cast<int>(is_prompt) + static_cast<int>(is_pi0) + static_cast<int>(is_eta);
          if (((*prompt_flag)[cluster] && (*dominant_fraction)[cluster] == kMajorityThreshold) ||
              ((*pi0_valid)[cluster] && (*pi0_fraction)[cluster] == kMajorityThreshold) || eta_contribution == kMajorityThreshold) ++half_boundary_count;
          std::size_t category = other;
          if (majority_count > 1) ++overlap_count;
          else if (is_prompt) category = prompt;
          else if (is_pi0) category = pi0_category((*topology)[cluster]);
          else if (is_eta) category = eta;
          if (category == category_count) return 6;
          histograms.fill(category, cluster_et, weight);
          if (category == prompt) ++prompt_count;
          else if (category >= pi0_separated && category <= pi0_other) ++pi0_count;
          else if (category == eta) ++eta_count;
          else ++other_count;
        }
      }
    }
  }

  if (sample_names.empty() || selected_count != prompt_count + pi0_count + eta_count + other_count ||
      !valid_partition(histograms.counts) || !valid_partition(histograms.weighted)) return 7;
  std::array<std::unique_ptr<TH1D>, category_count> fractions;
  for (std::size_t index = 1; index < category_count; ++index)
  {
    const std::string name = index == prompt ? "h_photon_candidate_purity"
                                             : std::string("h_candidate_") + candidate_composition::kKeys[index] + "_fraction";
    fractions[index] = fraction_histogram(*histograms.weighted[index], *histograms.weighted[denominator], name);
  }
  draw_stack(fractions, resolved_output_base + "_category_fraction_stack.pdf", family, selection, min_cluster_energy, false);
  draw_stack(fractions, resolved_output_base + "_category_fraction_stack_detailed.pdf", family, selection, min_cluster_energy, true);

  TFile output((resolved_output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) return 8;
  for (std::size_t index = 0; index < category_count; ++index)
  {
    histograms.counts[index]->Write();
    histograms.weighted[index]->Write();
    if (index > 0) fractions[index]->Write();
  }
  int schema_version = 1, source_schema_version = 4, signal_embedding_id = kSignalEmbeddingId;
  double majority_threshold = kMajorityThreshold;
  std::string majority_comparison = "strictly_greater_than";
  std::string eta_definition = "sum_signal_embedding_g4_eta_or_generator_eta_decay_photon_contributor_fraction";
  std::string other_definition = "no_unique_prompt_pi0_eta_strict_majority_including_exact_0p5";
  std::string weight_definition = "cross_section_pb_times_generator_weight_divided_by_full_sample_sum_generator_weight_processed";
  std::string metadata_family = family;
  std::string metadata_selection = selection;
  bool metadata_require_complete = require_complete;
  int metadata_n_bins = n_bins;
  double metadata_et_max = et_max;
  TTree metadata("metadata", "Photon-candidate composition reduce metadata");
  metadata.Branch("schema_version", &schema_version);
  metadata.Branch("source_map_schema_version", &source_schema_version);
  metadata.Branch("family", &metadata_family);
  metadata.Branch("selection", &metadata_selection);
  metadata.Branch("map_root", &normalized_map_root);
  metadata.Branch("require_complete", &metadata_require_complete);
  metadata.Branch("n_bins", &metadata_n_bins);
  metadata.Branch("et_max", &metadata_et_max);
  metadata.Branch("min_cluster_energy", &min_cluster_energy);
  metadata.Branch("partner_diagnostic_min_cluster_energy", &partner_diagnostic_min_cluster_energy);
  metadata.Branch("meson_partner_min_energy", &meson_partner_min_energy);
  metadata.Branch("pi0_topology_algorithm_version", &pi0_topology_algorithm_version);
  metadata.Branch("signal_embedding_id", &signal_embedding_id);
  metadata.Branch("majority_threshold", &majority_threshold);
  metadata.Branch("majority_comparison", &majority_comparison);
  metadata.Branch("eta_definition", &eta_definition);
  metadata.Branch("other_definition", &other_definition);
  metadata.Branch("weight_definition", &weight_definition);
  metadata.Branch("analysis_release", &analysis_release);
  metadata.Branch("model_sha256", &model_sha256);
  metadata.Branch("sample_names", &sample_names);
  metadata.Branch("sample_map_counts", &sample_map_counts);
  metadata.Branch("sample_sum_generator_weights", &sample_sum_generator_weights);
  metadata.Branch("selected_cluster_count", &selected_count);
  metadata.Branch("prompt_cluster_count", &prompt_count);
  metadata.Branch("pi0_cluster_count", &pi0_count);
  metadata.Branch("eta_cluster_count", &eta_count);
  metadata.Branch("other_cluster_count", &other_count);
  metadata.Branch("overlap_cluster_count", &overlap_count);
  metadata.Branch("half_boundary_cluster_count", &half_boundary_count);
  metadata.Branch("invalid_truth_cluster_count", &invalid_truth_count);
  metadata.Fill();
  metadata.Write();
  output.Close();
  std::cout << "ReducePythiaPhotonCandidateComposition - family/selection/selected/prompt/pi0/eta/other/overlap/half/output = "
            << family << "/" << selection << "/" << selected_count << "/" << prompt_count << "/" << pi0_count << "/" << eta_count << "/"
            << other_count << "/" << overlap_count << "/" << half_boundary_count << "/" << resolved_output_base << std::endl;
  return overlap_count == 0ULL ? 0 : 9;
}
