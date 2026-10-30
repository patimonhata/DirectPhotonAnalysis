#ifndef PHOTON_CANDIDATE_CANDIDATECOMPOSITION_H
#define PHOTON_CANDIDATE_CANDIDATECOMPOSITION_H

#include "CandidateSelectionPlots.h"

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
  pi0_not_applicable,
  pi0_unknown,
  eta,
  other,
  category_count
};

constexpr double kMajorityThreshold = 0.5;
constexpr int kSignalEmbeddingId = 1;
constexpr std::array<const char*, category_count> kKeys = {
    "denominator", "prompt", "pi0_separated", "pi0_merged", "pi0_single_contaminated", "pi0_missing", "pi0_other", "pi0_not_applicable", "pi0_unknown", "eta", "other"};
constexpr std::array<const char*, category_count> kLabels = {
    "Selected candidates", "Prompt-#gamma", "#pi^{0}: separated", "#pi^{0}: merged", "#pi^{0}: single contaminated",
    "#pi^{0}: missing", "#pi^{0}: other", "#pi^{0}: topology N/A", "#pi^{0}: topology unknown", "#eta", "Other"};
constexpr std::array<int, category_count> kColors = {
    kBlack, kRed + 1, kAzure + 7, kMagenta + 1, kCyan + 2, kGreen + 2, kGray + 1, kViolet + 1, kYellow + 2, kOrange + 7, kGray + 2};

int required_shard_count(const std::string& sample_name)
{
  return sample_name == "jet12" ? 10 : 1;
}

struct Summary
{
  unsigned long long selected = 0;
  unsigned long long prompt = 0;
  unsigned long long pi0 = 0;
  unsigned long long eta = 0;
  unsigned long long other = 0;
  unsigned long long overlap = 0;
  unsigned long long half_boundary = 0;
  unsigned long long invalid_truth = 0;
};

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
      counts[index]->SetDirectory(nullptr);
      weighted[index]->SetDirectory(nullptr);
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
    double sum = 0.0, sumw2 = 0.0;
    for (std::size_t index = 1; index < category_count; ++index) { sum += histograms[index]->GetBinContent(bin); sumw2 += std::pow(histograms[index]->GetBinError(bin), 2); }
    const double scale = std::max({1.0, std::abs(histograms[denominator]->GetBinContent(bin)), std::abs(sum)});
    if (std::abs(histograms[denominator]->GetBinContent(bin) - sum) > 1e-9 * scale || std::abs(sumw2 - std::pow(histograms[denominator]->GetBinError(bin), 2)) > 1e-9 * std::max(1., sumw2)) return false;
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


std::size_t pi0_category(int topology)
{
  if (topology == 1) return pi0_separated;
  if (topology == 2) return pi0_merged;
  if (topology == 4) return pi0_single_contaminated;
  if (topology == 3) return pi0_missing;
  return topology == 0 ? pi0_other : category_count;
}

std::size_t origin_category(bool is_prompt, bool is_pi0, bool is_eta, int topology)
{
  if (static_cast<int>(is_prompt) + static_cast<int>(is_pi0) + static_cast<int>(is_eta) > 1) return other;
  if (is_prompt) return prompt;
  if (is_pi0) return pi0_category(topology);
  return is_eta ? eta : other;
}

void draw_composition_components(const std::vector<TH1D*>& components, const std::vector<std::string>& labels,
                                 const std::vector<int>& colors, const std::string& output, const std::string& detail, const std::string& fraction_axis = "Fraction")
{
  SetsPhenixStyle();
  const bool superdetailed = detail == "superdetailed";
  constexpr int default_canvas_width = 1600;
  const int canvas_width = superdetailed ? 2100 : default_canvas_width;
  constexpr double plot_height = kPlotPadHeight;
  const double plot_width = kCompositionPlotWidth * default_canvas_width / canvas_width;
  TCanvas canvas(("c_" + detail + "_photon_candidate_composition").c_str(), "", canvas_width, kCanvasHeight);
  auto plot_pad = make_plot_pad(detail + "_photon_candidate_composition_pad", false, plot_height, plot_width);
  THStack stack(("stack_" + detail + "_photon_candidate_composition").c_str(), "");
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    components[i]->SetFillColor(colors[i]);
    components[i]->SetLineColor(kBlack);
    stack.Add(components[i]);
  }
  stack.SetMinimum(0.0);
  stack.SetMaximum(1.05);
  auto frame = make_plot_frame(*components.front(), "h_candidate_composition_frame", fraction_axis.c_str(), 0.0, 1.05);
  scale_axes_to_canvas(frame.get(), plot_height);
  if (superdetailed)
  {
    frame->GetXaxis()->SetTitleOffset(1.0);
    frame->GetYaxis()->SetTitleOffset(1.0);
  }
  frame->Draw("AXIS");
  stack.Draw("HIST SAME");
  std::vector<std::unique_ptr<TLegend>> legends;
  canvas.cd();
  if (superdetailed)
  {
    legends.push_back(std::make_unique<TLegend>(0.55, 0.05, 0.99, 0.95));
    legends[0]->SetTextSize(std::min(0.04, 0.82 / std::max(std::size_t{1}, components.size())));
    for (std::size_t i = components.size(); i-- > 0;) add_readable_legend_entry(*legends[0], components[i], labels[i], "f");
  }
  else if (detail == "detailed")
  {
    legends.push_back(std::make_unique<TLegend>(0.72, 0.20, 0.99, 0.80));
    for (std::size_t i = components.size(); i-- > 0;)
    {
      const char* label = i == 3 ? "#pi^{0}: contaminated" : labels[i].c_str();
      legends[0]->AddEntry(components[i], label, "f");
    }
  }
  else
  {
    legends.push_back(std::make_unique<TLegend>(0.72, 0.32, 0.99, 0.68));
    for (std::size_t i = components.size(); i-- > 0;) legends[0]->AddEntry(components[i], labels[i].c_str(), "f");
  }
  for (auto& legend : legends)
  {
    legend->SetFillStyle(0);
    legend->Draw();
  }
  plot_pad->RedrawAxis();
  canvas.cd();
  draw_annotations();
  // ROOT's unembedded Symbol mu can be substituted incorrectly by PDF viewers.
  const bool embed_fonts = output.size() >= 4 && output.compare(output.size() - 4, 4, ".pdf") == 0 &&
      std::any_of(labels.begin(), labels.end(), [](const std::string& label) { return label.find("#mu") != std::string::npos; });
  if (embed_fonts) canvas.Print(output.c_str(), "EmbedFonts");
  else canvas.SaveAs(output.c_str());
}

void draw_stack(std::array<std::unique_ptr<TH1D>, category_count>& fractions, const std::string& output, bool detailed)
{
  SetsPhenixStyle();
  std::vector<TH1D*> components;
  std::vector<std::string> labels;
  std::vector<int> colors;
  const auto add = [&](std::size_t i) { components.push_back(fractions[i].get()); labels.emplace_back(kLabels[i]); colors.push_back(kColors[i]); };
  std::unique_ptr<TH1D> pi0_fraction;
  if (detailed)
  {
    for (std::size_t i = 1; i < category_count; ++i) add(i);
  }
  else
  {
    pi0_fraction.reset(static_cast<TH1D*>(fractions[pi0_separated]->Clone("h_candidate_pi0_fraction")));
    pi0_fraction->Reset("ICES");
    pi0_fraction->SetDirectory(nullptr);
    for (std::size_t i = pi0_separated; i <= pi0_unknown; ++i) pi0_fraction->Add(fractions[i].get());
    add(prompt);
    components.push_back(pi0_fraction.get()); labels.emplace_back("#pi^{0}"); colors.push_back(kColors[pi0_separated]);
    add(eta);
    add(other);
  }
  draw_composition_components(components, labels, colors, output, detailed ? "detailed" : "summary");
}

constexpr std::array<std::size_t, 6> kMissingSpectrumIndices = {6, 7, 9, 10, 8, 11};

// These totals accumulate the same weights in different orders (cluster order vs ET-bin order).
// Use the histogram partition tolerance here; keep same_double strict for metadata and binwise comparisons.
bool same_weighted_integral(double left, double right)
{
  return std::isfinite(left) && std::isfinite(right) &&
      std::abs(left - right) <= 1e-8 * std::max({1.0, std::abs(left), std::abs(right)});
}

bool consistent_pi0_bdt(const candidate_composition::Histograms& composition, const Spectra& topology, const pi0_bdt::Histograms& bdt)
{
  if (!bdt.valid()) return false;
  constexpr std::array<int, 5> indices = {2, 3, 4, 5, 12};
  for (int k = 0; k < 5; ++k)
  {
    const auto c = candidate_composition::pi0_separated + k;
    for (int bin = 0; bin <= composition.counts[c]->GetNbinsX() + 1; ++bin)
      if (!same_double(composition.counts[c]->GetBinContent(bin), topology.counts[indices[k]]->GetBinContent(bin)) ||
          !same_double(composition.weighted[c]->GetBinContent(bin), topology.weighted_pb[indices[k]]->GetBinContent(bin))) return false;
    const auto& summary = bdt.cells[0][k]->summary;
    const double composition_count = composition.counts[c]->Integral(0, composition.counts[c]->GetNbinsX() + 1);
    const double composition_weight = composition.weighted[c]->Integral(0, composition.weighted[c]->GetNbinsX() + 1);
    if (!same_double(summary.GetBinContent(1), composition_count) || !same_weighted_integral(summary.GetBinContent(3), composition_weight))
    {
      std::cerr << std::setprecision(17) << "Pi0 BDT total mismatch for " << pi0_bdt::keys[k]
                << ": selected count=" << summary.GetBinContent(1) << ", composition count=" << composition_count
                << ", selected weight=" << summary.GetBinContent(3) << ", composition weight=" << composition_weight << std::endl;
      return false;
    }
  }
  return true;
}

std::vector<std::unique_ptr<TH1D>> composition_missing_fractions(const Histograms& composition, const Spectra& topology)
{
  // Both populations use the composition pi0 definition; verify their missing histograms agree.
  for (int bin = 0; bin <= composition.counts[pi0_missing]->GetNbinsX() + 1; ++bin)
  {
    if (!same_double(composition.counts[pi0_missing]->GetBinContent(bin), topology.counts[5]->GetBinContent(bin)) ||
        !same_double(composition.weighted[pi0_missing]->GetBinContent(bin), topology.weighted_pb[5]->GetBinContent(bin)) ||
        !same_double(composition.weighted[pi0_missing]->GetBinError(bin), topology.weighted_pb[5]->GetBinError(bin))) return {};
  }
  std::vector<std::unique_ptr<TH1D>> result;
  for (std::size_t index : kMissingSpectrumIndices)
    result.push_back(fraction_histogram(*topology.weighted_pb[index], *composition.weighted[denominator], std::string("h_candidate_pi0_") + ::kKeys[index] + "_fraction"));
  return result;
}

void draw_superdetailed_stack(std::array<std::unique_ptr<TH1D>, category_count>& fractions,
                              const std::vector<std::unique_ptr<TH1D>>& missing, const std::string& output)
{
  SetsPhenixStyle();
  std::vector<TH1D*> components;
  std::vector<std::string> labels;
  std::vector<int> colors;
  for (std::size_t i = 1; i < category_count; ++i)
  {
    if (i == pi0_missing)
    {
      for (std::size_t j = 0; j < missing.size(); ++j)
      {
        const std::size_t index = kMissingSpectrumIndices[j];
        components.push_back(missing[j].get());
        labels.push_back("#pi^{0}: " + ::kLabels[index]);
        colors.push_back(::kColors[index]);
      }
    }
    else
    {
      components.push_back(fractions[i].get()); labels.emplace_back(kLabels[i]); colors.push_back(kColors[i]);
    }
  }
  draw_composition_components(components, labels, colors, output, "superdetailed");
}

}

#endif
