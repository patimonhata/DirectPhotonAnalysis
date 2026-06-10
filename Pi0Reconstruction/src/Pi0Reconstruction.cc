#include "Pi0Reconstruction.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

Pi0Reconstruction::Pi0Reconstruction(const std::string &name)
  : SubsysReco(name)
{
}

Pi0Reconstruction::~Pi0Reconstruction()
{
  if (output_file_)
  {
    output_file_->Close();
    delete output_file_;
    output_file_ = nullptr;
  }
}

void Pi0Reconstruction::set_output_file_name(const std::string &output_file_name)
{
  output_file_name_ = output_file_name;
}

void Pi0Reconstruction::set_cluster_node_name(const std::string &cluster_node_name)
{
  cluster_node_name_ = cluster_node_name;
}

void Pi0Reconstruction::set_vertex_node_name(const std::string &vertex_node_name)
{
  vertex_node_name_ = vertex_node_name;
}

void Pi0Reconstruction::set_vertex_mode(VertexMode vertex_mode)
{
  vertex_mode_ = vertex_mode;
}

void Pi0Reconstruction::set_abort_on_missing_cluster_node(bool abort_on_missing_cluster_node)
{
  abort_on_missing_cluster_node_ = abort_on_missing_cluster_node;
}

void Pi0Reconstruction::set_abort_on_missing_vertex_node(bool abort_on_missing_vertex_node)
{
  abort_on_missing_vertex_node_ = abort_on_missing_vertex_node;
}

void Pi0Reconstruction::set_min_cluster_energy(double min_cluster_energy)
{
  min_cluster_energy_ = min_cluster_energy;
}

void Pi0Reconstruction::set_mass_histogram_bins(int nbins, double min, double max)
{
  if (nbins <= 0 || min >= max)
  {
    std::cout << "Pi0Reconstruction::set_mass_histogram_bins - invalid binning requested; keeping current settings" << std::endl;
    return;
  }

  mass_histogram_nbins_ = nbins;
  mass_histogram_min_ = min;
  mass_histogram_max_ = max;
}

int Pi0Reconstruction::Init(PHCompositeNode * /*topNode*/)
{
  create_output_directory();

  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cout << "Pi0Reconstruction::Init - failed to open output file: " << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  h_m_gg_ = new TH1D("h_m_gg", "CEMC cluster pair invariant mass;M_{#gamma#gamma} [GeV];Pairs", mass_histogram_nbins_, mass_histogram_min_, mass_histogram_max_);
  h_ncluster_ = new TH1D("h_ncluster", "CEMC clusters per event;N_{cluster};Events", 100, 0.0, 100.0);
  h_cluster_e_ = new TH1D("h_cluster_e", "CEMC cluster energy;E_{cluster} [GeV];Clusters", 200, 0.0, 20.0);
  h_pair_e_asym_ = new TH1D("h_pair_e_asym", "CEMC cluster pair energy asymmetry;(|E_{1}-E_{2}|)/(E_{1}+E_{2});Pairs", 100, -1.0, 1.0);

  std::cout << "Pi0Reconstruction::Init - writing output to " << output_file_name_ << std::endl;
  std::cout << "Pi0Reconstruction::Init - cluster node: " << cluster_node_name_ << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::InitRun(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::process_event(PHCompositeNode *topNode)
{
  if (event_counter_ % 200 == 0)
  {
    std::cout << "Pi0Reconstruction::process_event - event " << event_counter_ << std::endl;
  }
  ++event_counter_;

  RawClusterContainer *cluster_container = findNode::getClass<RawClusterContainer>(topNode, cluster_node_name_);
  if (!cluster_container)
  {
    if (missing_cluster_node_warnings_ < 5)
    {
      std::cout << "Pi0Reconstruction::process_event - missing required cluster node: " << cluster_node_name_ << std::endl;
    }
    ++missing_cluster_node_warnings_;

    return abort_on_missing_cluster_node_ ? Fun4AllReturnCodes::ABORTRUN : Fun4AllReturnCodes::EVENT_OK;
  }

  std::array<double, 3> vertex = {0.0, 0.0, 0.0};
  if (!get_event_vertex(topNode, vertex))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  std::vector<PhotonCandidate> photons;
  photons.reserve(cluster_container->size());

  RawClusterContainer::ConstRange cluster_range = cluster_container->getClusters();
  for (RawClusterContainer::ConstIterator iter = cluster_range.first; iter != cluster_range.second; ++iter)
  {
    const RawCluster *cluster = iter->second;
    if (!cluster)
    {
      continue;
    }

    const double energy = cluster->get_energy();
    if (!std::isfinite(energy) || energy < min_cluster_energy_)
    {
      continue;
    }

    h_cluster_e_->Fill(energy);

    PhotonCandidate candidate;
    if (build_photon_candidate(energy, cluster->get_x(), cluster->get_y(), cluster->get_z(), vertex, candidate))
    {
      photons.push_back(candidate);
    }
  }

  h_ncluster_->Fill(static_cast<double>(photons.size()));

  for (std::size_t i = 0; i < photons.size(); ++i)
  {
    for (std::size_t j = i + 1; j < photons.size(); ++j)
    {
      const PhotonCandidate &first = photons[i];
      const PhotonCandidate &second = photons[j];

      const double total_energy = first.energy + second.energy;
      const double px = first.momentum[0] + second.momentum[0];
      const double py = first.momentum[1] + second.momentum[1];
      const double pz = first.momentum[2] + second.momentum[2];
      const double mass2 = total_energy * total_energy - px * px - py * py - pz * pz;
      const double mass = std::sqrt(std::max(0.0, mass2));

      h_m_gg_->Fill(mass);
      if (total_energy > 0.0)
      {
        h_pair_e_asym_->Fill( abs(first.energy - second.energy) / total_energy);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::ResetEvent(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::Reset(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::EndRun(const int /*runnumber*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::End(PHCompositeNode * /*topNode*/)
{
  if (!output_file_)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  output_file_->cd();
  h_m_gg_->Write();
  h_ncluster_->Write();
  h_cluster_e_->Write();
  h_pair_e_asym_->Write();
  output_file_->Close();
  delete output_file_;
  output_file_ = nullptr;

  std::cout << "Pi0Reconstruction::End - processed " << event_counter_ << " events" << std::endl;
  std::cout << "Pi0Reconstruction::End - wrote " << output_file_name_ << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

bool Pi0Reconstruction::get_event_vertex(PHCompositeNode *topNode, std::array<double, 3> &vertex)
{
  vertex = {0.0, 0.0, 0.0};

  if (vertex_mode_ == VertexMode::Origin)
  {
    return true;
  }

  GlobalVertexMap *vertex_map = findNode::getClass<GlobalVertexMap>(topNode, vertex_node_name_);
  if (!vertex_map || vertex_map->empty())
  {
    if (missing_vertex_node_warnings_ < 5)
    {
      std::cout << "Pi0Reconstruction::get_event_vertex - missing or empty vertex node: " << vertex_node_name_ << "; using origin" << std::endl;
    }
    ++missing_vertex_node_warnings_;
    return !abort_on_missing_vertex_node_;
  }

  const GlobalVertex *global_vertex = vertex_map->begin()->second;
  if (!global_vertex)
  {
    return !abort_on_missing_vertex_node_;
  }

  vertex[0] = global_vertex->get_x();
  vertex[1] = global_vertex->get_y();
  vertex[2] = global_vertex->get_z();

  if (!std::isfinite(vertex[0]) || !std::isfinite(vertex[1]) || !std::isfinite(vertex[2]))
  {
    vertex = {0.0, 0.0, 0.0};
    return !abort_on_missing_vertex_node_;
  }

  return true;
}

bool Pi0Reconstruction::build_photon_candidate(double energy, double x, double y, double z, const std::array<double, 3> &vertex, PhotonCandidate &candidate) const
{
  if (!std::isfinite(energy) || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
  {
    return false;
  }

  const double dx = x - vertex[0];
  const double dy = y - vertex[1];
  const double dz = z - vertex[2];
  const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (distance <= std::numeric_limits<double>::epsilon())
  {
    return false;
  }

  candidate.energy = energy;
  candidate.momentum = {energy * dx / distance, energy * dy / distance, energy * dz / distance};
  return true;
}

void Pi0Reconstruction::create_output_directory() const
{
  const std::string::size_type slash_position = output_file_name_.find_last_of('/');
  if (slash_position == std::string::npos)
  {
    return;
  }

  const std::string directory = output_file_name_.substr(0, slash_position);
  if (!directory.empty())
  {
    gSystem->mkdir(directory.c_str(), true);
  }
}
