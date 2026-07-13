#include "TowerClusterEnergyAudit.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>

TowerClusterEnergyAudit::TowerClusterEnergyAudit(const std::string &name)
  : SubsysReco(name)
{
}

TowerClusterEnergyAudit::~TowerClusterEnergyAudit()
{
  if (output_file_)
  {
    output_file_->Close();
    delete output_file_;
    output_file_ = nullptr;
  }
}

void TowerClusterEnergyAudit::set_output_file_name(const std::string &output_file_name)
{
  output_file_name_ = output_file_name;
}

void TowerClusterEnergyAudit::set_tower_node_name(const std::string &tower_node_name)
{
  tower_node_name_ = tower_node_name;
}

void TowerClusterEnergyAudit::set_split_cluster_node_name(const std::string &cluster_node_name)
{
  split_cluster_node_name_ = cluster_node_name;
}

void TowerClusterEnergyAudit::set_no_split_cluster_node_name(const std::string &cluster_node_name)
{
  no_split_cluster_node_name_ = cluster_node_name;
}

void TowerClusterEnergyAudit::set_tower_energy_threshold(double threshold)
{
  tower_energy_threshold_ = threshold;
}

void TowerClusterEnergyAudit::set_cluster_energy_threshold(double threshold)
{
  cluster_energy_threshold_ = threshold;
}

void TowerClusterEnergyAudit::set_abort_on_missing_nodes(bool abort_on_missing_nodes)
{
  abort_on_missing_nodes_ = abort_on_missing_nodes;
}

int TowerClusterEnergyAudit::Init(PHCompositeNode * /*topNode*/)
{
  create_output_directory();
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cout << "TowerClusterEnergyAudit::Init - failed to open " << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  tree_ = new TTree("energy_audit_tree", "Tower and cluster energy audit tree");
  create_tree_branches();

  std::cout << "TowerClusterEnergyAudit::Init - output: " << output_file_name_ << std::endl;
  std::cout << "TowerClusterEnergyAudit::Init - tower node: " << tower_node_name_ << std::endl;
  std::cout << "TowerClusterEnergyAudit::Init - split cluster node: " << split_cluster_node_name_ << std::endl;
  std::cout << "TowerClusterEnergyAudit::Init - no-split cluster node: " << no_split_cluster_node_name_ << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int TowerClusterEnergyAudit::process_event(PHCompositeNode *topNode)
{
  reset_tree_variables();
  EnergySummary tower_summary;
  EnergySummary split_cluster_summary;
  EnergySummary no_split_cluster_summary;

  const bool have_towers = fill_tower_summary(topNode, tower_summary);
  const bool have_split_clusters = fill_cluster_summary(topNode, split_cluster_node_name_, split_cluster_summary);
  const bool have_no_split_clusters = fill_cluster_summary(topNode, no_split_cluster_node_name_, no_split_cluster_summary);

  if ((!have_towers || !have_split_clusters || !have_no_split_clusters) && abort_on_missing_nodes_)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  copy_summary_to_tree(tower_summary,
                       tower_count_,
                       tower_count_above_threshold_,
                       tower_energy_sum_,
                       tower_energy_sum_above_threshold_,
                       tower_max_energy_,
                       tower_energies_);

  copy_summary_to_tree(split_cluster_summary,
                       split_cluster_count_,
                       split_cluster_count_above_threshold_,
                       split_cluster_energy_sum_,
                       split_cluster_energy_sum_above_threshold_,
                       split_cluster_max_energy_,
                       split_cluster_energies_);

  copy_summary_to_tree(no_split_cluster_summary,
                       no_split_cluster_count_,
                       no_split_cluster_count_above_threshold_,
                       no_split_cluster_energy_sum_,
                       no_split_cluster_energy_sum_above_threshold_,
                       no_split_cluster_max_energy_,
                       no_split_cluster_energies_);

  if (tree_)
  {
    tree_->Fill();
  }
  ++event_;

  return Fun4AllReturnCodes::EVENT_OK;
}

int TowerClusterEnergyAudit::End(PHCompositeNode * /*topNode*/)
{
  if (!output_file_)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  output_file_->cd();
  if (tree_)
  {
    tree_->Write();
  }
  output_file_->Close();
  delete output_file_;
  output_file_ = nullptr;

  std::cout << "TowerClusterEnergyAudit::End - wrote " << output_file_name_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool TowerClusterEnergyAudit::fill_tower_summary(PHCompositeNode *topNode, EnergySummary &summary)
{
  TowerInfoContainer *tower_container = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name_);
  if (!tower_container)
  {
    if (missing_node_warnings_ < 10)
    {
      std::cout << "TowerClusterEnergyAudit::fill_tower_summary - missing tower node " << tower_node_name_ << std::endl;
    }
    ++missing_node_warnings_;
    return false;
  }

  const unsigned int ntowers = tower_container->size();
  summary.energies.reserve(ntowers);
  for (unsigned int channel = 0; channel < ntowers; ++channel)
  {
    TowerInfo *tower = tower_container->get_tower_at_channel(channel);
    if (!tower)
    {
      continue;
    }

    const double energy = tower->get_energy();
    if (!std::isfinite(energy))
    {
      continue;
    }

    summary.count++;
    summary.energy_sum += energy;
    summary.max_energy = std::max(summary.max_energy, energy);
    if (energy >= tower_energy_threshold_)
    {
      summary.count_above_threshold++;
      summary.energy_sum_above_threshold += energy;
      summary.energies.push_back(energy);
    }
  }

  return true;
}

bool TowerClusterEnergyAudit::fill_cluster_summary(PHCompositeNode *topNode, const std::string &node_name, EnergySummary &summary)
{
  RawClusterContainer *cluster_container = findNode::getClass<RawClusterContainer>(topNode, node_name);
  if (!cluster_container)
  {
    if (missing_node_warnings_ < 10)
    {
      std::cout << "TowerClusterEnergyAudit::fill_cluster_summary - missing cluster node " << node_name << std::endl;
    }
    ++missing_node_warnings_;
    return false;
  }

  summary.energies.reserve(cluster_container->size());
  RawClusterContainer::ConstRange cluster_range = cluster_container->getClusters();
  for (RawClusterContainer::ConstIterator iter = cluster_range.first; iter != cluster_range.second; ++iter)
  {
    const RawCluster *cluster = iter->second;
    if (!cluster)
    {
      continue;
    }

    const double energy = cluster->get_energy();
    if (!std::isfinite(energy))
    {
      continue;
    }

    summary.count++;
    summary.energy_sum += energy;
    summary.max_energy = std::max(summary.max_energy, energy);
    if (energy >= cluster_energy_threshold_)
    {
      summary.count_above_threshold++;
      summary.energy_sum_above_threshold += energy;
      summary.energies.push_back(energy);
    }
  }

  return true;
}

void TowerClusterEnergyAudit::copy_summary_to_tree(const EnergySummary &summary,
                                                   unsigned int &count,
                                                   unsigned int &count_above_threshold,
                                                   double &energy_sum,
                                                   double &energy_sum_above_threshold,
                                                   double &max_energy,
                                                   std::vector<double> &energies)
{
  count = summary.count;
  count_above_threshold = summary.count_above_threshold;
  energy_sum = summary.energy_sum;
  energy_sum_above_threshold = summary.energy_sum_above_threshold;
  max_energy = summary.max_energy;
  energies = summary.energies;
}

void TowerClusterEnergyAudit::reset_tree_variables()
{
  tower_count_ = 0;
  tower_count_above_threshold_ = 0;
  tower_energy_sum_ = 0.0;
  tower_energy_sum_above_threshold_ = 0.0;
  tower_max_energy_ = 0.0;
  tower_energies_.clear();

  split_cluster_count_ = 0;
  split_cluster_count_above_threshold_ = 0;
  split_cluster_energy_sum_ = 0.0;
  split_cluster_energy_sum_above_threshold_ = 0.0;
  split_cluster_max_energy_ = 0.0;
  split_cluster_energies_.clear();

  no_split_cluster_count_ = 0;
  no_split_cluster_count_above_threshold_ = 0;
  no_split_cluster_energy_sum_ = 0.0;
  no_split_cluster_energy_sum_above_threshold_ = 0.0;
  no_split_cluster_max_energy_ = 0.0;
  no_split_cluster_energies_.clear();
}

void TowerClusterEnergyAudit::create_tree_branches()
{
  if (!tree_)
  {
    return;
  }

  tree_->Branch("event", &event_);
  tree_->Branch("tower_count", &tower_count_);
  tree_->Branch("tower_count_above_threshold", &tower_count_above_threshold_);
  tree_->Branch("tower_energy_sum", &tower_energy_sum_);
  tree_->Branch("tower_energy_sum_above_threshold", &tower_energy_sum_above_threshold_);
  tree_->Branch("tower_max_energy", &tower_max_energy_);
  tree_->Branch("tower_energies", &tower_energies_);

  tree_->Branch("split_cluster_count", &split_cluster_count_);
  tree_->Branch("split_cluster_count_above_threshold", &split_cluster_count_above_threshold_);
  tree_->Branch("split_cluster_energy_sum", &split_cluster_energy_sum_);
  tree_->Branch("split_cluster_energy_sum_above_threshold", &split_cluster_energy_sum_above_threshold_);
  tree_->Branch("split_cluster_max_energy", &split_cluster_max_energy_);
  tree_->Branch("split_cluster_energies", &split_cluster_energies_);

  tree_->Branch("no_split_cluster_count", &no_split_cluster_count_);
  tree_->Branch("no_split_cluster_count_above_threshold", &no_split_cluster_count_above_threshold_);
  tree_->Branch("no_split_cluster_energy_sum", &no_split_cluster_energy_sum_);
  tree_->Branch("no_split_cluster_energy_sum_above_threshold", &no_split_cluster_energy_sum_above_threshold_);
  tree_->Branch("no_split_cluster_max_energy", &no_split_cluster_max_energy_);
  tree_->Branch("no_split_cluster_energies", &no_split_cluster_energies_);
}

void TowerClusterEnergyAudit::create_output_directory() const
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
