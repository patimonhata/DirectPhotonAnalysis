#include "TopoClusterHCalTree.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerDefs.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>

#include <cmath>
#include <iostream>

TopoClusterHCalTree::TopoClusterHCalTree(const std::string &name)
  : SubsysReco(name)
{
}

TopoClusterHCalTree::~TopoClusterHCalTree()
{
  close_output_file();
}

int TopoClusterHCalTree::Init(PHCompositeNode * /*topNode*/)
{
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cerr << Name() << "::Init - could not create " << output_file_name_ << std::endl;
    close_output_file();
    return Fun4AllReturnCodes::ABORTRUN;
  }

  tree_ = new TTree("topocluster_tree", "All-calorimeter TopoCluster energy components");
  create_branches();

  std::cout << Name() << "::Init - input: " << input_file_name_ << '\n'
            << Name() << "::Init - node: " << topocluster_node_name_ << '\n'
            << Name() << "::Init - output: " << output_file_name_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TopoClusterHCalTree::process_event(PHCompositeNode *topNode)
{
  reset_event();

  auto *clusters = findNode::getClass<RawClusterContainer>(topNode, topocluster_node_name_);
  if (!clusters)
  {
    std::cerr << Name() << "::process_event - missing RawClusterContainer node "
              << topocluster_node_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const std::size_t capacity = clusters->size();
  topocluster_id_.reserve(capacity);
  emcal_energy_.reserve(capacity);
  hcalin_energy_.reserve(capacity);
  hcalout_energy_.reserve(capacity);
  hcal_total_energy_.reserve(capacity);
  other_calo_energy_.reserve(capacity);
  topocluster_energy_.reserve(capacity);
  energy_residual_.reserve(capacity);
  emcal_ntower_.reserve(capacity);
  hcalin_ntower_.reserve(capacity);
  hcalout_ntower_.reserve(capacity);
  other_calo_ntower_.reserve(capacity);

  const RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
  for (auto cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; ++cluster_iter)
  {
    const RawCluster *cluster = cluster_iter->second;
    if (!cluster)
    {
      std::cerr << Name() << "::process_event - null cluster pointer in "
                << topocluster_node_name_ << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    float emcal_energy = 0.0F;
    float hcalin_energy = 0.0F;
    float hcalout_energy = 0.0F;
    float other_calo_energy = 0.0F;
    unsigned int emcal_ntower = 0;
    unsigned int hcalin_ntower = 0;
    unsigned int hcalout_ntower = 0;
    unsigned int other_calo_ntower = 0;

    const RawCluster::TowerConstRange tower_range = cluster->get_towers();
    for (auto tower_iter = tower_range.first; tower_iter != tower_range.second; ++tower_iter)
    {
      const float contribution = tower_iter->second;
      if (!std::isfinite(contribution))
      {
        std::cerr << Name() << "::process_event - non-finite tower contribution in cluster "
                  << cluster->get_id() << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      switch (RawTowerDefs::decode_caloid(tower_iter->first))
      {
      case RawTowerDefs::CEMC:
        emcal_energy += contribution;
        ++emcal_ntower;
        break;
      case RawTowerDefs::HCALIN:
        hcalin_energy += contribution;
        ++hcalin_ntower;
        break;
      case RawTowerDefs::HCALOUT:
        hcalout_energy += contribution;
        ++hcalout_ntower;
        break;
      default:
        other_calo_energy += contribution;
        ++other_calo_ntower;
        break;
      }
    }

    const float hcal_total_energy = hcalin_energy + hcalout_energy;
    const float component_sum = emcal_energy + hcal_total_energy + other_calo_energy;

    topocluster_id_.push_back(static_cast<unsigned int>(cluster->get_id()));
    emcal_energy_.push_back(emcal_energy);
    hcalin_energy_.push_back(hcalin_energy);
    hcalout_energy_.push_back(hcalout_energy);
    hcal_total_energy_.push_back(hcal_total_energy);
    other_calo_energy_.push_back(other_calo_energy);
    topocluster_energy_.push_back(cluster->get_energy());
    energy_residual_.push_back(cluster->get_energy() - component_sum);
    emcal_ntower_.push_back(emcal_ntower);
    hcalin_ntower_.push_back(hcalin_ntower);
    hcalout_ntower_.push_back(hcalout_ntower);
    other_calo_ntower_.push_back(other_calo_ntower);
  }

  n_topocluster_ = static_cast<unsigned int>(topocluster_id_.size());
  event_uid_ = (static_cast<std::uint64_t>(job_index_) << 32U) |
               static_cast<std::uint64_t>(event_);
  tree_->Fill();
  ++event_;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TopoClusterHCalTree::End(PHCompositeNode * /*topNode*/)
{
  if (!output_file_)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  output_file_->cd();
  tree_->Write();
  TNamed("input_file", input_file_name_.c_str()).Write();
  TNamed("sample", sample_name_.c_str()).Write();
  TNamed("topocluster_node", topocluster_node_name_.c_str()).Write();
  close_output_file();

  std::cout << Name() << "::End - wrote " << event_ << " events to "
            << output_file_name_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void TopoClusterHCalTree::create_branches()
{
  tree_->Branch("sample_id", &sample_id_);
  tree_->Branch("job_index", &job_index_);
  tree_->Branch("process_id", &process_id_);
  tree_->Branch("event", &event_);
  tree_->Branch("event_uid", &event_uid_);
  tree_->Branch("n_topocluster", &n_topocluster_);

  tree_->Branch("topocluster_id", &topocluster_id_);
  tree_->Branch("emcal_energy", &emcal_energy_);
  tree_->Branch("hcalin_energy", &hcalin_energy_);
  tree_->Branch("hcalout_energy", &hcalout_energy_);
  tree_->Branch("hcal_total_energy", &hcal_total_energy_);
  tree_->Branch("other_calo_energy", &other_calo_energy_);
  tree_->Branch("topocluster_energy", &topocluster_energy_);
  tree_->Branch("energy_residual", &energy_residual_);

  tree_->Branch("emcal_ntower", &emcal_ntower_);
  tree_->Branch("hcalin_ntower", &hcalin_ntower_);
  tree_->Branch("hcalout_ntower", &hcalout_ntower_);
  tree_->Branch("other_calo_ntower", &other_calo_ntower_);
}

void TopoClusterHCalTree::reset_event()
{
  n_topocluster_ = 0;
  topocluster_id_.clear();
  emcal_energy_.clear();
  hcalin_energy_.clear();
  hcalout_energy_.clear();
  hcal_total_energy_.clear();
  other_calo_energy_.clear();
  topocluster_energy_.clear();
  energy_residual_.clear();
  emcal_ntower_.clear();
  hcalin_ntower_.clear();
  hcalout_ntower_.clear();
  other_calo_ntower_.clear();
}

void TopoClusterHCalTree::close_output_file()
{
  if (output_file_)
  {
    if (output_file_->IsOpen())
    {
      output_file_->Close();
    }
    delete output_file_;
    output_file_ = nullptr;
    tree_ = nullptr;
  }
}
