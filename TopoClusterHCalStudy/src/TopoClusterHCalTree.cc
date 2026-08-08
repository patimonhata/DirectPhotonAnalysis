#include "TopoClusterHCalTree.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerDefs.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
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

  std::cout << Name() << "::Init - node: " << topocluster_node_name_ << '\n'
            << Name() << "::Init - output: " << output_file_name_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TopoClusterHCalTree::process_event(PHCompositeNode *topNode)
{
  reset_event();

  auto *clusters = findNode::getClass<RawClusterContainer>(topNode, topocluster_node_name_);
  if (!clusters)
  {
    std::cerr << Name() << "::process_event - missing RawClusterContainer node " << topocluster_node_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto *truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truth)
  {
    std::cerr << Name() << "::process_event - missing PHG4TruthInfoContainer node G4TruthInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const auto primary_range = truth->GetPrimaryParticleRange();
  if (primary_range.first == primary_range.second || !primary_range.first->second)
  {
    std::cerr << Name() << "::process_event - missing primary truth particle" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const PHG4Particle *primary = primary_range.first->second;
  truth_pt_ = std::hypot(primary->get_px(), primary->get_py());
  if (!std::isfinite(truth_pt_))
  {
    std::cerr << Name() << "::process_event - non-finite primary truth pT" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const std::size_t capacity = clusters->size();
  emcal_energy_.reserve(capacity);
  hcalin_energy_.reserve(capacity);
  hcalout_energy_.reserve(capacity);
  hcal_total_energy_.reserve(capacity);

  const RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
  for (auto cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; ++cluster_iter)
  {
    const RawCluster *cluster = cluster_iter->second;
    if (!cluster)
    {
      std::cerr << Name() << "::process_event - null cluster pointer in " << topocluster_node_name_ << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    float emcal_energy = 0.0F;
    float hcalin_energy = 0.0F;
    float hcalout_energy = 0.0F;

    const RawCluster::TowerConstRange tower_range = cluster->get_towers();
    for (auto tower_iter = tower_range.first; tower_iter != tower_range.second; ++tower_iter)
    {
      const float contribution = tower_iter->second;
      if (!std::isfinite(contribution))
      {
        std::cerr << Name() << "::process_event - non-finite tower contribution in cluster " << cluster->get_id() << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      switch (RawTowerDefs::decode_caloid(tower_iter->first))
      {
      case RawTowerDefs::CEMC:
        emcal_energy += contribution;
        break;
      case RawTowerDefs::HCALIN:
        hcalin_energy += contribution;
        break;
      case RawTowerDefs::HCALOUT:
        hcalout_energy += contribution;
        break;
      default:
        break;
      }
    }

    emcal_energy_.push_back(emcal_energy);
    hcalin_energy_.push_back(hcalin_energy);
    hcalout_energy_.push_back(hcalout_energy);
    hcal_total_energy_.push_back(hcalin_energy + hcalout_energy);
  }

  n_topocluster_ = static_cast<unsigned int>(emcal_energy_.size());
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
  close_output_file();

  std::cout << Name() << "::End - wrote " << event_ << " events to "
            << output_file_name_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void TopoClusterHCalTree::create_branches()
{
  tree_->Branch("process_id", &process_id_);
  tree_->Branch("event", &event_);
  tree_->Branch("n_topocluster", &n_topocluster_);
  tree_->Branch("truth_pt", &truth_pt_);

  tree_->Branch("emcal_energy", &emcal_energy_);
  tree_->Branch("hcalin_energy", &hcalin_energy_);
  tree_->Branch("hcalout_energy", &hcalout_energy_);
  tree_->Branch("hcal_total_energy", &hcal_total_energy_);
}

void TopoClusterHCalTree::reset_event()
{
  n_topocluster_ = 0;
  truth_pt_ = 0.0F;
  emcal_energy_.clear();
  hcalin_energy_.clear();
  hcalout_energy_.clear();
  hcal_total_energy_.clear();
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
