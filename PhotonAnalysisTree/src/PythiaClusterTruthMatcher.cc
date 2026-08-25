#include "PythiaClusterTruthMatcher.h"

#include "HepMCPhotonClassifier.h"

#include <calobase/RawCluster.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <g4eval/CaloEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>
#include <g4eval/CaloTruthEval.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>

#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

namespace photon_tree
{
PythiaClusterTruthMatcher::PythiaClusterTruthMatcher() = default;
PythiaClusterTruthMatcher::~PythiaClusterTruthMatcher() = default;

bool PythiaClusterTruthMatcher::begin_event(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    return false;
  }
  if (!initialized_)
  {
    eval_stack_ = std::make_unique<CaloEvalStack>(topNode, "CEMC");
    eval_stack_->set_strict(false);
    eval_stack_->set_verbosity(verbosity_);
    eval_stack_->get_rawcluster_eval()->set_usetowerinfo(true);
    initialized_ = true;
  }
  else
  {
    eval_stack_->next_event(topNode);
  }
  return eval_stack_ && eval_stack_->get_truth_eval()->has_reduced_node_pointers();
}

ClusterTruthMatch PythiaClusterTruthMatcher::match(
    const RawCluster* cluster,
    TowerInfoContainer* towers,
    RawTowerContainer* raw_truth_towers,
    PHG4TruthInfoContainer* truth,
    const PHHepMCGenEventMap* hepmc_event_map,
    bool allocate_split_tower_energy) const
{
  ClusterTruthMatch result;
  if (!cluster || !towers || !truth || !eval_stack_)
  {
    return result;
  }

  std::map<int, std::pair<PHG4Particle*, float>> contribution_by_primary;
  bool saw_shower_provenance = false;
  const auto accumulate_shower = [&](int shower_id, float shower_edep, float allocation) {
    if (!(shower_edep > 0.0F) || !std::isfinite(shower_edep))
    {
      return;
    }
    PHG4Shower* shower = truth->GetShower(shower_id);
    PHG4Particle* primary = shower
                                ? eval_stack_->get_truth_eval()->get_primary_particle(shower)
                                : nullptr;
    if (!primary)
    {
      return;
    }
    auto& accumulated = contribution_by_primary[primary->get_track_id()];
    accumulated.first = primary;
    accumulated.second += allocation * shower_edep;
  };

  const auto tower_range = cluster->get_towers();
  for (auto tower_iter = tower_range.first; tower_iter != tower_range.second; ++tower_iter)
  {
    const unsigned int raw_key = tower_iter->first;
    const int ieta = static_cast<int>(RawTowerDefs::decode_index1(raw_key));
    const int iphi = static_cast<int>(RawTowerDefs::decode_index2(raw_key));
    const unsigned int tower_info_key = TowerInfoDefs::encode_emcal(static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));
    TowerInfo* tower = towers->get_tower_at_key(static_cast<int>(tower_info_key));
    if (!tower)
    {
      return result;
    }

    float allocation = 1.0F;
    if (allocate_split_tower_energy)
    {
      const float tower_energy = tower->get_energy();
      const float cluster_tower_energy = tower_iter->second;
      if (!(tower_energy > 0.0F) || !std::isfinite(tower_energy) ||
          !std::isfinite(cluster_tower_energy) || cluster_tower_energy < 0.0F)
      {
        return result;
      }
      allocation = std::clamp(cluster_tower_energy / tower_energy, 0.0F, 1.0F);
    }

    const TowerInfo::ShowerEdepMap& towerinfo_showers = tower->get_showerEdepMap();
    if (!towerinfo_showers.empty())
    {
      saw_shower_provenance = true;
      for (const auto& [shower_id, shower_edep] : towerinfo_showers)
      {
        accumulate_shower(shower_id, shower_edep, allocation);
      }
      continue;
    }

    RawTower* raw_tower = raw_truth_towers ? raw_truth_towers->getTower(raw_key) : nullptr;
    if (raw_tower)
    {
      const auto raw_shower_range = raw_tower->get_g4showers();
      if (raw_shower_range.first != raw_shower_range.second)
      {
        saw_shower_provenance = true;
      }
      for (auto shower_iter = raw_shower_range.first;
           shower_iter != raw_shower_range.second; ++shower_iter)
      {
        accumulate_shower(shower_iter->first, shower_iter->second, allocation);
      }
    }
  }

  result.valid = cluster->getNTowers() == 0U || saw_shower_provenance;
  if (!result.valid)
  {
    return result;
  }
  for (const auto& [track_id, accumulated] : contribution_by_primary)
  {
    PHG4Particle* primary = accumulated.first;
    const float edep = accumulated.second;
    if (!primary || !(edep > 0.0F) || !std::isfinite(edep))
    {
      continue;
    }
    TruthContributor contributor;
    contributor.g4_track_id = track_id;
    contributor.g4_pdg_id = primary->get_pid();
    contributor.embedding_id = eval_stack_->get_truth_eval()->get_embed(primary);
    contributor.hepmc_barcode = primary->get_barcode();
    contributor.edep = edep;

    const PHHepMCGenEvent* subevent = hepmc_event_map
                                         ? hepmc_event_map->get(contributor.embedding_id)
                                         : nullptr;
    const HepMC::GenEvent* event = subevent ? subevent->getEvent() : nullptr;
    const HepMC::GenParticle* hepmc_particle = event
                                                   ? event->barcode_to_particle(contributor.hepmc_barcode)
                                                   : nullptr;
    if (hepmc_particle)
    {
      contributor.hepmc_valid = 1U;
      contributor.hepmc_pdg_id = hepmc_particle->pdg_id();
      const HepMCPhotonClassification classification = HepMCPhotonClassifier().classify(hepmc_particle);
      if (classification.valid)
      {
        contributor.photon_category = classification.category;
        contributor.photon_source = classification.source;
        contributor.immediate_parent_pdg = classification.immediate_parent_pdg;
        contributor.classification_parent_pdg = classification.classification_parent_pdg;
      }
    }
    result.total_edep += edep;
    result.contributors.push_back(contributor);
  }

  std::sort(result.contributors.begin(), result.contributors.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.edep != rhs.edep)
    {
      return lhs.edep > rhs.edep;
    }
    return lhs.g4_track_id < rhs.g4_track_id;
  });
  if (result.total_edep > 0.0F)
  {
    for (TruthContributor& contributor : result.contributors)
    {
      contributor.fraction = contributor.edep / result.total_edep;
    }
  }
  return result;
}

void ClusterTruthCollection::clear()
{
  *this = ClusterTruthCollection{};
  contributor_offset.push_back(0U);
}

void ClusterTruthCollection::append(const ClusterTruthMatch& match)
{
  if (contributor_offset.empty())
  {
    contributor_offset.push_back(0U);
  }
  const unsigned int cluster_index = static_cast<unsigned int>(valid.size());
  valid.push_back(match.valid ? 1U : 0U);
  total_edep.push_back(match.total_edep);
  n_contributor.push_back(static_cast<unsigned int>(match.contributors.size()));

  const TruthContributor empty;
  const TruthContributor& dominant = match.contributors.empty() ? empty : match.contributors.front();
  dominant_g4_track_id.push_back(dominant.g4_track_id);
  dominant_g4_pdg_id.push_back(dominant.g4_pdg_id);
  dominant_embedding_id.push_back(dominant.embedding_id);
  dominant_hepmc_barcode.push_back(dominant.hepmc_barcode);
  dominant_edep.push_back(dominant.edep);
  dominant_fraction.push_back(dominant.fraction);
  dominant_hepmc_valid.push_back(dominant.hepmc_valid);
  dominant_hepmc_pdg_id.push_back(dominant.hepmc_pdg_id);
  dominant_photon_category.push_back(dominant.photon_category);
  dominant_photon_source.push_back(dominant.photon_source);
  dominant_immediate_parent_pdg.push_back(dominant.immediate_parent_pdg);
  dominant_classification_parent_pdg.push_back(dominant.classification_parent_pdg);

  for (const TruthContributor& contributor : match.contributors)
  {
    contributor_cluster_index.push_back(cluster_index);
    contributor_g4_track_id.push_back(contributor.g4_track_id);
    contributor_g4_pdg_id.push_back(contributor.g4_pdg_id);
    contributor_embedding_id.push_back(contributor.embedding_id);
    contributor_hepmc_barcode.push_back(contributor.hepmc_barcode);
    contributor_edep.push_back(contributor.edep);
    contributor_fraction.push_back(contributor.fraction);
    contributor_hepmc_valid.push_back(contributor.hepmc_valid);
    contributor_hepmc_pdg_id.push_back(contributor.hepmc_pdg_id);
    contributor_photon_category.push_back(contributor.photon_category);
    contributor_photon_source.push_back(contributor.photon_source);
    contributor_immediate_parent_pdg.push_back(contributor.immediate_parent_pdg);
    contributor_classification_parent_pdg.push_back(contributor.classification_parent_pdg);
  }
  contributor_offset.push_back(static_cast<unsigned int>(contributor_g4_track_id.size()));
}

void ClusterTruthCollection::create_branches(TTree* tree, const std::string& prefix)
{
  const auto name = [&prefix](const char* suffix) { return prefix + "_cluster_truth_" + suffix; };
  tree->Branch(name("valid").c_str(), &valid);
  tree->Branch(name("total_edep").c_str(), &total_edep);
  tree->Branch(name("n_contributor").c_str(), &n_contributor);
  tree->Branch(name("dominant_g4_track_id").c_str(), &dominant_g4_track_id);
  tree->Branch(name("dominant_g4_pdg_id").c_str(), &dominant_g4_pdg_id);
  tree->Branch(name("dominant_embedding_id").c_str(), &dominant_embedding_id);
  tree->Branch(name("dominant_hepmc_barcode").c_str(), &dominant_hepmc_barcode);
  tree->Branch(name("dominant_edep").c_str(), &dominant_edep);
  tree->Branch(name("dominant_fraction").c_str(), &dominant_fraction);
  tree->Branch(name("dominant_hepmc_valid").c_str(), &dominant_hepmc_valid);
  tree->Branch(name("dominant_hepmc_pdg_id").c_str(), &dominant_hepmc_pdg_id);
  tree->Branch(name("dominant_photon_category").c_str(), &dominant_photon_category);
  tree->Branch(name("dominant_photon_source").c_str(), &dominant_photon_source);
  tree->Branch(name("dominant_immediate_parent_pdg").c_str(), &dominant_immediate_parent_pdg);
  tree->Branch(name("dominant_classification_parent_pdg").c_str(), &dominant_classification_parent_pdg);
  tree->Branch(name("contributor_offset").c_str(), &contributor_offset);
  tree->Branch(name("contributor_cluster_index").c_str(), &contributor_cluster_index);
  tree->Branch(name("contributor_g4_track_id").c_str(), &contributor_g4_track_id);
  tree->Branch(name("contributor_g4_pdg_id").c_str(), &contributor_g4_pdg_id);
  tree->Branch(name("contributor_embedding_id").c_str(), &contributor_embedding_id);
  tree->Branch(name("contributor_hepmc_barcode").c_str(), &contributor_hepmc_barcode);
  tree->Branch(name("contributor_edep").c_str(), &contributor_edep);
  tree->Branch(name("contributor_fraction").c_str(), &contributor_fraction);
  tree->Branch(name("contributor_hepmc_valid").c_str(), &contributor_hepmc_valid);
  tree->Branch(name("contributor_hepmc_pdg_id").c_str(), &contributor_hepmc_pdg_id);
  tree->Branch(name("contributor_photon_category").c_str(), &contributor_photon_category);
  tree->Branch(name("contributor_photon_source").c_str(), &contributor_photon_source);
  tree->Branch(name("contributor_immediate_parent_pdg").c_str(), &contributor_immediate_parent_pdg);
  tree->Branch(name("contributor_classification_parent_pdg").c_str(), &contributor_classification_parent_pdg);
}
}
