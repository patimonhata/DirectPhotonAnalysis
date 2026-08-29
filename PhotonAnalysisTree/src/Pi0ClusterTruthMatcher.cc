#include "Pi0ClusterTruthMatcher.h"

#include <calobase/RawCluster.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>

namespace photon_tree
{
const char* pi0_cluster_truth_match_status_name(Pi0ClusterTruthMatchStatus value)
{
  switch (value)
  {
  case Pi0ClusterTruthMatchStatus::invalid: return "invalid";
  case Pi0ClusterTruthMatchStatus::partial: return "partial";
  case Pi0ClusterTruthMatchStatus::complete: return "complete";
  }
  return "unknown";
}

const char* pi0_cluster_truth_match_failure_name(Pi0ClusterTruthMatchFailure value)
{
  switch (value)
  {
  case Pi0ClusterTruthMatchFailure::none: return "none";
  case Pi0ClusterTruthMatchFailure::missing_input: return "missing_input";
  case Pi0ClusterTruthMatchFailure::missing_tower_info: return "missing_tower_info";
  case Pi0ClusterTruthMatchFailure::invalid_tower_energy: return "invalid_tower_energy";
  case Pi0ClusterTruthMatchFailure::missing_truth_tower: return "missing_truth_tower";
  case Pi0ClusterTruthMatchFailure::missing_g4_cell: return "missing_g4_cell";
  case Pi0ClusterTruthMatchFailure::missing_g4_hit: return "missing_g4_hit";
  case Pi0ClusterTruthMatchFailure::invalid_hit_edep: return "invalid_hit_edep";
  }
  return "unknown";
}

int Pi0ClusterTruthMatcher::direct_gamma_index(int track_id, PHG4TruthInfoContainer* truth, const std::array<int, 2>& direct_gamma_track_ids)
{
  std::set<int> visited;
  while (track_id != 0 && visited.insert(track_id).second)
  {
    for (std::size_t gamma = 0; gamma < direct_gamma_track_ids.size(); ++gamma)
    {
      if (track_id == direct_gamma_track_ids[gamma])
      {
        return static_cast<int>(gamma);
      }
    }
    PHG4Particle* particle = truth ? truth->GetParticle(track_id) : nullptr;
    if (!particle)
    {
      break;
    }
    track_id = particle->get_parent_id();
  }
  return -1;
}

Pi0ClusterTruthMatch Pi0ClusterTruthMatcher::match(
    const RawCluster* cluster,
    TowerInfoContainer* towers,
    RawTowerContainer* raw_truth_towers,
    PHG4CellContainer* cells,
    PHG4HitContainer* hits,
    PHG4TruthInfoContainer* truth,
    const std::array<int, 2>& direct_gamma_track_ids,
    bool allocate_split_tower_energy) const
{
  Pi0ClusterTruthMatch result;
  if (!cluster || !towers || !raw_truth_towers || !cells || !hits || !truth)
  {
    result.failure = Pi0ClusterTruthMatchFailure::missing_input;
    return result;
  }

  const auto tower_range = cluster->get_towers();
  for (auto tower_iter = tower_range.first; tower_iter != tower_range.second; ++tower_iter)
  {
    const unsigned int raw_key = tower_iter->first;
    const int ieta = static_cast<int>(RawTowerDefs::decode_index1(raw_key));
    const int iphi = static_cast<int>(RawTowerDefs::decode_index2(raw_key));
    const float cluster_tower_energy = tower_iter->second;
    ++result.tower_count;
    if (std::isfinite(cluster_tower_energy) && cluster_tower_energy >= 0.0F)
    {
      result.cluster_member_energy += cluster_tower_energy;
    }
    const auto record_failure = [&](Pi0ClusterTruthMatchFailure failure, unsigned long long cell_id = 0, unsigned long long hit_id = 0) {
      if (result.failure == Pi0ClusterTruthMatchFailure::none)
      {
        result.failure = failure;
        result.failure_ieta = ieta;
        result.failure_iphi = iphi;
        result.failure_tower_key = raw_key;
        result.failure_cell_id = cell_id;
        result.failure_hit_id = hit_id;
      }
    };
    const unsigned int tower_info_key = TowerInfoDefs::encode_emcal(static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));
    TowerInfo* tower = towers->get_tower_at_key(static_cast<int>(tower_info_key));
    if (!tower)
    {
      record_failure(Pi0ClusterTruthMatchFailure::missing_tower_info);
      continue;
    }

    float allocation = 1.0F;
    if (allocate_split_tower_energy)
    {
      const float tower_energy = tower->get_energy();
      if (!(tower_energy > 0.0F) || !std::isfinite(tower_energy) ||
          !std::isfinite(cluster_tower_energy) || cluster_tower_energy < 0.0F)
      {
        record_failure(Pi0ClusterTruthMatchFailure::invalid_tower_energy);
        continue;
      }
      allocation = std::clamp(cluster_tower_energy / tower_energy, 0.0F, 1.0F);
    }

    RawTower* raw_tower = raw_truth_towers->getTower(raw_key);
    if (!raw_tower)
    {
      record_failure(Pi0ClusterTruthMatchFailure::missing_truth_tower);
      continue;
    }

    Pi0ClusterTruthMatch tower_result;
    bool tower_valid = true;
    const auto cell_range = raw_tower->get_g4cells();
    for (auto cell_iter = cell_range.first; cell_iter != cell_range.second; ++cell_iter)
    {
      PHG4Cell* cell = cells->findCell(cell_iter->first);
      if (!cell)
      {
        record_failure(Pi0ClusterTruthMatchFailure::missing_g4_cell, cell_iter->first);
        tower_valid = false;
        break;
      }
      const auto hit_range = cell->get_g4hits();
      for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
      {
        PHG4Hit* hit = hits->findHit(hit_iter->first);
        const float hit_edep = hit_iter->second;
        if (!hit)
        {
          record_failure(Pi0ClusterTruthMatchFailure::missing_g4_hit, cell_iter->first, hit_iter->first);
          tower_valid = false;
          break;
        }
        if (!std::isfinite(hit_edep) || hit_edep < 0.0F)
        {
          record_failure(Pi0ClusterTruthMatchFailure::invalid_hit_edep, cell_iter->first, hit_iter->first);
          tower_valid = false;
          break;
        }
        const float allocated_edep = allocation * hit_edep;
        const int gamma = direct_gamma_index(hit->get_trkid(), truth, direct_gamma_track_ids);
        if (gamma >= 0)
        {
          tower_result.gamma_edep[static_cast<std::size_t>(gamma)] += allocated_edep;
        }
        else
        {
          tower_result.other_edep += allocated_edep;
        }
        tower_result.total_edep += allocated_edep;
      }
      if (!tower_valid) break;
    }
    if (!tower_valid)
    {
      continue;
    }
    ++result.matched_tower_count;
    if (std::isfinite(cluster_tower_energy) && cluster_tower_energy >= 0.0F)
    {
      result.matched_cluster_member_energy += cluster_tower_energy;
    }
    result.total_edep += tower_result.total_edep;
    result.gamma_edep[0] += tower_result.gamma_edep[0];
    result.gamma_edep[1] += tower_result.gamma_edep[1];
    result.other_edep += tower_result.other_edep;
  }

  result.cluster_member_energy_coverage = result.cluster_member_energy > 0.0F
      ? result.matched_cluster_member_energy / result.cluster_member_energy
      : (result.matched_tower_count == result.tower_count ? 1.0F : 0.0F);
  if (result.matched_tower_count == result.tower_count)
  {
    result.status = Pi0ClusterTruthMatchStatus::complete;
    result.failure = Pi0ClusterTruthMatchFailure::none;
    result.valid = true;
    result.usable = true;
  }
  else if (result.matched_tower_count > 0U)
  {
    result.status = Pi0ClusterTruthMatchStatus::partial;
  }
  return result;
}

std::vector<Pi0ClusterTruthMatch> Pi0ClusterTruthMatcher::match_many(
    const RawCluster* cluster,
    TowerInfoContainer* towers,
    RawTowerContainer* raw_truth_towers,
    PHG4CellContainer* cells,
    PHG4HitContainer* hits,
    PHG4TruthInfoContainer* truth,
    const std::vector<std::array<int, 2>>& direct_gamma_track_ids,
    bool allocate_split_tower_energy) const
{
  std::vector<Pi0ClusterTruthMatch> results(direct_gamma_track_ids.size());
  if (results.empty()) return results;
  if (!cluster || !towers || !raw_truth_towers || !cells || !hits || !truth)
  {
    for (auto& result : results) result.failure = Pi0ClusterTruthMatchFailure::missing_input;
    return results;
  }

  std::map<int, std::vector<std::pair<std::size_t, std::size_t>>> owners;
  for (std::size_t candidate = 0; candidate < direct_gamma_track_ids.size(); ++candidate)
    for (std::size_t photon = 0; photon < 2U; ++photon)
    {
      const int track_id = direct_gamma_track_ids[candidate][photon];
      if (track_id != 0 && track_id != -999) owners[track_id].push_back({candidate, photon});
    }

  Pi0ClusterTruthMatch common;
  std::vector<std::array<float, 2>> gamma_edep(results.size(), {0.0F, 0.0F});
  const auto tower_range = cluster->get_towers();
  for (auto tower_iter = tower_range.first; tower_iter != tower_range.second; ++tower_iter)
  {
    const unsigned int raw_key = tower_iter->first;
    const int ieta = static_cast<int>(RawTowerDefs::decode_index1(raw_key));
    const int iphi = static_cast<int>(RawTowerDefs::decode_index2(raw_key));
    const float cluster_tower_energy = tower_iter->second;
    ++common.tower_count;
    if (std::isfinite(cluster_tower_energy) && cluster_tower_energy >= 0.0F) common.cluster_member_energy += cluster_tower_energy;
    const auto record_failure = [&](Pi0ClusterTruthMatchFailure failure, unsigned long long cell_id = 0, unsigned long long hit_id = 0) {
      if (common.failure != Pi0ClusterTruthMatchFailure::none) return;
      common.failure = failure;
      common.failure_ieta = ieta;
      common.failure_iphi = iphi;
      common.failure_tower_key = raw_key;
      common.failure_cell_id = cell_id;
      common.failure_hit_id = hit_id;
    };

    const unsigned int tower_info_key = TowerInfoDefs::encode_emcal(static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));
    TowerInfo* tower = towers->get_tower_at_key(static_cast<int>(tower_info_key));
    if (!tower)
    {
      record_failure(Pi0ClusterTruthMatchFailure::missing_tower_info);
      continue;
    }
    float allocation = 1.0F;
    if (allocate_split_tower_energy)
    {
      const float tower_energy = tower->get_energy();
      if (!(tower_energy > 0.0F) || !std::isfinite(tower_energy) || !std::isfinite(cluster_tower_energy) || cluster_tower_energy < 0.0F)
      {
        record_failure(Pi0ClusterTruthMatchFailure::invalid_tower_energy);
        continue;
      }
      allocation = std::clamp(cluster_tower_energy / tower_energy, 0.0F, 1.0F);
    }

    RawTower* raw_tower = raw_truth_towers->getTower(raw_key);
    if (!raw_tower)
    {
      record_failure(Pi0ClusterTruthMatchFailure::missing_truth_tower);
      continue;
    }
    bool tower_valid = true;
    float tower_total_edep = 0.0F;
    std::vector<std::array<float, 2>> tower_gamma_edep(results.size(), {0.0F, 0.0F});
    const auto cell_range = raw_tower->get_g4cells();
    for (auto cell_iter = cell_range.first; cell_iter != cell_range.second; ++cell_iter)
    {
      PHG4Cell* cell = cells->findCell(cell_iter->first);
      if (!cell)
      {
        record_failure(Pi0ClusterTruthMatchFailure::missing_g4_cell, cell_iter->first);
        tower_valid = false;
        break;
      }
      const auto hit_range = cell->get_g4hits();
      for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
      {
        PHG4Hit* hit = hits->findHit(hit_iter->first);
        const float hit_edep = hit_iter->second;
        if (!hit)
        {
          record_failure(Pi0ClusterTruthMatchFailure::missing_g4_hit, cell_iter->first, hit_iter->first);
          tower_valid = false;
          break;
        }
        if (!std::isfinite(hit_edep) || hit_edep < 0.0F)
        {
          record_failure(Pi0ClusterTruthMatchFailure::invalid_hit_edep, cell_iter->first, hit_iter->first);
          tower_valid = false;
          break;
        }
        const float allocated_edep = allocation * hit_edep;
        tower_total_edep += allocated_edep;
        int track_id = hit->get_trkid();
        std::set<int> visited;
        while (track_id != 0 && visited.insert(track_id).second)
        {
          const auto owner = owners.find(track_id);
          if (owner != owners.end())
          {
            for (const auto& [candidate, photon] : owner->second) tower_gamma_edep[candidate][photon] += allocated_edep;
            break;
          }
          PHG4Particle* particle = truth->GetParticle(track_id);
          if (!particle) break;
          track_id = particle->get_parent_id();
        }
      }
      if (!tower_valid) break;
    }
    if (!tower_valid) continue;
    ++common.matched_tower_count;
    if (std::isfinite(cluster_tower_energy) && cluster_tower_energy >= 0.0F) common.matched_cluster_member_energy += cluster_tower_energy;
    common.total_edep += tower_total_edep;
    for (std::size_t candidate = 0; candidate < results.size(); ++candidate)
      for (std::size_t photon = 0; photon < 2U; ++photon)
        gamma_edep[candidate][photon] += tower_gamma_edep[candidate][photon];
  }

  common.cluster_member_energy_coverage = common.cluster_member_energy > 0.0F
      ? common.matched_cluster_member_energy / common.cluster_member_energy
      : (common.matched_tower_count == common.tower_count ? 1.0F : 0.0F);
  if (common.matched_tower_count == common.tower_count)
  {
    common.status = Pi0ClusterTruthMatchStatus::complete;
    common.failure = Pi0ClusterTruthMatchFailure::none;
    common.valid = true;
    common.usable = true;
  }
  else if (common.matched_tower_count > 0U)
  {
    common.status = Pi0ClusterTruthMatchStatus::partial;
  }
  for (std::size_t candidate = 0; candidate < results.size(); ++candidate)
  {
    results[candidate] = common;
    results[candidate].gamma_edep = gamma_edep[candidate];
    results[candidate].other_edep = std::max(0.0F, common.total_edep - gamma_edep[candidate][0] - gamma_edep[candidate][1]);
  }
  return results;
}

void Pi0ClusterTruthCollection::clear()
{
  *this = Pi0ClusterTruthCollection{};
}

void Pi0ClusterTruthCollection::append(const Pi0ClusterTruthMatch& match,
                                       float truth_gamma0_energy,
                                       float truth_gamma1_energy)
{
  match_valid.push_back(match.valid ? 1U : 0U);
  total_edep.push_back(match.total_edep);
  gamma0_edep.push_back(match.gamma_edep[0]);
  gamma1_edep.push_back(match.gamma_edep[1]);
  other_edep.push_back(match.other_edep);
  const float inverse_total = match.total_edep > 0.0F ? 1.0F / match.total_edep : 0.0F;
  gamma0_fraction.push_back(match.gamma_edep[0] * inverse_total);
  gamma1_fraction.push_back(match.gamma_edep[1] * inverse_total);
  other_fraction.push_back(match.other_edep * inverse_total);
  gamma0_recovery.push_back(
      truth_gamma0_energy > 0.0F ? match.gamma_edep[0] / truth_gamma0_energy : 0.0F);
  gamma1_recovery.push_back(
      truth_gamma1_energy > 0.0F ? match.gamma_edep[1] / truth_gamma1_energy : 0.0F);
}

void Pi0ClusterTruthCollection::create_branches(TTree* tree, const std::string& prefix)
{
  const auto name = [&prefix](const char* suffix) {
    return prefix + "_cluster_truth_" + suffix;
  };
  tree->Branch(name("match_valid").c_str(), &match_valid);
  tree->Branch(name("total_edep").c_str(), &total_edep);
  tree->Branch(name("gamma0_edep").c_str(), &gamma0_edep);
  tree->Branch(name("gamma1_edep").c_str(), &gamma1_edep);
  tree->Branch(name("other_edep").c_str(), &other_edep);
  tree->Branch(name("gamma0_fraction").c_str(), &gamma0_fraction);
  tree->Branch(name("gamma1_fraction").c_str(), &gamma1_fraction);
  tree->Branch(name("other_fraction").c_str(), &other_fraction);
  tree->Branch(name("gamma0_recovery").c_str(), &gamma0_recovery);
  tree->Branch(name("gamma1_recovery").c_str(), &gamma1_recovery);
}
}
