#include "ShowerShapeCalculator.h"

#include <calobase/RawCluster.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <cmath>
#include <set>

namespace
{
int wrap_phi_index(int iphi)
{
  iphi %= ShowerShapeCalculator::kCemcPhiBins;
  if (iphi < 0)
  {
    iphi += ShowerShapeCalculator::kCemcPhiBins;
  }
  return iphi;
}

float wrap_phi_coordinate(float iphi)
{
  iphi = std::fmod(iphi, static_cast<float>(ShowerShapeCalculator::kCemcPhiBins));
  if (iphi < 0.0F)
  {
    iphi += static_cast<float>(ShowerShapeCalculator::kCemcPhiBins);
  }
  return iphi;
}

float wrapped_delta_phi(float iphi, float reference_iphi)
{
  float delta = iphi - reference_iphi;
  const float half_range = 0.5F * static_cast<float>(ShowerShapeCalculator::kCemcPhiBins);
  if (delta > half_range)
  {
    delta -= static_cast<float>(ShowerShapeCalculator::kCemcPhiBins);
  }
  else if (delta < -half_range)
  {
    delta += static_cast<float>(ShowerShapeCalculator::kCemcPhiBins);
  }
  return delta;
}

std::size_t patch_index(int eta_offset, int phi_offset)
{
  return static_cast<std::size_t>((eta_offset + 3) * ShowerShapeCalculator::kPatchSide + (phi_offset + 3));
}
}  // namespace

ShowerShapeCalculator::ShowerShapeCalculator(const Config &config)
  : config_(config)
{
}

ShowerShapeCalculator::Result ShowerShapeCalculator::calculate(
    const RawCluster &cluster,
    TowerInfoContainer &cemc_towers) const
{
  Result result;
  result.full_containment = true;
  result.tower_data_complete = true;

  const RawCluster::TowerMap &tower_map = cluster.get_towermap();
  if (tower_map.empty())
  {
    result.full_containment = false;
    result.tower_data_complete = false;
    return result;
  }

  RawTowerDefs::keytype lead_key = 0;
  float lead_energy = 0.0F;
  for (const auto &[tower_key, tower_energy] : tower_map)
  {
    if (std::isfinite(tower_energy) && tower_energy > lead_energy)
    {
      lead_key = tower_key;
      lead_energy = tower_energy;
    }
  }
  if (!(lead_energy > 0.0F))
  {
    return result;
  }

  const int lead_ieta = RawTowerDefs::decode_index1(lead_key);
  const int lead_iphi = RawTowerDefs::decode_index2(lead_key);
  float weighted_deta = 0.0F;
  float weighted_dphi = 0.0F;
  std::set<unsigned int> owned_tower_keys;

  for (const auto &[tower_key, tower_energy] : tower_map)
  {
    const int ieta = RawTowerDefs::decode_index1(tower_key);
    const int iphi = wrap_phi_index(RawTowerDefs::decode_index2(tower_key));
    if (ieta >= 0 && ieta < kCemcEtaBins)
    {
      owned_tower_keys.insert(TowerInfoDefs::encode_emcal(static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi)));
    }

    if (!std::isfinite(tower_energy) || tower_energy <= config_.min_tower_energy)
    {
      continue;
    }

    result.cluster_energy_above_threshold += tower_energy;
    weighted_deta += tower_energy * static_cast<float>(ieta - lead_ieta);
    weighted_dphi += tower_energy * wrapped_delta_phi(static_cast<float>(iphi), static_cast<float>(lead_iphi));
  }

  if (!(result.cluster_energy_above_threshold > 0.0F))
  {
    return result;
  }

  result.cog_ieta = static_cast<float>(lead_ieta) + weighted_deta / result.cluster_energy_above_threshold;
  result.cog_iphi = wrap_phi_coordinate(
      static_cast<float>(lead_iphi) + weighted_dphi / result.cluster_energy_above_threshold);

  const int center_ieta = static_cast<int>(std::floor(result.cog_ieta + 0.5F));
  const int center_iphi_unwrapped = static_cast<int>(std::floor(result.cog_iphi + 0.5F));
  const int center_iphi = wrap_phi_index(center_iphi_unwrapped);
  const float eta_shift = result.cog_ieta - static_cast<float>(center_ieta);
  const float phi_shift = wrapped_delta_phi(result.cog_iphi, static_cast<float>(center_iphi));

  const int eta_direction = eta_shift < 0.0F ? -1 : 1;
  const int phi_direction_2x2 = phi_shift < 0.0F ? -1 : 1;

  // Intentional local implementation instead of blindly consuming the
  // positional vector returned by RawCluster::get_shower_shapes(). This keeps
  // the 2x2 definition explicit, validates the denominator above, and applies
  // phi wrapping before the cluster-tower energy lookup at iphi=0/255.
  const auto cluster_tower_energy = [&](int ieta, int iphi) -> float
  {
    if (ieta < 0 || ieta >= kCemcEtaBins)
    {
      return 0.0F;
    }
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, wrap_phi_index(iphi));
    const auto iter = tower_map.find(key);
    if (iter == tower_map.end() || !std::isfinite(iter->second) || iter->second <= config_.min_tower_energy)
    {
      return 0.0F;
    }
    return iter->second;
  };

  // Match RawClusterv1::get_shower_shapes and the reference CaloAna24 code:
  // E1...E4 come from the cluster towermap (including split fractions), while
  // the 7x7 window variables below use full calibrated TowerInfo energies.
  const float e1 = cluster_tower_energy(center_ieta, center_iphi);
  const float e2 = cluster_tower_energy(center_ieta, center_iphi + phi_direction_2x2);
  const float e3 = cluster_tower_energy(center_ieta + eta_direction, center_iphi + phi_direction_2x2);
  const float e4 = cluster_tower_energy(center_ieta + eta_direction, center_iphi);
  const float inverse_cluster_energy = 1.0F / result.cluster_energy_above_threshold;
  result.et1 = (e1 + e2 + e3 + e4) * inverse_cluster_energy;
  result.et2 = (e1 + e2 - e3 - e4) * inverse_cluster_energy;
  result.et3 = (e1 - e2 - e3 + e4) * inverse_cluster_energy;
  result.et4 = e3 * inverse_cluster_energy;

  std::array<float, kPatchSize> selected_energy = {};
  for (int eta_offset = -3; eta_offset <= 3; ++eta_offset)
  {
    for (int phi_offset = -3; phi_offset <= 3; ++phi_offset)
    {
      const std::size_t index = patch_index(eta_offset, phi_offset);
      const int ieta = center_ieta + eta_offset;
      const int iphi = wrap_phi_index(center_iphi + phi_offset);

      // Intentional extension beyond CaloAna24: the reference skips clusters
      // whose 7x7 patch crosses an eta edge. We zero-pad the non-instrumented
      // cells instead so every input cluster keeps the same TTree vector index.
      if (ieta < 0 || ieta >= kCemcEtaBins)
      {
        result.full_containment = false;
        result.edge_padded = true;
        continue;
      }

      const unsigned int tower_key = TowerInfoDefs::encode_emcal(static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));
      result.patch_owned[index] = owned_tower_keys.contains(tower_key) ? 1U : 0U;

      TowerInfo *tower = cemc_towers.get_tower_at_key(tower_key);
      if (!tower)
      {
        result.tower_data_complete = false;
        continue;
      }

      const float energy = tower->get_energy();
      result.patch_energy[index] = energy;
      const bool is_good = tower->get_isGood() && std::isfinite(energy);
      result.patch_good[index] = is_good ? 1U : 0U;
      if (is_good && energy > config_.min_tower_energy)
      {
        selected_energy[index] = energy;
      }
    }
  }

  float w_eta_cogx_sum = 0.0F;
  float w_phi_cogx_sum = 0.0F;

  // CaloAna24 chooses the second phi column on the CoG side. At exactly
  // zero phi shift its strict comparison selects the negative-phi column.
  const int phi_direction_3x2 = phi_shift > 0.0F ? 1 : -1;
  for (int eta_offset = -3; eta_offset <= 3; ++eta_offset)
  {
    for (int phi_offset = -3; phi_offset <= 3; ++phi_offset)
    {
      const std::size_t index = patch_index(eta_offset, phi_offset);
      const float energy = selected_energy[index];
      const int abs_eta_offset = std::abs(eta_offset);
      const int abs_phi_offset = std::abs(phi_offset);

      // Match the reference ownership semantics. A tower appearing in more
      // than one cluster is owned by every such cluster; no exclusive claim or
      // fractional reweighting of the calibrated TowerInfo energy is applied.
      if (result.patch_owned[index] != 0U)
      {
        result.owned_patch_energy += energy;
        if (eta_offset != 0 || phi_offset != 0)
        {
          const float deta = static_cast<float>(eta_offset) - eta_shift;
          const float dphi = static_cast<float>(phi_offset) - phi_shift;
          w_eta_cogx_sum += energy * deta * deta;
          w_phi_cogx_sum += energy * dphi * dphi;
        }
      }

      // These fixed-window sums intentionally do not require cluster
      // ownership, matching CaloAna24's E77-based implementation.
      if (eta_offset == 0 && phi_offset == 0)
      {
        result.e11 = energy;
      }
      if (abs_eta_offset <= 1 && abs_phi_offset <= 1)
      {
        result.e33 += energy;
      }
      if (abs_eta_offset <= 1 && abs_phi_offset <= 2)
      {
        result.e35 += energy;
      }
      if (abs_eta_offset <= 1 && (phi_offset == 0 || phi_offset == phi_direction_3x2))
      {
        result.e32 += energy;
      }
    }
  }

  if (result.owned_patch_energy > 0.0F)
  {
    result.w_eta_cogx = w_eta_cogx_sum / result.owned_patch_energy;
    result.w_phi_cogx = w_phi_cogx_sum / result.owned_patch_energy;
  }
  if (result.e33 > 0.0F)
  {
    result.e11_over_e33 = result.e11 / result.e33;
  }
  if (result.e35 > 0.0F)
  {
    result.e32_over_e35 = result.e32 / result.e35;
  }

  result.valid = result.tower_data_complete &&
                 result.owned_patch_energy > 0.0F &&
                 result.e33 > 0.0F &&
                 result.e35 > 0.0F;
  return result;
}
