#include "PhotonTreeCommon.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>

namespace photon_tree
{
void ClusterCollection::clear()
{
  *this = ClusterCollection{};
}

void FilledCollection::clear()
{
  data.clear();
  ordered_clusters.clear();
}

PhotonTreeCommon::PhotonTreeCommon() = default;

bool PhotonTreeCommon::initialize()
{
  if (min_cluster_energy_ < 0.0 || shower_shape_min_tower_energy_ < 0.0)
  {
    return false;
  }
  ShowerShapeCalculator::Config config;
  config.min_tower_energy = static_cast<float>(shower_shape_min_tower_energy_);
  shower_shape_calculator_ = ShowerShapeCalculator(config);
  return true;
}

void PhotonTreeCommon::clear_event()
{
  split_.clear();
  nosplit_.clear();
}

bool PhotonTreeCommon::fill_collection(RawClusterContainer* clusters,
                                       TowerInfoContainer* towers,
                                       RawTowerGeomContainer* geometry,
                                       const EventVertex& vertex,
                                       bool include_towers,
                                       bool require_unique_tower_keys,
                                       FilledCollection& filled)
{
  filled.clear();
  if (!clusters || !towers || !geometry || !vertex.valid)
  {
    return false;
  }

  ClusterCollection& output = filled.data;
  const auto range = clusters->getClusters();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    const RawCluster* cluster = iter->second;
    if (!cluster)
    {
      continue;
    }
    if (!std::isfinite(cluster->get_energy()) || !std::isfinite(cluster->get_x()) ||
        !std::isfinite(cluster->get_y()) || !std::isfinite(cluster->get_z()))
    {
      return false;
    }
    const bool passes_energy = min_cluster_energy_inclusive_
        ? cluster->get_energy() >= min_cluster_energy_
        : cluster->get_energy() > min_cluster_energy_;
    if (passes_energy)
    {
      filled.ordered_clusters.push_back(cluster);
    }
  }
  std::sort(filled.ordered_clusters.begin(), filled.ordered_clusters.end(),
            [](const RawCluster* lhs, const RawCluster* rhs) {
              if (lhs->get_energy() != rhs->get_energy())
              {
                return lhs->get_energy() > rhs->get_energy();
              }
              return lhs->get_id() < rhs->get_id();
            });

  std::set<unsigned int> seen_tower_keys;
  for (std::size_t cluster_index = 0; cluster_index < filled.ordered_clusters.size(); ++cluster_index)
  {
    const RawCluster* cluster = filled.ordered_clusters[cluster_index];
    const double dx = cluster->get_x() - vertex.x;
    const double dy = cluster->get_y() - vertex.y;
    const double dz = cluster->get_z() - vertex.z;
    const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!(distance > std::numeric_limits<double>::epsilon()))
    {
      return false;
    }

    const double energy = cluster->get_energy();
    const double px = energy * dx / distance;
    const double py = energy * dy / distance;
    const double pz = energy * dz / distance;
    const double pt = std::hypot(px, py);
    if (!(pt > 0.0))
    {
      return false;
    }

    output.cluster_id.push_back(cluster->get_id());
    output.cluster_ntower.push_back(static_cast<int>(cluster->getNTowers()));
    output.cluster_e.push_back(energy);
    output.cluster_et.push_back(pt);
    output.cluster_eta.push_back(std::asinh(pz / pt));
    output.cluster_phi.push_back(std::atan2(py, px));
    output.cluster_x.push_back(cluster->get_x());
    output.cluster_y.push_back(cluster->get_y());
    output.cluster_z.push_back(cluster->get_z());
    output.cluster_px.push_back(px);
    output.cluster_py.push_back(py);
    output.cluster_pz.push_back(pz);
    append_shower_shape(shower_shape_calculator_.calculate(*cluster, *towers), output);

    if (!include_towers)
    {
      continue;
    }
    const auto tower_range = cluster->get_towers();
    for (auto tower_iter = tower_range.first; tower_iter != tower_range.second; ++tower_iter)
    {
      const unsigned int raw_key = tower_iter->first;
      if (require_unique_tower_keys && !seen_tower_keys.insert(raw_key).second)
      {
        std::cout << "PhotonTreeCommon::fill_collection - duplicate tower key "
                  << raw_key << std::endl;
        return false;
      }
      const int ieta = static_cast<int>(RawTowerDefs::decode_index1(raw_key));
      const int iphi = static_cast<int>(RawTowerDefs::decode_index2(raw_key));
      const unsigned int tower_info_key = TowerInfoDefs::encode_emcal(static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));
      TowerInfo* tower = towers->get_tower_at_key(static_cast<int>(tower_info_key));
      RawTowerGeom* tower_geometry = geometry->get_tower_geometry(raw_key);
      if (!tower || !tower_geometry || !std::isfinite(tower->get_energy()))
      {
        return false;
      }
      const double x = tower_geometry->get_center_x();
      const double y = tower_geometry->get_center_y();
      const double z = tower_geometry->get_center_z();
      const double tower_radius = radius(x, y);
      const double eta = eta_from_xyz(x, y, z);
      const double phi = phi_from_xy(x, y);
      if (tower_radius == kInvalidDouble || eta == kInvalidDouble || phi == kInvalidDouble)
      {
        return false;
      }

      output.tower_cluster_index.push_back(static_cast<int>(cluster_index));
      output.tower_key.push_back(raw_key);
      output.tower_ieta.push_back(ieta);
      output.tower_iphi.push_back(iphi);
      output.tower_x.push_back(x);
      output.tower_y.push_back(y);
      output.tower_z.push_back(z);
      output.tower_r.push_back(tower_radius);
      output.tower_eta.push_back(eta);
      output.tower_phi.push_back(phi);
      output.tower_energy.push_back(tower->get_energy());
      output.tower_cluster_value.push_back(
          std::isfinite(tower_iter->second) ? tower_iter->second : kInvalidDouble);
      output.tower_time.push_back(
          std::isfinite(tower->get_time()) ? tower->get_time() : kInvalidDouble);
      output.tower_is_good.push_back(tower->get_isGood() ? 1 : 0);
      output.tower_status.push_back(static_cast<int>(tower->get_status()));
    }
  }

  output.ncluster = static_cast<unsigned int>(output.cluster_id.size());
  output.ntower = static_cast<unsigned int>(output.tower_key.size());
  if (store_cluster_pairs_)
  {
    for (std::size_t i = 0; i < output.ncluster; ++i)
    {
      for (std::size_t j = i + 1; j < output.ncluster; ++j)
      {
        const double total_e = output.cluster_e[i] + output.cluster_e[j];
        const double px = output.cluster_px[i] + output.cluster_px[j];
        const double py = output.cluster_py[i] + output.cluster_py[j];
        const double pz = output.cluster_pz[i] + output.cluster_pz[j];
        const double mass2 = total_e * total_e - px * px - py * py - pz * pz;
        output.pair_cluster_i.push_back(static_cast<unsigned int>(i));
        output.pair_cluster_j.push_back(static_cast<unsigned int>(j));
        output.pair_m_gg.push_back(std::sqrt(std::max(0.0, mass2)));
        output.pair_e_asym.push_back(total_e > 0.0 ? std::abs(output.cluster_e[i] - output.cluster_e[j]) / total_e : kInvalidDouble);
      }
    }
  }
  return true;
}

void PhotonTreeCommon::append_shower_shape(const ShowerShapeCalculator::Result& result,
                                           ClusterCollection& output) const
{
  output.shower_valid.push_back(result.valid ? 1U : 0U);
  output.shower_full_containment.push_back(result.full_containment ? 1U : 0U);
  output.shower_edge_padded.push_back(result.edge_padded ? 1U : 0U);
  output.shower_tower_data_complete.push_back(result.tower_data_complete ? 1U : 0U);
  output.shower_cog_ieta.push_back(result.cog_ieta);
  output.shower_cog_iphi.push_back(result.cog_iphi);
  output.shower_cluster_e_thresholded.push_back(result.cluster_energy_above_threshold);
  output.shower_owned_patch_e.push_back(result.owned_patch_energy);
  output.shower_w_eta_cogx.push_back(result.w_eta_cogx);
  output.shower_w_phi_cogx.push_back(result.w_phi_cogx);
  output.shower_e11.push_back(result.e11);
  output.shower_e33.push_back(result.e33);
  output.shower_e32.push_back(result.e32);
  output.shower_e35.push_back(result.e35);
  output.shower_e11_over_e33.push_back(result.e11_over_e33);
  output.shower_e32_over_e35.push_back(result.e32_over_e35);
  output.shower_et1.push_back(result.et1);
  output.shower_et2.push_back(result.et2);
  output.shower_et3.push_back(result.et3);
  output.shower_et4.push_back(result.et4);
  if (store_shower_shape_tower_patch_)
  {
    output.shower_patch_e.insert(output.shower_patch_e.end(), result.patch_energy.begin(), result.patch_energy.end());
    output.shower_patch_good.insert(output.shower_patch_good.end(), result.patch_good.begin(), result.patch_good.end());
    output.shower_patch_owned.insert(output.shower_patch_owned.end(), result.patch_owned.begin(), result.patch_owned.end());
  }
}

void PhotonTreeCommon::create_collection_branches(TTree* tree,
                                                  const std::string& prefix,
                                                  ClusterCollection& c,
                                                  bool include_towers) const
{
  const auto name = [&prefix](const char* suffix) { return prefix + "_" + suffix; };
  tree->Branch(name("ncluster").c_str(), &c.ncluster);
  tree->Branch(name("cluster_id").c_str(), &c.cluster_id);
  tree->Branch(name("cluster_ntower").c_str(), &c.cluster_ntower);
  tree->Branch(name("cluster_e").c_str(), &c.cluster_e);
  tree->Branch(name("cluster_et").c_str(), &c.cluster_et);
  tree->Branch(name("cluster_eta").c_str(), &c.cluster_eta);
  tree->Branch(name("cluster_phi").c_str(), &c.cluster_phi);
  tree->Branch(name("cluster_x").c_str(), &c.cluster_x);
  tree->Branch(name("cluster_y").c_str(), &c.cluster_y);
  tree->Branch(name("cluster_z").c_str(), &c.cluster_z);
  tree->Branch(name("cluster_px").c_str(), &c.cluster_px);
  tree->Branch(name("cluster_py").c_str(), &c.cluster_py);
  tree->Branch(name("cluster_pz").c_str(), &c.cluster_pz);
  tree->Branch(name("cluster_shower_valid").c_str(), &c.shower_valid);
  tree->Branch(name("cluster_shower_full_containment").c_str(), &c.shower_full_containment);
  tree->Branch(name("cluster_shower_edge_padded").c_str(), &c.shower_edge_padded);
  tree->Branch(name("cluster_shower_tower_data_complete").c_str(), &c.shower_tower_data_complete);
  tree->Branch(name("cluster_shower_cog_ieta").c_str(), &c.shower_cog_ieta);
  tree->Branch(name("cluster_shower_cog_iphi").c_str(), &c.shower_cog_iphi);
  tree->Branch(name("cluster_shower_cluster_e_thresholded").c_str(), &c.shower_cluster_e_thresholded);
  tree->Branch(name("cluster_shower_owned_patch_e").c_str(), &c.shower_owned_patch_e);
  tree->Branch(name("cluster_shower_w_eta_cogx").c_str(), &c.shower_w_eta_cogx);
  tree->Branch(name("cluster_shower_w_phi_cogx").c_str(), &c.shower_w_phi_cogx);
  tree->Branch(name("cluster_shower_e11").c_str(), &c.shower_e11);
  tree->Branch(name("cluster_shower_e33").c_str(), &c.shower_e33);
  tree->Branch(name("cluster_shower_e32").c_str(), &c.shower_e32);
  tree->Branch(name("cluster_shower_e35").c_str(), &c.shower_e35);
  tree->Branch(name("cluster_shower_e11_over_e33").c_str(), &c.shower_e11_over_e33);
  tree->Branch(name("cluster_shower_e32_over_e35").c_str(), &c.shower_e32_over_e35);
  tree->Branch(name("cluster_shower_et1").c_str(), &c.shower_et1);
  tree->Branch(name("cluster_shower_et2").c_str(), &c.shower_et2);
  tree->Branch(name("cluster_shower_et3").c_str(), &c.shower_et3);
  tree->Branch(name("cluster_shower_et4").c_str(), &c.shower_et4);
  tree->Branch(name("cluster_shower_patch_e").c_str(), &c.shower_patch_e);
  tree->Branch(name("cluster_shower_patch_good").c_str(), &c.shower_patch_good);
  tree->Branch(name("cluster_shower_patch_owned").c_str(), &c.shower_patch_owned);
  tree->Branch(name("pair_cluster_i").c_str(), &c.pair_cluster_i);
  tree->Branch(name("pair_cluster_j").c_str(), &c.pair_cluster_j);
  tree->Branch(name("pair_m_gg").c_str(), &c.pair_m_gg);
  tree->Branch(name("pair_e_asym").c_str(), &c.pair_e_asym);

  if (include_towers)
  {
    tree->Branch(name("ntower").c_str(), &c.ntower);
    tree->Branch(name("tower_cluster_index").c_str(), &c.tower_cluster_index);
    tree->Branch(name("tower_key").c_str(), &c.tower_key);
    tree->Branch(name("tower_ieta").c_str(), &c.tower_ieta);
    tree->Branch(name("tower_iphi").c_str(), &c.tower_iphi);
    tree->Branch(name("tower_x").c_str(), &c.tower_x);
    tree->Branch(name("tower_y").c_str(), &c.tower_y);
    tree->Branch(name("tower_z").c_str(), &c.tower_z);
    tree->Branch(name("tower_r").c_str(), &c.tower_r);
    tree->Branch(name("tower_eta").c_str(), &c.tower_eta);
    tree->Branch(name("tower_phi").c_str(), &c.tower_phi);
    tree->Branch(name("tower_energy").c_str(), &c.tower_energy);
    tree->Branch(name("tower_cluster_value").c_str(), &c.tower_cluster_value);
    tree->Branch(name("tower_time").c_str(), &c.tower_time);
    tree->Branch(name("tower_is_good").c_str(), &c.tower_is_good);
    tree->Branch(name("tower_status").c_str(), &c.tower_status);
  }
}

double PhotonTreeCommon::radius(double x, double y)
{
  return std::isfinite(x) && std::isfinite(y) ? std::hypot(x, y) : kInvalidDouble;
}

double PhotonTreeCommon::eta_from_xyz(double x, double y, double z)
{
  const double r = radius(x, y);
  return r > 0.0 && std::isfinite(z) ? std::asinh(z / r) : kInvalidDouble;
}

double PhotonTreeCommon::phi_from_xy(double x, double y)
{
  return std::isfinite(x) && std::isfinite(y) ? std::atan2(y, x) : kInvalidDouble;
}
}
