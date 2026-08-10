#ifndef RYOTARO_PHOTONTREECOMMON_H_20260810
#define RYOTARO_PHOTONTREECOMMON_H_20260810

#include "ShowerShapeCalculator.h"

#include <string>
#include <vector>

class RawCluster;
class RawClusterContainer;
class RawTowerGeomContainer;
class TowerInfoContainer;
class TTree;

namespace photon_tree
{
constexpr double kInvalidDouble = -999.0;
constexpr int kInvalidInt = -999;

struct EventVertex
{
  bool valid = false;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  int source = 0;
};

struct ClusterCollection
{
  unsigned int ncluster = 0;
  unsigned int ntower = 0;

  std::vector<unsigned int> cluster_id;
  std::vector<int> cluster_ntower;
  std::vector<double> cluster_e;
  std::vector<double> cluster_et;
  std::vector<double> cluster_eta;
  std::vector<double> cluster_phi;
  std::vector<double> cluster_x;
  std::vector<double> cluster_y;
  std::vector<double> cluster_z;
  std::vector<double> cluster_px;
  std::vector<double> cluster_py;
  std::vector<double> cluster_pz;

  std::vector<unsigned char> shower_valid;
  std::vector<unsigned char> shower_full_containment;
  std::vector<unsigned char> shower_edge_padded;
  std::vector<unsigned char> shower_tower_data_complete;
  std::vector<float> shower_cog_ieta;
  std::vector<float> shower_cog_iphi;
  std::vector<float> shower_cluster_e_thresholded;
  std::vector<float> shower_owned_patch_e;
  std::vector<float> shower_w_eta_cogx;
  std::vector<float> shower_w_phi_cogx;
  std::vector<float> shower_e11;
  std::vector<float> shower_e33;
  std::vector<float> shower_e32;
  std::vector<float> shower_e35;
  std::vector<float> shower_e11_over_e33;
  std::vector<float> shower_e32_over_e35;
  std::vector<float> shower_et1;
  std::vector<float> shower_et2;
  std::vector<float> shower_et3;
  std::vector<float> shower_et4;
  std::vector<float> shower_patch_e;
  std::vector<unsigned char> shower_patch_good;
  std::vector<unsigned char> shower_patch_owned;

  std::vector<unsigned int> pair_cluster_i;
  std::vector<unsigned int> pair_cluster_j;
  std::vector<double> pair_m_gg;
  std::vector<double> pair_e_asym;

  std::vector<int> tower_cluster_index;
  std::vector<unsigned int> tower_key;
  std::vector<int> tower_ieta;
  std::vector<int> tower_iphi;
  std::vector<double> tower_x;
  std::vector<double> tower_y;
  std::vector<double> tower_z;
  std::vector<double> tower_r;
  std::vector<double> tower_eta;
  std::vector<double> tower_phi;
  std::vector<double> tower_energy;
  std::vector<double> tower_cluster_value;
  std::vector<double> tower_time;
  std::vector<int> tower_is_good;
  std::vector<int> tower_status;

  void clear();
};

struct FilledCollection
{
  ClusterCollection data;
  std::vector<const RawCluster*> ordered_clusters;

  void clear();
};

class PhotonTreeCommon
{
 public:
  PhotonTreeCommon();

  void set_min_cluster_energy(double value) { min_cluster_energy_ = value; }
  void set_shower_shape_min_tower_energy(double value) { shower_shape_min_tower_energy_ = value; }
  void set_store_shower_shape_tower_patch(bool value) { store_shower_shape_tower_patch_ = value; }

  double min_cluster_energy() const { return min_cluster_energy_; }
  double shower_shape_min_tower_energy() const { return shower_shape_min_tower_energy_; }
  int shower_shape_algorithm_version() const { return shower_shape_algorithm_version_; }
  int shower_shape_patch_side() const { return shower_shape_patch_side_; }
  bool store_shower_shape_tower_patch() const { return store_shower_shape_tower_patch_; }

  bool initialize();
  void clear_event();
  bool fill_collection(RawClusterContainer* clusters,
                       TowerInfoContainer* towers,
                       RawTowerGeomContainer* geometry,
                       const EventVertex& vertex,
                       bool include_towers,
                       bool require_unique_tower_keys,
                       FilledCollection& output);
  void create_collection_branches(TTree* tree,
                                  const std::string& prefix,
                                  ClusterCollection& collection,
                                  bool include_towers) const;

  FilledCollection& split() { return split_; }
  FilledCollection& nosplit() { return nosplit_; }
  const FilledCollection& split() const { return split_; }
  const FilledCollection& nosplit() const { return nosplit_; }

 private:
  void append_shower_shape(const ShowerShapeCalculator::Result& result,
                           ClusterCollection& output) const;
  static double radius(double x, double y);
  static double eta_from_xyz(double x, double y, double z);
  static double phi_from_xy(double x, double y);

  double min_cluster_energy_ = 0.0;
  double shower_shape_min_tower_energy_ = 0.070;
  int shower_shape_algorithm_version_ = ShowerShapeCalculator::kAlgorithmVersion;
  int shower_shape_patch_side_ = ShowerShapeCalculator::kPatchSide;
  bool store_shower_shape_tower_patch_ = true;
  ShowerShapeCalculator shower_shape_calculator_;
  FilledCollection split_;
  FilledCollection nosplit_;
};
}

#endif
