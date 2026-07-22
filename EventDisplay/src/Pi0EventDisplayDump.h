#ifndef RYOTARO_PI0EVENTDISPLAYDUMP_H_20260611
#define RYOTARO_PI0EVENTDISPLAYDUMP_H_20260611

#include <fun4all/SubsysReco.h>

#include <calobase/RawTowerDefs.h>

#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class RawTowerGeom;
class RawTowerGeomContainer;
class TowerInfo;
class TowerInfoContainer;
class TFile;
class TTree;

class Pi0EventDisplayDump : public SubsysReco
{
 public:
  explicit Pi0EventDisplayDump(const std::string& name = "Pi0EventDisplayDump");
  ~Pi0EventDisplayDump() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void set_output_file(const std::string& filename){ m_outputFileName = filename; };

  void set_truth_node(const std::string& name){ m_truthNode = name; };
  void set_cluster_node(const std::string& name){ m_clusterNode = name; };
  void set_tower_node(const std::string& name){ m_towerNode = name; };
  void set_tower_geom_node(const std::string& name){ m_towerGeomNode = name; };
  void set_cemc_hit_node(const std::string& name){ m_cemcHitNode = name; };

  void set_cluster_energy_min(double e){ m_clusterEnergyMin = e; };
  void set_tower_energy_min(double e){ m_towerEnergyMin = e; };
  void set_event_min(int event){ m_eventMin = event; };
  void set_event_max(int event){ m_eventMax = event; };

  void set_write_hits(bool value){ m_writeHits = value; };
  void set_verbosity(int v){ m_verbosity = v; };

 private:
  static constexpr double kInvalidDouble = -999.0;
  static constexpr int kInvalidInt = -999;
  static constexpr double kDefaultDisplayRadius = 95.0;

  struct TruthProjection
  {
    int track_id = kInvalidInt;
    double eta = kInvalidDouble;
    double phi = kInvalidDouble;
    double x = kInvalidDouble;
    double y = kInvalidDouble;
    double z = kInvalidDouble;
  };

  struct TowerLookup
  {
    TowerInfo* tower = nullptr;
    unsigned int raw_key = 0;
    int ieta = kInvalidInt;
    int iphi = kInvalidInt;
    int channel = kInvalidInt;
  };

  void create_output_directory() const;
  void create_trees();
  void close_output_file();

  void fill_truth(PHG4TruthInfoContainer* truth, RawTowerGeomContainer* geom, std::vector<TruthProjection>& gamma_projections);
  void fill_truth_particle(PHG4TruthInfoContainer* truth, PHG4Particle* particle);
  void fill_truth_segment(PHG4TruthInfoContainer* truth, PHG4Particle* particle, double display_radius);
  void fill_clusters(PHCompositeNode* topNode, const std::vector<TruthProjection>& gamma_projections, const std::map<unsigned int, TowerLookup>& tower_lookup);
  void fill_towers(TowerInfoContainer* towers, RawTowerGeomContainer* geom, std::map<unsigned int, TowerLookup>& tower_lookup);
  void fill_hits(PHCompositeNode* topNode, PHG4TruthInfoContainer* truth);

  int find_ancestor(PHG4TruthInfoContainer* truth, int track_id, int pid, int* generation = nullptr) const;
  bool has_direct_parent_pid(PHG4TruthInfoContainer* truth, const PHG4Particle* particle, int pid) const;
  bool get_vertex(PHG4TruthInfoContainer* truth, int vtx_id, double& x, double& y, double& z) const;
  bool project_to_radius(double x0, double y0, double z0, double px, double py, double pz, double radius, double& x1, double& y1, double& z1) const;
  bool find_first_daughter_vertex(PHG4TruthInfoContainer* truth, int track_id, double& x, double& y, double& z) const;
  double get_display_radius(RawTowerGeomContainer* geom) const;
  RawTowerGeom* get_tower_geom(RawTowerGeomContainer* geom, unsigned int raw_key) const;
  static double radius(double x, double y);
  static double eta_from_xyz(double x, double y, double z);
  static double phi_from_xy(double x, double y);
  static double delta_phi(double a, double b);
  static double finite_or_invalid(double value);
  static int finite_or_invalid_int(double value);

  std::string m_outputFileName = "event_display.root";
  std::string m_truthNode = "G4TruthInfo";
  std::string m_clusterNode = "CLUSTERINFO_CEMC";
  std::string m_towerNode = "TOWERINFO_CALIB_CEMC";
  std::string m_towerGeomNode = "TOWERGEOM_CEMC";
  std::string m_cemcHitNode = "G4HIT_CEMC";

  double m_clusterEnergyMin = 0.0;
  double m_towerEnergyMin = 0.0;
  int m_eventMin = 0;
  int m_eventMax = -1;
  bool m_writeHits = true;
  int m_verbosity = 0;

  TFile* m_outputFile = nullptr;
  TTree* m_eventsTree = nullptr;
  TTree* m_truthParticlesTree = nullptr;
  TTree* m_truthSegmentsTree = nullptr;
  TTree* m_cemcClustersTree = nullptr;
  TTree* m_cemcClusterTowersTree = nullptr;
  TTree* m_cemcTowersTree = nullptr;
  TTree* m_cemcHitsTree = nullptr;

  int m_eventCounter = 0;
  int m_missingTruthWarnings = 0;
  int m_missingClusterWarnings = 0;
  int m_missingTowerWarnings = 0;
  int m_missingGeomWarnings = 0;
  int m_missingHitWarnings = 0;

  int b_event = 0;
  int b_n_truth_particles = 0;
  int b_n_truth_pi0 = 0;
  int b_n_truth_gamma = 0;
  int b_n_clusters = 0;
  int b_n_towers_above_threshold = 0;
  int b_n_cemc_hits = 0;

  int b_track_id = kInvalidInt;
  int b_pid = kInvalidInt;
  int b_parent_id = kInvalidInt;
  int b_primary_id = kInvalidInt;
  int b_vtx_id = kInvalidInt;
  double b_px = kInvalidDouble;
  double b_py = kInvalidDouble;
  double b_pz = kInvalidDouble;
  double b_e = kInvalidDouble;
  double b_vx = kInvalidDouble;
  double b_vy = kInvalidDouble;
  double b_vz = kInvalidDouble;
  double b_pt = kInvalidDouble;
  double b_p = kInvalidDouble;
  double b_eta = kInvalidDouble;
  double b_phi = kInvalidDouble;
  int b_is_primary = 0;
  int b_is_pi0 = 0;
  int b_is_gamma = 0;
  int b_is_pi0_daughter = 0;
  int b_ancestor_pi0 = kInvalidInt;
  int b_ancestor_gamma = kInvalidInt;
  int b_generation = kInvalidInt;

  double b_x0 = kInvalidDouble;
  double b_y0 = kInvalidDouble;
  double b_z0 = kInvalidDouble;
  double b_x1 = kInvalidDouble;
  double b_y1 = kInvalidDouble;
  double b_z1 = kInvalidDouble;
  double b_r0 = kInvalidDouble;
  double b_r1 = kInvalidDouble;
  double b_eta0 = kInvalidDouble;
  double b_phi0 = kInvalidDouble;
  double b_eta1 = kInvalidDouble;
  double b_phi1 = kInvalidDouble;
  int b_segment_type = kInvalidInt;

  std::string b_cluster_node;
  unsigned int b_cluster_id = 0;
  double b_energy = kInvalidDouble;
  double b_ecore = kInvalidDouble;
  double b_chi2 = kInvalidDouble;
  double b_prob = kInvalidDouble;
  double b_merged_cluster_prob = kInvalidDouble;
  double b_x = kInvalidDouble;
  double b_y = kInvalidDouble;
  double b_z = kInvalidDouble;
  double b_r = kInvalidDouble;
  int b_ntowers = 0;
  unsigned int b_lead_tower_key = 0;
  int b_lead_tower_ieta = kInvalidInt;
  int b_lead_tower_iphi = kInvalidInt;
  double b_lead_tower_energy = kInvalidDouble;
  int b_nearest_truth_gamma_track_id = kInvalidInt;
  double b_nearest_truth_gamma_delta_eta = kInvalidDouble;
  double b_nearest_truth_gamma_delta_phi = kInvalidDouble;
  double b_nearest_truth_gamma_delta_r = kInvalidDouble;
  double b_nearest_truth_gamma_delta_tower = kInvalidDouble;

  unsigned int b_tower_key = 0;
  int b_ieta = kInvalidInt;
  int b_iphi = kInvalidInt;
  double b_cluster_tower_value = kInvalidDouble;
  double b_tower_energy = kInvalidDouble;

  std::string b_tower_node;
  int b_channel = kInvalidInt;
  double b_time = kInvalidDouble;
  int b_is_good = kInvalidInt;
  int b_status = kInvalidInt;

  int b_trkid = kInvalidInt;
  double b_edep = kInvalidDouble;
  double b_eion = kInvalidDouble;
  double b_light_yield = kInvalidDouble;
};

#endif
