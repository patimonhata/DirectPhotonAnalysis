// Fun4All analysis module for CEMC cluster dumping, made by Ryotaro Koike 20250808
#ifndef EMCALETVIEWER_SRC_CEMCCLUSTERDUMPER_H_
#define EMCALETVIEWER_SRC_CEMCCLUSTERDUMPER_H_

#include <string>
#include <utility>
#include <vector>

#include <fun4all/SubsysReco.h>

class PHCompositeNode;
class TFile;
class TTree;

class CemcClusterDumper : public SubsysReco
{
 public:
  CemcClusterDumper(const std::string& name, int run, std::string job_index, bool save_tree);
  ~CemcClusterDumper() override;

  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int ResetEvent(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode* topNode) override;

  void set_output_dir(std::string output_dir) { m_output_dir = std::move(output_dir); };
  void set_output_file(std::string output_file) { m_output_file = std::move(output_file); };
  void set_cluster_node_name(std::string cluster_node_name) { m_cluster_node_name = std::move(cluster_node_name); };

 private:
  void resetEventBuffers();
  bool getVertexXyz(PHCompositeNode* topNode, float& vx, float& vy, float& vz) const;

  int m_run = -1;
  std::string m_job_index;
  bool m_save_tree = false;
  std::string m_output_dir;
  std::string m_output_file;
  std::string m_cluster_node_name = "CLUSTERINFO_CEMC";

  TFile* m_out_file = nullptr;
  TTree* m_tree = nullptr;

  int m_event = 0;
  int m_has_vertex = 0;
  float m_vtx_x = 0.0f;
  float m_vtx_y = 0.0f;
  float m_vtx_z = 0.0f;

  std::vector<unsigned int> m_cluster_id;
  std::vector<float> m_cluster_energy;
  std::vector<float> m_cluster_ecore;
  std::vector<float> m_cluster_chi2;
  std::vector<float> m_cluster_prob;
  std::vector<float> m_cluster_merged_prob;
  std::vector<float> m_cluster_et_iso;
  std::vector<float> m_cluster_mean_time;

  std::vector<float> m_cluster_x;
  std::vector<float> m_cluster_y;
  std::vector<float> m_cluster_z;
  std::vector<float> m_cluster_r;
  std::vector<float> m_cluster_phi_det;

  std::vector<float> m_cluster_eta0;
  std::vector<float> m_cluster_phi0;
  std::vector<float> m_cluster_et0;
  std::vector<float> m_cluster_etavtx;
  std::vector<float> m_cluster_phivtx;
  std::vector<float> m_cluster_etvtx;

  std::vector<int> m_cluster_n_towers;
  std::vector<int> m_cluster_lead_tower_ieta;
  std::vector<int> m_cluster_lead_tower_iphi;

  std::vector<float> m_cluster_x_tower_raw;
  std::vector<float> m_cluster_y_tower_raw;
  std::vector<float> m_cluster_x_tower_corr;
  std::vector<float> m_cluster_y_tower_corr;

  std::vector<int> m_member_cluster_index;
  std::vector<unsigned int> m_member_cluster_id;
  std::vector<unsigned int> m_member_tower_key;
  std::vector<int> m_member_tower_caloid;
  std::vector<int> m_member_tower_ieta;
  std::vector<int> m_member_tower_iphi;
  std::vector<float> m_member_energy;
  std::vector<float> m_member_energy_fraction;
};

#endif  // EMCALETVIEWER_SRC_CEMCCLUSTERDUMPER_H_
