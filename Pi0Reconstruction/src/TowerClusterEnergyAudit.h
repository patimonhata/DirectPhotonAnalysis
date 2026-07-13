#ifndef RYOTARO_TOWERCLUSTERENERGYAUDIT_H_20260627
#define RYOTARO_TOWERCLUSTERENERGYAUDIT_H_20260627

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;

class TowerClusterEnergyAudit : public SubsysReco
{
 public:
  explicit TowerClusterEnergyAudit(const std::string &name = "TowerClusterEnergyAudit");
  ~TowerClusterEnergyAudit() override;

  void set_output_file_name(const std::string &output_file_name);
  void set_tower_node_name(const std::string &tower_node_name);
  void set_split_cluster_node_name(const std::string &cluster_node_name);
  void set_no_split_cluster_node_name(const std::string &cluster_node_name);
  void set_tower_energy_threshold(double threshold);
  void set_cluster_energy_threshold(double threshold);
  void set_abort_on_missing_nodes(bool abort_on_missing_nodes);

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  struct EnergySummary
  {
    unsigned int count = 0;
    unsigned int count_above_threshold = 0;
    double energy_sum = 0.0;
    double energy_sum_above_threshold = 0.0;
    double max_energy = 0.0;
    std::vector<double> energies;
  };

  void reset_tree_variables();
  void create_tree_branches();
  void create_output_directory() const;
  bool fill_tower_summary(PHCompositeNode *topNode, EnergySummary &summary);
  bool fill_cluster_summary(PHCompositeNode *topNode, const std::string &node_name, EnergySummary &summary);
  void copy_summary_to_tree(const EnergySummary &summary,
                            unsigned int &count,
                            unsigned int &count_above_threshold,
                            double &energy_sum,
                            double &energy_sum_above_threshold,
                            double &max_energy,
                            std::vector<double> &energies);

  std::string output_file_name_ = "tower_cluster_energy_audit.root";
  std::string tower_node_name_ = "TOWERINFO_CALIBryotaro_CEMC";
  std::string split_cluster_node_name_ = "CLUSTERINFO_CEMC";
  std::string no_split_cluster_node_name_ = "CLUSTERINFO_CEMC_NO_SPLIT";

  double tower_energy_threshold_ = 0.070;
  double cluster_energy_threshold_ = 0.070;
  bool abort_on_missing_nodes_ = false;

  TFile *output_file_ = nullptr;
  TTree *tree_ = nullptr;

  unsigned int event_ = 0;
  unsigned int tower_count_ = 0;
  unsigned int tower_count_above_threshold_ = 0;
  double tower_energy_sum_ = 0.0;
  double tower_energy_sum_above_threshold_ = 0.0;
  double tower_max_energy_ = 0.0;
  std::vector<double> tower_energies_;

  unsigned int split_cluster_count_ = 0;
  unsigned int split_cluster_count_above_threshold_ = 0;
  double split_cluster_energy_sum_ = 0.0;
  double split_cluster_energy_sum_above_threshold_ = 0.0;
  double split_cluster_max_energy_ = 0.0;
  std::vector<double> split_cluster_energies_;

  unsigned int no_split_cluster_count_ = 0;
  unsigned int no_split_cluster_count_above_threshold_ = 0;
  double no_split_cluster_energy_sum_ = 0.0;
  double no_split_cluster_energy_sum_above_threshold_ = 0.0;
  double no_split_cluster_max_energy_ = 0.0;
  std::vector<double> no_split_cluster_energies_;

  unsigned int missing_node_warnings_ = 0;
};

#endif
