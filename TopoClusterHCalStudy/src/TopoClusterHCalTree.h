#ifndef TOPOCLUSTERHCALSTUDY_TOPOCLUSTERHCALTREE_H
#define TOPOCLUSTERHCALSTUDY_TOPOCLUSTERHCALTREE_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;

class TopoClusterHCalTree : public SubsysReco
{
 public:
  explicit TopoClusterHCalTree(const std::string &name = "TopoClusterHCalTree");
  ~TopoClusterHCalTree() override;

  void set_output_file_name(const std::string &value) { output_file_name_ = value; }
  void set_topocluster_node_name(const std::string &value) { topocluster_node_name_ = value; }
  void set_process_id(unsigned int value) { process_id_ = value; }

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  void create_branches();
  void reset_event();
  void close_output_file();

  std::string output_file_name_ = "topocluster_hcal_tree.root";
  std::string topocluster_node_name_ = "TOPOCLUSTER_ALLCALO";

  unsigned int process_id_ = 0;
  unsigned int event_ = 0;
  unsigned int n_topocluster_ = 0;
  float truth_pt_ = 0.0F;
  float truth_eta_ = 0.0F;
  float truth_vertex_z_ = 0.0F;
  float truth_energy_asymmetry_ = -1.0F;

  std::vector<float> emcal_energy_;
  std::vector<float> hcalin_energy_;
  std::vector<float> hcalout_energy_;
  std::vector<float> hcal_total_energy_;

  TFile *output_file_ = nullptr;
  TTree *tree_ = nullptr;
};

#endif
