// Fun4All analysis module for INTT hit carryover monitor, made by Ryotaro Koike 20250808
#ifndef EMCALETVIEWER_SRC_DISPLAYLEGOPLOT_H_
#define EMCALETVIEWER_SRC_DISPLAYLEGOPLOT_H_

#include <string>
#include <vector>

#include <fun4all/SubsysReco.h>

class PHCompositeNode;
class TFile;
class TTree;

class DisplayLegoPlot : public SubsysReco
{
 public:
  DisplayLegoPlot(const std::string& name, int run, std::string job_index, bool save_tree);
  ~DisplayLegoPlot() override;

  // mandatory methods for Fun4All analysis module
  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int ResetEvent(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode* topNode) override;
  // void Print(const std::string &what) const override;

 private:
  void resetEventBuffers();
  bool getVertexXyz(PHCompositeNode* topNode, float& vx, float& vy, float& vz) const;

  int m_run = -1;
  std::string m_job_index;
  bool m_save_tree = false;
  std::string m_output_dir;
  std::string m_output_file;

  TFile* m_out_file = nullptr;
  TTree* m_tree = nullptr;

  int m_event = 0;
  float m_vtx_x = 0.0f;
  float m_vtx_y = 0.0f;
  float m_vtx_z = 0.0f;

  std::vector<float> m_eta0;
  std::vector<float> m_etavtx;
  std::vector<float> m_phi0;
  std::vector<float> m_phivtx;
  std::vector<float> m_et0;
  std::vector<float> m_etvtx;
  std::vector<float> m_energy;
  std::vector<int> m_ieta;
  std::vector<int> m_iphi;
};

#endif  // EMCALETVIEWER_SRC_DISPLAYLEGOPLOT_H_
