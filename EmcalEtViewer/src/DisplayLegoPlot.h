// Fun4All analysis module for INTT hit carryover monitor, made by Ryotaro Koike 20250808
#ifndef RYOTARO_DISPLAY_LEGO_PLOT_H_20260224
#define RYOTARO_DISPLAY_LEGO_PLOT_H_20260224

#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TTree.h>


#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <globalvertex/GlobalVertexMap.h>
#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/getClass.h>

class PHCompositeNode;
class CaloEvalStack;
class CaloRawClusterEval;
class CaloTruthEval;
class RawTowerGeom;
class RawTowerGeomContainer;
class TowerInfoContainer;
class PHG4TruthInfoContainer;
class PHG4Particle;
class CaloEvalStack;

namespace HepMC
{
  class GenEvent;
}

class DisplayLegoPlot : public SubsysReco {
  public:
    DisplayLegoPlot(const std::string &name, int run, std::string job_index, bool save_tree);
    ~DisplayLegoPlot() override;
    
    // mandatory methods for Fun4All analysis module
    int Init(PHCompositeNode *topNode) override;
    int InitRun(PHCompositeNode *topNode) override;
    int process_event(PHCompositeNode *topNode) override;
    int ResetEvent(PHCompositeNode *topNode) override;
    int Reset(PHCompositeNode *topNode) override;
    int EndRun(const int runnumber) override;
    int End(PHCompositeNode *topNode) override;
    // void Print(const std::string &what) const override;



  private:
    void reset_event_buffers();
    bool get_vertex_xyz(PHCompositeNode *topNode, float &vx, float &vy, float &vz) const;

    int run_ = -1;
    std::string job_index_;
    bool save_tree_ = false;
    std::string output_dir_;
    std::string output_file_;

    TFile *out_file_ = nullptr;
    TTree *tree_ = nullptr;

    int event_ = 0;
    float vtx_x_ = 0.0f;
    float vtx_y_ = 0.0f;
    float vtx_z_ = 0.0f;

    std::vector<float> eta0_;
    std::vector<float> etavtx_;
    std::vector<float> phi0_;
    std::vector<float> phivtx_;
    std::vector<float> et0_;
    std::vector<float> etvtx_;
    std::vector<float> energy_;
    std::vector<int> ieta_;
    std::vector<int> iphi_;
};



#endif // RYOTARO_DISPLAY_LEGO_PLOT_H_20260224
