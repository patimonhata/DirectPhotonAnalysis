// Fun4All analysis module for INTT hit carryover monitor, made by Ryotaro Koike 20250808
#ifndef RYOTARO_Pi0Reconstruction_H_20260210
#define RYOTARO_Pi0Reconstruction_H_20260210

#include <vector>
#include <unordered_map>
#include <set>
#include <bitset>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>


#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainerv2.h>
#include <ffarawobjects/Gl1Packetv3.h>
#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/getClass.h>

class Pi0Reconstruction : public SubsysReco {
  public:
    Pi0Reconstruction(const std::string &name, const int run, const std::string arg_file_path_to_hot_channel_map, std::string job_index, bool save_tree);
    ~Pi0Reconstruction();
    
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

};



#endif // RYOTARO_Pi0Reconstruction_H_20260210
