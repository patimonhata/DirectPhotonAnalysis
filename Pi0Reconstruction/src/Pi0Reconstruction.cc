#include "Pi0Reconstruction.h"

Pi0Reconstruction::Pi0Reconstruction(const std::string &name, const int run, const std::string arg_file_path_to_hot_channel_map, std::string job_index, bool save_tree):
{

};

Pi0Reconstruction::~Pi0Reconstruction(){
  std::cout << "Pi0Reconstruction::~Pi0Reconstruction() Calling the deconstructor" << std::endl;
};


int Pi0Reconstruction::Init(PHCompositeNode *topNode){
  std::cout << "Pi0Reconstruction::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
};


int Pi0Reconstruction::InitRun(PHCompositeNode *topNode){
  std::cout << "Pi0Reconstruction::InitRun(PHCompositeNode *topNode) Initializing for Run XXX... " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
};

int Pi0Reconstruction::process_event(PHCompositeNode *topNode){
  if ( my_event_counter_ %200 == 0 ) {
    std::cout<<"eID: " << my_event_counter_ << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
};


int Pi0Reconstruction::ResetEvent(PHCompositeNode *topNode){
  // std::cout << "Pi0Reconstruction::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
};

int Pi0Reconstruction::Reset(PHCompositeNode *topNode){
  std::cout << "Pi0Reconstruction::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
};

int Pi0Reconstruction::EndRun(const int runnumber){
  std::cout << "Pi0Reconstruction::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
};

int Pi0Reconstruction::End(PHCompositeNode *topNode){
  std::cout << "Pi0Reconstruction::End(PHCompositeNode *topNode) This is the End... " << std::endl;


  return Fun4AllReturnCodes::EVENT_OK;
};

// void Pi0Reconstruction::Print(const std::string &what) const {
//   std::cout << "Pi0Reconstruction::Print(const std::string &what) const Printing info for " << what << std::endl;
// }

