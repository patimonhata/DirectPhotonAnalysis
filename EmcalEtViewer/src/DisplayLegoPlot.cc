#include "DisplayLegoPlot.h"

#include <TSystem.h>
#include <TVector3.h>

#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfo.h>
#include <globalvertex/GlobalVertex.h>

#include <cmath>
#include <iostream>

DisplayLegoPlot::DisplayLegoPlot(const std::string &name, int run, std::string job_index, bool save_tree)
  : SubsysReco(name)
  , run_(run)
  , job_index_(std::move(job_index))
  , save_tree_(save_tree)
{
  output_dir_ = std::string("/sphenix/user/ryotaro/DisplayLegoPlot/output/") + std::to_string(run_);
  output_file_ = output_dir_ + "/" + std::to_string(run_) + "-" + job_index_ + ".root";

};

DisplayLegoPlot::~DisplayLegoPlot(){
  std::cout << "DisplayLegoPlot::~DisplayLegoPlot() Calling the deconstructor" << std::endl;
};


int DisplayLegoPlot::Init(PHCompositeNode *topNode){
  std::cout << "DisplayLegoPlot::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  if (save_tree_) {
    gSystem->mkdir(output_dir_.c_str(), true);
    out_file_ = new TFile(output_file_.c_str(), "RECREATE");
    tree_ = new TTree("cemc_et", "CEMC E_T per event");

    tree_->Branch("event", &event_);
    tree_->Branch("vtx_x", &vtx_x_);
    tree_->Branch("vtx_y", &vtx_y_);
    tree_->Branch("vtx_z", &vtx_z_);

    tree_->Branch("eta0", &eta0_);
    tree_->Branch("etavtx", &etavtx_);
    tree_->Branch("phi0", &phi0_);
    tree_->Branch("phivtx", &phivtx_);
    tree_->Branch("et0", &et0_);
    tree_->Branch("etvtx", &etvtx_);
    tree_->Branch("energy", &energy_);
    tree_->Branch("ieta", &ieta_);
    tree_->Branch("iphi", &iphi_);
  }

  return Fun4AllReturnCodes::EVENT_OK;
};


int DisplayLegoPlot::InitRun(PHCompositeNode *topNode){
  std::cout << "DisplayLegoPlot::InitRun(PHCompositeNode *topNode) Initializing for Run XXX... " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
};

int DisplayLegoPlot::process_event(PHCompositeNode *topNode){
  reset_event_buffers();
  event_++;
  if (event_ %20 == 0) {
    std::cout << "event: " << event_ << std::endl;
  }

  auto *towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  auto *geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");

  if (!towers || !geom) {
    std::cout << "DisplayLegoPlot::process_event Missing TOWERINFO_CALIB_CEMC or TOWERGEOM_CEMC" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  float vx = 0.0f;
  float vy = 0.0f;
  float vz = 0.0f;
  const bool has_vtx = get_vertex_xyz(topNode, vx, vy, vz);
  vtx_x_ = has_vtx ? vx : -999.0f;
  vtx_y_ = has_vtx ? vy : -999.0f;
  vtx_z_ = has_vtx ? vz : -999.0f;

  const size_t n_towers = towers->size();
  if (n_towers == 0) {
    std::cout << "event: " << event_ << " n_towers was 0." << std::endl;
  }
  eta0_.reserve(n_towers);
  etavtx_.reserve(n_towers);
  phi0_.reserve(n_towers);
  phivtx_.reserve(n_towers);
  et0_.reserve(n_towers);
  etvtx_.reserve(n_towers);
  energy_.reserve(n_towers);
  ieta_.reserve(n_towers);
  iphi_.reserve(n_towers);

  for (unsigned int ich = 0; ich < n_towers; ++ich) {
    TowerInfo *tower = towers->get_tower_at_channel(ich);
    if (!tower) {
      std::cout << "tower not found" << std::endl;
      continue;
    }

    const float energy = tower->get_energy();
    if (!(energy > 0.0f)) {
      // std::cout << "energy is zero" << std::endl;
      continue;
    }

    const unsigned int tower_key = towers->encode_key(ich);  /* TowerInfo の方の key */
    int ieta = towers->getTowerEtaBin(tower_key); /* tower がカバーする eta 領域のビン番号 */
    int iphi = towers->getTowerPhiBin(tower_key);
    RawTowerDefs::CalorimeterId caloid = RawTowerDefs::CalorimeterId::CEMC;
    RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(caloid, ieta, iphi); /* RawTower の方の key */
    RawTowerGeom *tower_geom = geom->get_tower_geometry(key);
    if (!tower_geom) {
      continue;
    }

    const double x = tower_geom->get_center_x();
    const double y = tower_geom->get_center_y();
    const double z = tower_geom->get_center_z();

    const TVector3 pos0(x, y, z);
    const TVector3 posv(x - vtx_x_, y - vtx_y_, z - vtx_z_);

    const float eta0 = pos0.Eta();
    const float etavtx = posv.Eta();
    const float phi0 = pos0.Phi();
    const float phivtx = posv.Phi();

    const float et0 = energy / std::cosh(eta0);
    const float etvtx = energy / std::cosh(etavtx);

    eta0_.push_back(eta0);
    etavtx_.push_back(etavtx);
    phi0_.push_back(phi0);
    phivtx_.push_back(phivtx);
    et0_.push_back(et0);
    etvtx_.push_back(etvtx);
    energy_.push_back(energy);
    ieta_.push_back(tower_geom->get_bineta());
    iphi_.push_back(tower_geom->get_binphi());
  }

  if (save_tree_ && tree_) {
    tree_->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
};


int DisplayLegoPlot::ResetEvent(PHCompositeNode *topNode){
  // std::cout << "DisplayLegoPlot::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
};

int DisplayLegoPlot::Reset(PHCompositeNode *topNode){
  std::cout << "DisplayLegoPlot::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
};

int DisplayLegoPlot::EndRun(const int runnumber){
  std::cout << "DisplayLegoPlot::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
};

int DisplayLegoPlot::End(PHCompositeNode *topNode){
  std::cout << "DisplayLegoPlot::End(PHCompositeNode *topNode) This is the End... " << std::endl;

  if (save_tree_ && out_file_) {
    out_file_->cd();
    if (tree_) {
      tree_->Write();
    }
    out_file_->Close();
    delete out_file_;
    out_file_ = nullptr;
    tree_ = nullptr;
  }

  return Fun4AllReturnCodes::EVENT_OK;
};

void DisplayLegoPlot::reset_event_buffers(){
  eta0_.clear();
  etavtx_.clear();
  phi0_.clear();
  phivtx_.clear();
  et0_.clear();
  etvtx_.clear();
  energy_.clear();
  ieta_.clear();
  iphi_.clear();
}

bool DisplayLegoPlot::get_vertex_xyz(PHCompositeNode *topNode, float &vx, float &vy, float &vz) const{
  auto *gvtx_map = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!gvtx_map || gvtx_map->empty()) {
    std::cout << "In DisplayLegoPlot::get_vertex_xyz, gvtx_map not found." << std::endl;
    return false;
  }

  std::vector<GlobalVertex::VTXTYPE> prefer_types = {
    GlobalVertex::MBD,
    GlobalVertex::SVTX_MBD,
    GlobalVertex::SVTX,
    GlobalVertex::MBD_CALO,
    GlobalVertex::CALO,
    GlobalVertex::TRUTH,
    GlobalVertex::SMEARED
  };

  std::vector<GlobalVertex *> candidates = gvtx_map->get_gvtxs_with_type(prefer_types);
  const GlobalVertex *vtx = nullptr;
  if (!candidates.empty()) {
    vtx = candidates.front();
  } else {
    vtx = gvtx_map->begin() != gvtx_map->end() ? gvtx_map->begin()->second : nullptr;
  }

  if (!vtx) {
    return false;
  }

  vx = vtx->get_x();
  vy = vtx->get_y();
  vz = vtx->get_z();
  return std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz);
}

// void DisplayLegoPlot::Print(const std::string &what) const {
//   std::cout << "DisplayLegoPlot::Print(const std::string &what) const Printing info for " << what << std::endl;
// }
