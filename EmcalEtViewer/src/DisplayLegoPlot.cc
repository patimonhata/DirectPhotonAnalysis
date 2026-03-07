#include "DisplayLegoPlot.h"

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>

DisplayLegoPlot::DisplayLegoPlot(const std::string& name, int run, std::string job_index, bool save_tree)
  : SubsysReco(name)
  , m_run(run)
  , m_job_index(std::move(job_index))
  , m_save_tree(save_tree)
{
  m_output_dir = std::string("/sphenix/user/ryotaro/DisplayLegoPlot/output/") + std::to_string(m_run);
  m_output_file = m_output_dir + "/" + std::to_string(m_run) + "-" + m_job_index + ".root";

}

DisplayLegoPlot::~DisplayLegoPlot()
{
  std::cout << "DisplayLegoPlot::~DisplayLegoPlot() Calling the deconstructor" << std::endl;
}


int DisplayLegoPlot::Init(PHCompositeNode* topNode)
{
  std::cout << "DisplayLegoPlot::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  if (m_save_tree) {
    gSystem->mkdir(m_output_dir.c_str(), true);
    m_out_file = new TFile(m_output_file.c_str(), "RECREATE");
    m_tree = new TTree("cemc_et", "CEMC E_T per event");

    m_tree->Branch("event", &m_event);
    m_tree->Branch("vtx_x", &m_vtx_x);
    m_tree->Branch("vtx_y", &m_vtx_y);
    m_tree->Branch("vtx_z", &m_vtx_z);

    m_tree->Branch("eta0", &m_eta0);
    m_tree->Branch("etavtx", &m_etavtx);
    m_tree->Branch("phi0", &m_phi0);
    m_tree->Branch("phivtx", &m_phivtx);
    m_tree->Branch("et0", &m_et0);
    m_tree->Branch("etvtx", &m_etvtx);
    m_tree->Branch("energy", &m_energy);
    m_tree->Branch("ieta", &m_ieta);
    m_tree->Branch("iphi", &m_iphi);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int DisplayLegoPlot::InitRun(PHCompositeNode* topNode)
{
  std::cout << "DisplayLegoPlot::InitRun(PHCompositeNode *topNode) Initializing for Run XXX... " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int DisplayLegoPlot::process_event(PHCompositeNode* topNode)
{
  resetEventBuffers();
  m_event++;
  if (m_event % 20 == 0) {
    std::cout << "event: " << m_event << std::endl;
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
  const bool has_vtx = getVertexXyz(topNode, vx, vy, vz);
  m_vtx_x = has_vtx ? vx : -999.0f;
  m_vtx_y = has_vtx ? vy : -999.0f;
  m_vtx_z = has_vtx ? vz : -999.0f;

  const std::size_t n_towers = towers->size();
  if (n_towers == 0) {
    std::cout << "event: " << m_event << " n_towers was 0." << std::endl;
  }
  m_eta0.reserve(n_towers);
  m_etavtx.reserve(n_towers);
  m_phi0.reserve(n_towers);
  m_phivtx.reserve(n_towers);
  m_et0.reserve(n_towers);
  m_etvtx.reserve(n_towers);
  m_energy.reserve(n_towers);
  m_ieta.reserve(n_towers);
  m_iphi.reserve(n_towers);

  for (std::size_t ich = 0; ich < n_towers; ++ich) {
    auto* tower = towers->get_tower_at_channel(static_cast<unsigned int>(ich));
    if (!tower) {
      std::cout << "tower not found" << std::endl;
      continue;
    }

    const float energy = tower->get_energy();
    if (!(energy > 0.0f)) {
      // std::cout << "energy is zero" << std::endl;
      continue;
    }

    const unsigned int tower_key = towers->encode_key(static_cast<unsigned int>(ich));  /* TowerInfo の方の key */
    int ieta = towers->getTowerEtaBin(tower_key); /* tower がカバーする eta 領域のビン番号 */
    int iphi = towers->getTowerPhiBin(tower_key);
    RawTowerDefs::CalorimeterId caloid = RawTowerDefs::CalorimeterId::CEMC;
    RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(caloid, ieta, iphi); /* RawTower の方の key */
    RawTowerGeom* tower_geom = geom->get_tower_geometry(key);
    if (!tower_geom) {
      continue;
    }

    const double x = tower_geom->get_center_x();
    const double y = tower_geom->get_center_y();
    const double z = tower_geom->get_center_z();

    const TVector3 pos0(x, y, z);
    const TVector3 posv(x - m_vtx_x, y - m_vtx_y, z - m_vtx_z);

    const float eta0 = pos0.Eta();
    const float etavtx = posv.Eta();
    const float phi0 = pos0.Phi();
    const float phivtx = posv.Phi();

    const float et0 = energy / std::cosh(eta0);
    const float etvtx = energy / std::cosh(etavtx);

    m_eta0.push_back(eta0);
    m_etavtx.push_back(etavtx);
    m_phi0.push_back(phi0);
    m_phivtx.push_back(phivtx);
    m_et0.push_back(et0);
    m_etvtx.push_back(etvtx);
    m_energy.push_back(energy);
    m_ieta.push_back(tower_geom->get_bineta());
    m_iphi.push_back(tower_geom->get_binphi());
  }

  if (m_save_tree && m_tree) {
    m_tree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int DisplayLegoPlot::ResetEvent(PHCompositeNode* topNode)
{
  // std::cout << "DisplayLegoPlot::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int DisplayLegoPlot::Reset(PHCompositeNode* topNode)
{
  std::cout << "DisplayLegoPlot::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int DisplayLegoPlot::EndRun(const int runnumber)
{
  std::cout << "DisplayLegoPlot::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int DisplayLegoPlot::End(PHCompositeNode* topNode)
{
  std::cout << "DisplayLegoPlot::End(PHCompositeNode *topNode) This is the End... " << std::endl;

  if (m_save_tree && m_out_file) {
    m_out_file->cd();
    if (m_tree) {
      m_tree->Write();
    }
    m_out_file->Close();
    delete m_out_file;
    m_out_file = nullptr;
    m_tree = nullptr;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void DisplayLegoPlot::resetEventBuffers()
{
  m_eta0.clear();
  m_etavtx.clear();
  m_phi0.clear();
  m_phivtx.clear();
  m_et0.clear();
  m_etvtx.clear();
  m_energy.clear();
  m_ieta.clear();
  m_iphi.clear();
}

bool DisplayLegoPlot::getVertexXyz(PHCompositeNode* topNode, float& vx, float& vy, float& vz) const
{
  auto* gvtx_map = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!gvtx_map || gvtx_map->empty()) {
    std::cout << "In DisplayLegoPlot::GetVertexXYZ, gvtx_map not found." << std::endl;
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

  std::vector<GlobalVertex*> candidates = gvtx_map->get_gvtxs_with_type(prefer_types);
  const GlobalVertex* vtx = nullptr;
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
