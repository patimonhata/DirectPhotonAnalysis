#include "CemcClusterDumper.h"

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerDefs.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>

namespace
{
  constexpr float invalid_float = -999.0f;
}

CemcClusterDumper::CemcClusterDumper(const std::string& name, int run, std::string job_index, bool save_tree)
  : SubsysReco(name)
  , m_run(run)
  , m_job_index(std::move(job_index))
  , m_save_tree(save_tree)
{
  m_output_dir = std::string("/sphenix/user/ryotaro/CemcClusterDumper/output/") + std::to_string(m_run);
  m_output_file = m_output_dir + "/" + std::to_string(m_run) + "-" + m_job_index + ".root";
}

CemcClusterDumper::~CemcClusterDumper()
{
  std::cout << "CemcClusterDumper::~CemcClusterDumper() Calling the deconstructor" << std::endl;
}

int CemcClusterDumper::Init(PHCompositeNode* topNode)
{
  std::cout << "CemcClusterDumper::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  if (m_save_tree) {
    gSystem->mkdir(m_output_dir.c_str(), true);
    m_out_file = new TFile(m_output_file.c_str(), "RECREATE");
    m_tree = new TTree("cemc_clusters", "CEMC clusters per event");

    m_tree->Branch("event", &m_event);
    m_tree->Branch("has_vertex", &m_has_vertex);
    m_tree->Branch("vtx_x", &m_vtx_x);
    m_tree->Branch("vtx_y", &m_vtx_y);
    m_tree->Branch("vtx_z", &m_vtx_z);

    m_tree->Branch("cluster_id", &m_cluster_id);
    m_tree->Branch("cluster_energy", &m_cluster_energy);
    m_tree->Branch("cluster_ecore", &m_cluster_ecore);
    m_tree->Branch("cluster_chi2", &m_cluster_chi2);
    m_tree->Branch("cluster_prob", &m_cluster_prob);
    m_tree->Branch("cluster_merged_prob", &m_cluster_merged_prob);
    m_tree->Branch("cluster_et_iso", &m_cluster_et_iso);
    m_tree->Branch("cluster_mean_time", &m_cluster_mean_time);

    m_tree->Branch("cluster_x", &m_cluster_x);
    m_tree->Branch("cluster_y", &m_cluster_y);
    m_tree->Branch("cluster_z", &m_cluster_z);
    m_tree->Branch("cluster_r", &m_cluster_r);
    m_tree->Branch("cluster_phi_det", &m_cluster_phi_det);

    m_tree->Branch("cluster_eta0", &m_cluster_eta0);
    m_tree->Branch("cluster_phi0", &m_cluster_phi0);
    m_tree->Branch("cluster_et0", &m_cluster_et0);
    m_tree->Branch("cluster_etavtx", &m_cluster_etavtx);
    m_tree->Branch("cluster_phivtx", &m_cluster_phivtx);
    m_tree->Branch("cluster_etvtx", &m_cluster_etvtx);

    m_tree->Branch("cluster_n_towers", &m_cluster_n_towers);
    m_tree->Branch("cluster_lead_tower_ieta", &m_cluster_lead_tower_ieta);
    m_tree->Branch("cluster_lead_tower_iphi", &m_cluster_lead_tower_iphi);

    m_tree->Branch("cluster_x_tower_raw", &m_cluster_x_tower_raw);
    m_tree->Branch("cluster_y_tower_raw", &m_cluster_y_tower_raw);
    m_tree->Branch("cluster_x_tower_corr", &m_cluster_x_tower_corr);
    m_tree->Branch("cluster_y_tower_corr", &m_cluster_y_tower_corr);

    m_tree->Branch("member_cluster_index", &m_member_cluster_index);
    m_tree->Branch("member_cluster_id", &m_member_cluster_id);
    m_tree->Branch("member_tower_key", &m_member_tower_key);
    m_tree->Branch("member_tower_caloid", &m_member_tower_caloid);
    m_tree->Branch("member_tower_ieta", &m_member_tower_ieta);
    m_tree->Branch("member_tower_iphi", &m_member_tower_iphi);
    m_tree->Branch("member_energy", &m_member_energy);
    m_tree->Branch("member_energy_fraction", &m_member_energy_fraction);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CemcClusterDumper::InitRun(PHCompositeNode* topNode)
{
  std::cout << "CemcClusterDumper::InitRun(PHCompositeNode *topNode) Initializing for Run XXX... " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CemcClusterDumper::process_event(PHCompositeNode* topNode)
{
  resetEventBuffers();
  m_event++;
  if (m_event % 20 == 0) {
    std::cout << "event: " << m_event << std::endl;
  }

  auto* clusters = findNode::getClass<RawClusterContainer>(topNode, m_cluster_node_name);
  if (!clusters) {
    std::cout << "CemcClusterDumper::process_event Missing " << m_cluster_node_name << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  float vx = 0.0f;
  float vy = 0.0f;
  float vz = 0.0f;
  const bool has_vtx = getVertexXyz(topNode, vx, vy, vz);
  m_has_vertex = has_vtx ? 1 : 0;
  m_vtx_x = has_vtx ? vx : invalid_float;
  m_vtx_y = has_vtx ? vy : invalid_float;
  m_vtx_z = has_vtx ? vz : invalid_float;

  const unsigned int n_clusters = clusters->size();
  if (n_clusters == 0) {
    std::cout << "event: " << m_event << " n_clusters was 0." << std::endl;
  }

  m_cluster_id.reserve(n_clusters);
  m_cluster_energy.reserve(n_clusters);
  m_cluster_ecore.reserve(n_clusters);
  m_cluster_chi2.reserve(n_clusters);
  m_cluster_prob.reserve(n_clusters);
  m_cluster_merged_prob.reserve(n_clusters);
  m_cluster_et_iso.reserve(n_clusters);
  m_cluster_mean_time.reserve(n_clusters);
  m_cluster_x.reserve(n_clusters);
  m_cluster_y.reserve(n_clusters);
  m_cluster_z.reserve(n_clusters);
  m_cluster_r.reserve(n_clusters);
  m_cluster_phi_det.reserve(n_clusters);
  m_cluster_eta0.reserve(n_clusters);
  m_cluster_phi0.reserve(n_clusters);
  m_cluster_et0.reserve(n_clusters);
  m_cluster_etavtx.reserve(n_clusters);
  m_cluster_phivtx.reserve(n_clusters);
  m_cluster_etvtx.reserve(n_clusters);
  m_cluster_n_towers.reserve(n_clusters);
  m_cluster_lead_tower_ieta.reserve(n_clusters);
  m_cluster_lead_tower_iphi.reserve(n_clusters);
  m_cluster_x_tower_raw.reserve(n_clusters);
  m_cluster_y_tower_raw.reserve(n_clusters);
  m_cluster_x_tower_corr.reserve(n_clusters);
  m_cluster_y_tower_corr.reserve(n_clusters);

  RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
  for (RawClusterContainer::ConstIterator iter = cluster_range.first; iter != cluster_range.second; ++iter) {
    const RawCluster* cluster = iter->second;
    if (!cluster) {
      continue;
    }

    const int cluster_index = static_cast<int>(m_cluster_id.size());
    const unsigned int cluster_id = cluster->get_id();
    const float energy = cluster->get_energy();
    const float x = cluster->get_x();
    const float y = cluster->get_y();
    const float z = cluster->get_z();

    const TVector3 pos0(x, y, z);
    const float eta0 = pos0.Eta();
    const float phi0 = pos0.Phi();
    const float et0 = energy / std::cosh(eta0);

    float etavtx = invalid_float;
    float phivtx = invalid_float;
    float etvtx = invalid_float;
    if (has_vtx) {
      const TVector3 posv(x - vx, y - vy, z - vz);
      etavtx = posv.Eta();
      phivtx = posv.Phi();
      etvtx = energy / std::cosh(etavtx);
    }

    const std::pair<int, int> lead_tower = cluster->get_lead_tower();

    m_cluster_id.push_back(cluster_id);
    m_cluster_energy.push_back(energy);
    m_cluster_ecore.push_back(cluster->get_ecore());
    m_cluster_chi2.push_back(cluster->get_chi2());
    m_cluster_prob.push_back(cluster->get_prob());
    m_cluster_merged_prob.push_back(cluster->get_merged_cluster_prob());
    m_cluster_et_iso.push_back(cluster->get_et_iso());
    m_cluster_mean_time.push_back(cluster->mean_time());

    m_cluster_x.push_back(x);
    m_cluster_y.push_back(y);
    m_cluster_z.push_back(z);
    m_cluster_r.push_back(cluster->get_r());
    m_cluster_phi_det.push_back(cluster->get_phi());

    m_cluster_eta0.push_back(eta0);
    m_cluster_phi0.push_back(phi0);
    m_cluster_et0.push_back(et0);
    m_cluster_etavtx.push_back(etavtx);
    m_cluster_phivtx.push_back(phivtx);
    m_cluster_etvtx.push_back(etvtx);

    m_cluster_n_towers.push_back(static_cast<int>(cluster->getNTowers()));
    m_cluster_lead_tower_ieta.push_back(lead_tower.first);
    m_cluster_lead_tower_iphi.push_back(lead_tower.second);

    m_cluster_x_tower_raw.push_back(cluster->x_tower_raw());
    m_cluster_y_tower_raw.push_back(cluster->y_tower_raw());
    m_cluster_x_tower_corr.push_back(cluster->x_tower_corr());
    m_cluster_y_tower_corr.push_back(cluster->y_tower_corr());

    for (const auto& tower_entry : cluster->get_towermap()) {
      const RawTowerDefs::keytype tower_key = tower_entry.first;
      const float member_energy = tower_entry.second;
      float member_fraction = 0.0f;
      if (std::isfinite(energy) && energy != 0.0f) {
        member_fraction = member_energy / energy;
      }

      m_member_cluster_index.push_back(cluster_index);
      m_member_cluster_id.push_back(cluster_id);
      m_member_tower_key.push_back(tower_key);
      m_member_tower_caloid.push_back(static_cast<int>(RawTowerDefs::decode_caloid(tower_key)));
      m_member_tower_ieta.push_back(static_cast<int>(RawTowerDefs::decode_index1(tower_key)));
      m_member_tower_iphi.push_back(static_cast<int>(RawTowerDefs::decode_index2(tower_key)));
      m_member_energy.push_back(member_energy);
      m_member_energy_fraction.push_back(member_fraction);
    }
  }

  if (m_save_tree && m_tree) {
    m_tree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CemcClusterDumper::ResetEvent(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CemcClusterDumper::Reset(PHCompositeNode* topNode)
{
  std::cout << "CemcClusterDumper::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CemcClusterDumper::EndRun(const int runnumber)
{
  std::cout << "CemcClusterDumper::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CemcClusterDumper::End(PHCompositeNode* topNode)
{
  std::cout << "CemcClusterDumper::End(PHCompositeNode *topNode) This is the End... " << std::endl;

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

void CemcClusterDumper::resetEventBuffers()
{
  m_has_vertex = 0;
  m_vtx_x = 0.0f;
  m_vtx_y = 0.0f;
  m_vtx_z = 0.0f;

  m_cluster_id.clear();
  m_cluster_energy.clear();
  m_cluster_ecore.clear();
  m_cluster_chi2.clear();
  m_cluster_prob.clear();
  m_cluster_merged_prob.clear();
  m_cluster_et_iso.clear();
  m_cluster_mean_time.clear();
  m_cluster_x.clear();
  m_cluster_y.clear();
  m_cluster_z.clear();
  m_cluster_r.clear();
  m_cluster_phi_det.clear();
  m_cluster_eta0.clear();
  m_cluster_phi0.clear();
  m_cluster_et0.clear();
  m_cluster_etavtx.clear();
  m_cluster_phivtx.clear();
  m_cluster_etvtx.clear();
  m_cluster_n_towers.clear();
  m_cluster_lead_tower_ieta.clear();
  m_cluster_lead_tower_iphi.clear();
  m_cluster_x_tower_raw.clear();
  m_cluster_y_tower_raw.clear();
  m_cluster_x_tower_corr.clear();
  m_cluster_y_tower_corr.clear();

  m_member_cluster_index.clear();
  m_member_cluster_id.clear();
  m_member_tower_key.clear();
  m_member_tower_caloid.clear();
  m_member_tower_ieta.clear();
  m_member_tower_iphi.clear();
  m_member_energy.clear();
  m_member_energy_fraction.clear();
}

bool CemcClusterDumper::getVertexXyz(PHCompositeNode* topNode, float& vx, float& vy, float& vz) const
{
  auto* gvtx_map = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!gvtx_map || gvtx_map->empty()) {
    std::cout << "In CemcClusterDumper::GetVertexXYZ, gvtx_map not found." << std::endl;
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
