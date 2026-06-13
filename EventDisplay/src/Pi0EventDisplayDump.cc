#include "Pi0EventDisplayDump.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <utility>

Pi0EventDisplayDump::Pi0EventDisplayDump(const std::string& name): SubsysReco(name) { }

Pi0EventDisplayDump::~Pi0EventDisplayDump() {
  close_output_file();
}

int Pi0EventDisplayDump::Init(PHCompositeNode* /*topNode*/) {
  create_output_directory();

  m_outputFile = TFile::Open(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie()) {
    std::cout << "Pi0EventDisplayDump::Init - failed to open output file " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  create_trees();

  std::cout << "Pi0EventDisplayDump::Init - writing " << m_outputFileName << std::endl;
  std::cout << "Pi0EventDisplayDump::Init - truth node " << m_truthNode << ", cluster node " << m_clusterNode
            << ", tower node " << m_towerNode << ", geometry node " << m_towerGeomNode << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0EventDisplayDump::process_event(PHCompositeNode* topNode) {
  const int event = m_eventCounter;
  ++m_eventCounter;

  if (event < m_eventMin) {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  if (m_eventMax >= 0 && event > m_eventMax) {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  b_event = event;
  b_n_truth_particles = 0;
  b_n_truth_pi0 = 0;
  b_n_truth_gamma = 0;
  b_n_clusters = 0;
  b_n_towers_above_threshold = 0;
  b_n_cemc_hits = 0;

  PHG4TruthInfoContainer* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_truthNode);
  RawTowerGeomContainer* geom = findNode::getClass<RawTowerGeomContainer>(topNode, m_towerGeomNode);

  if (!truth) {
    if (m_missingTruthWarnings < 5) { /* print out warning just five times */
      std::cout << "Pi0EventDisplayDump::process_event - missing truth node " << m_truthNode << std::endl;
    }
    ++m_missingTruthWarnings;
  }
  if (!geom) {
    if (m_missingGeomWarnings < 5) { /* print out warning just five times */
      std::cout << "Pi0EventDisplayDump::process_event - missing CEMC geometry node " << m_towerGeomNode << std::endl;
    }
    ++m_missingGeomWarnings;
  }

  std::vector<TruthProjection> gamma_projections;
  if (truth) {
    fill_truth(truth, geom, gamma_projections);
  }

  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, m_towerNode);
  std::map<unsigned int, TowerLookup> tower_lookup;
  if (towers) {
    fill_towers(towers, geom, tower_lookup);
  } else {
    if (m_missingTowerWarnings < 5) { /* print out warning just five times */
      std::cout << "Pi0EventDisplayDump::process_event - missing tower node " << m_towerNode << std::endl;
    }
    ++m_missingTowerWarnings;
  }

  fill_clusters(topNode, gamma_projections, tower_lookup);

  if (!truth && b_n_clusters == 0) {
    std::cout << "Pi0EventDisplayDump::process_event - warning: both truth and useful cluster information are missing for event " << event << std::endl;
  }

  if (m_writeHits) {
    fill_hits(topNode, truth);
  }

  if (m_eventsTree) {
    m_eventsTree->Fill();
  }

  if (m_verbosity > 0 && event < m_eventMin + 5) {
    std::cout << "Pi0EventDisplayDump::process_event - event " << event
              << " truth pi0/gamma " << b_n_truth_pi0 << "/" << b_n_truth_gamma
              << ", clusters " << b_n_clusters
              << ", towers " << b_n_towers_above_threshold
              << ", hits " << b_n_cemc_hits << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0EventDisplayDump::End(PHCompositeNode* /*topNode*/) {
  close_output_file();
  std::cout << "Pi0EventDisplayDump::End - processed " << m_eventCounter << " events" << std::endl;
  std::cout << "Pi0EventDisplayDump::End - wrote " << m_outputFileName << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}




























void Pi0EventDisplayDump::fill_truth(PHG4TruthInfoContainer* truth, RawTowerGeomContainer* geom, std::vector<TruthProjection>& gamma_projections) {
  const double display_radius = get_display_radius(geom);
  PHG4TruthInfoContainer::ConstRange range = truth->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter) {
    PHG4Particle* particle = iter->second;
    if (!particle) {
      continue;
    }

    fill_truth_particle(truth, particle);
    fill_truth_segment(truth, particle, display_radius);

    if (particle->get_pid() == 22 && has_direct_parent_pid(truth, particle, 111)) {
      TruthProjection projection;
      projection.track_id = particle->get_track_id();
      double vx = 0.0;
      double vy = 0.0;
      double vz = 0.0;
      if (get_vertex(truth, particle->get_vtx_id(), vx, vy, vz) &&
          project_to_radius(vx, vy, vz, particle->get_px(), particle->get_py(), particle->get_pz(), display_radius, projection.x, projection.y, projection.z))
      {
        projection.eta = eta_from_xyz(projection.x, projection.y, projection.z);
        projection.phi = phi_from_xy(projection.x, projection.y);
      }
      gamma_projections.push_back(projection);
    }
  }
}

void Pi0EventDisplayDump::fill_truth_particle(PHG4TruthInfoContainer* truth, PHG4Particle* particle) {
  b_track_id = particle->get_track_id();
  b_pid = particle->get_pid();
  b_parent_id = particle->get_parent_id();
  b_primary_id = particle->get_primary_id();
  b_vtx_id = particle->get_vtx_id();
  b_px = finite_or_invalid(particle->get_px());
  b_py = finite_or_invalid(particle->get_py());
  b_pz = finite_or_invalid(particle->get_pz());
  b_e = finite_or_invalid(particle->get_e());

  if (!get_vertex(truth, b_vtx_id, b_vx, b_vy, b_vz)) {
    b_vx = kInvalidDouble;
    b_vy = kInvalidDouble;
    b_vz = kInvalidDouble;
  }

  const double pt = std::sqrt(particle->get_px() * particle->get_px() + particle->get_py() * particle->get_py());
  const double p = std::sqrt(pt * pt + particle->get_pz() * particle->get_pz());
  b_pt = finite_or_invalid(pt);
  b_p = finite_or_invalid(p);
  b_eta = p > 0.0 ? finite_or_invalid(std::asinh(particle->get_pz() / pt)) : kInvalidDouble;
  if (pt <= 0.0)
  {
    b_eta = kInvalidDouble;
  }
  b_phi = phi_from_xy(particle->get_px(), particle->get_py());

  b_is_primary = truth->is_primary(particle) ? 1 : 0;
  b_is_pi0 = b_pid == 111 ? 1 : 0;
  b_is_gamma = b_pid == 22 ? 1 : 0;
  b_is_pi0_daughter = has_direct_parent_pid(truth, particle, 111) ? 1 : 0;
  b_ancestor_pi0 = find_ancestor(truth, b_track_id, 111, nullptr);
  b_ancestor_gamma = find_ancestor(truth, b_track_id, 22, &b_generation);

  if (b_pid == 111)
  {
    ++b_n_truth_pi0;
  }
  if (b_pid == 22)
  {
    ++b_n_truth_gamma;
  }
  ++b_n_truth_particles;

  if (m_truthParticlesTree)
  {
    m_truthParticlesTree->Fill();
  }

  if (m_verbosity > 1 && (b_is_pi0 || b_is_pi0_daughter || b_is_gamma))
  {
    std::cout << "  truth track=" << b_track_id << " pid=" << b_pid << " parent=" << b_parent_id
              << " E=" << b_e << " eta=" << b_eta << " phi=" << b_phi << std::endl;
  }
}

void Pi0EventDisplayDump::fill_truth_segment(PHG4TruthInfoContainer* truth, PHG4Particle* particle, double /*display_radius*/)
{
  // b_event = b_event; // これは何
  b_track_id = particle->get_track_id();
  b_pid = particle->get_pid();
  b_parent_id = particle->get_parent_id();
  b_ancestor_pi0 = find_ancestor(truth, b_track_id, 111, nullptr);
  b_ancestor_gamma = find_ancestor(truth, b_track_id, 22, nullptr);

  if (!get_vertex(truth, particle->get_vtx_id(), b_x0, b_y0, b_z0)) {
    return;
  }

  bool have_end = false;
  if (particle->get_pid() == 111) {
    have_end = find_first_daughter_vertex(truth, particle->get_track_id(), b_x1, b_y1, b_z1);
    b_segment_type = 1;
  } else if (particle->get_pid() == 22 && has_direct_parent_pid(truth, particle, 111)) {
    b_segment_type = 2;
  } else {
    b_segment_type = 3;
  }

  if (!have_end) {
    b_x1 = b_x0;
    b_y1 = b_y0;
    b_z1 = b_z0;
  }

  b_r0 = radius(b_x0, b_y0);
  b_r1 = radius(b_x1, b_y1);
  b_eta0 = eta_from_xyz(b_x0, b_y0, b_z0);
  b_phi0 = phi_from_xy(b_x0, b_y0);
  b_eta1 = eta_from_xyz(b_x1, b_y1, b_z1);
  b_phi1 = phi_from_xy(b_x1, b_y1);

  if (m_truthSegmentsTree) {
    m_truthSegmentsTree->Fill();
  }
}

void Pi0EventDisplayDump::fill_clusters(PHCompositeNode* topNode, const std::vector<TruthProjection>& gamma_projections, const std::map<unsigned int, TowerLookup>& tower_lookup)
{
  RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, m_clusterNode);
  if (!clusters)
  {
    if (m_missingClusterWarnings < 5)
    {
      std::cout << "Pi0EventDisplayDump::process_event - missing cluster node " << m_clusterNode << std::endl;
    }
    ++m_missingClusterWarnings;
    return;
  }

  RawClusterContainer::ConstRange range = clusters->getClusters();
  for (RawClusterContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    const RawCluster* cluster = iter->second;
    if (!cluster)
    {
      continue;
    }

    const double energy = cluster->get_energy();
    if (energy < m_clusterEnergyMin)
    {
      continue;
    }

    b_cluster_node = m_clusterNode;
    b_cluster_id = cluster->get_id();
    b_energy = finite_or_invalid(energy);
    b_ecore = finite_or_invalid(cluster->get_ecore());
    b_chi2 = finite_or_invalid(cluster->get_chi2());
    b_prob = finite_or_invalid(cluster->get_prob());
    b_merged_cluster_prob = finite_or_invalid(cluster->get_merged_cluster_prob());
    b_x = finite_or_invalid(cluster->get_x());
    b_y = finite_or_invalid(cluster->get_y());
    b_z = finite_or_invalid(cluster->get_z());
    b_r = finite_or_invalid(cluster->get_r());
    if (b_r == kInvalidDouble)
    {
      b_r = radius(cluster->get_x(), cluster->get_y());
    }
    b_phi = phi_from_xy(cluster->get_x(), cluster->get_y());
    b_eta = b_r > 0.0 ? finite_or_invalid(std::asinh(cluster->get_z() / b_r)) : kInvalidDouble;
    b_ntowers = static_cast<int>(cluster->getNTowers());

    b_lead_tower_key = 0;
    b_lead_tower_ieta = kInvalidInt;
    b_lead_tower_iphi = kInvalidInt;
    b_lead_tower_energy = kInvalidDouble;

    double lead_value = -std::numeric_limits<double>::max();
    RawCluster::TowerConstRange tower_range = cluster->get_towers();
    for (RawCluster::TowerConstIterator tower_iter = tower_range.first; tower_iter != tower_range.second; ++tower_iter)
    {
      const unsigned int raw_key = tower_iter->first;
      const double contribution = tower_iter->second;
      if (contribution > lead_value)
      {
        lead_value = contribution;
        b_lead_tower_key = raw_key;
        b_lead_tower_ieta = static_cast<int>(RawTowerDefs::decode_index1(raw_key));
        b_lead_tower_iphi = static_cast<int>(RawTowerDefs::decode_index2(raw_key));
        std::map<unsigned int, TowerLookup>::const_iterator lookup = tower_lookup.find(raw_key);
        b_lead_tower_energy = lookup != tower_lookup.end() && lookup->second.tower ? finite_or_invalid(lookup->second.tower->get_energy()) : kInvalidDouble;
      }
    }

    b_nearest_truth_gamma_track_id = kInvalidInt;
    b_nearest_truth_gamma_delta_eta = kInvalidDouble;
    b_nearest_truth_gamma_delta_phi = kInvalidDouble;
    b_nearest_truth_gamma_delta_r = kInvalidDouble;
    b_nearest_truth_gamma_delta_tower = kInvalidDouble;
    double best_delta_r = std::numeric_limits<double>::max();
    for (const TruthProjection& gamma : gamma_projections)
    {
      if (gamma.track_id == kInvalidInt || gamma.eta == kInvalidDouble || gamma.phi == kInvalidDouble || b_eta == kInvalidDouble || b_phi == kInvalidDouble)
      {
        continue;
      }
      const double deta = b_eta - gamma.eta;
      const double dphi = delta_phi(b_phi, gamma.phi);
      const double dr = std::sqrt(deta * deta + dphi * dphi);
      if (dr < best_delta_r)
      {
        best_delta_r = dr;
        b_nearest_truth_gamma_track_id = gamma.track_id;
        b_nearest_truth_gamma_delta_eta = deta;
        b_nearest_truth_gamma_delta_phi = dphi;
        b_nearest_truth_gamma_delta_r = dr;
      }
    }

    if (m_cemcClustersTree)
    {
      m_cemcClustersTree->Fill();
    }
    ++b_n_clusters;

    if (m_verbosity > 1)
    {
      std::cout << "  cluster id=" << b_cluster_id << " E=" << b_energy << " eta=" << b_eta << " phi=" << b_phi
                << " ntower=" << b_ntowers << std::endl;
    }

    tower_range = cluster->get_towers();
    for (RawCluster::TowerConstIterator tower_iter = tower_range.first; tower_iter != tower_range.second; ++tower_iter)
    {
      b_cluster_node = m_clusterNode;
      b_cluster_id = cluster->get_id();
      b_tower_key = tower_iter->first;
      b_ieta = static_cast<int>(RawTowerDefs::decode_index1(b_tower_key));
      b_iphi = static_cast<int>(RawTowerDefs::decode_index2(b_tower_key));
      b_eta = kInvalidDouble;
      b_phi = kInvalidDouble;
      b_cluster_tower_value = finite_or_invalid(tower_iter->second);
      b_tower_energy = kInvalidDouble;

      std::map<unsigned int, TowerLookup>::const_iterator lookup = tower_lookup.find(b_tower_key);
      if (lookup != tower_lookup.end())
      {
        b_tower_energy = lookup->second.tower ? finite_or_invalid(lookup->second.tower->get_energy()) : kInvalidDouble;
      }

      if (m_cemcClusterTowersTree)
      {
        m_cemcClusterTowersTree->Fill();
      }
    }
  }
}

void Pi0EventDisplayDump::fill_towers(TowerInfoContainer* towers, RawTowerGeomContainer* geom, std::map<unsigned int, TowerLookup>& tower_lookup)
{
  const std::size_t n_towers = towers->size();
  for (std::size_t channel = 0; channel < n_towers; ++channel)
  {
    TowerInfo* tower = towers->get_tower_at_channel(static_cast<int>(channel));
    if (!tower)
    {
      continue;
    }

    const unsigned int tower_info_key = towers->encode_key(static_cast<unsigned int>(channel));
    const int ieta = static_cast<int>(towers->getTowerEtaBin(tower_info_key));
    const int iphi = static_cast<int>(towers->getTowerPhiBin(tower_info_key));
    const unsigned int raw_key = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));

    TowerLookup lookup;
    lookup.tower = tower;
    lookup.raw_key = raw_key;
    lookup.ieta = ieta;
    lookup.iphi = iphi;
    lookup.channel = static_cast<int>(channel);
    tower_lookup[raw_key] = lookup;

    const double energy = tower->get_energy();
    if (energy < m_towerEnergyMin)
    {
      continue;
    }

    b_tower_node = m_towerNode;
    b_channel = static_cast<int>(channel);
    b_tower_key = raw_key;
    b_ieta = ieta;
    b_iphi = iphi;
    b_eta = kInvalidDouble;
    b_phi = kInvalidDouble;
    b_x = kInvalidDouble;
    b_y = kInvalidDouble;
    b_z = kInvalidDouble;
    b_energy = finite_or_invalid(energy);
    b_time = finite_or_invalid(tower->get_time());
    b_is_good = tower->get_isGood() ? 1 : 0;
    b_status = tower->get_status();

    RawTowerGeom* tower_geom = get_tower_geom(geom, raw_key);
    if (tower_geom)
    {
      b_x = finite_or_invalid(tower_geom->get_center_x());
      b_y = finite_or_invalid(tower_geom->get_center_y());
      b_z = finite_or_invalid(tower_geom->get_center_z());
      b_eta = eta_from_xyz(b_x, b_y, b_z);
      b_phi = phi_from_xy(b_x, b_y);
    }

    if (m_cemcTowersTree)
    {
      m_cemcTowersTree->Fill();
    }
    ++b_n_towers_above_threshold;
  }
}

void Pi0EventDisplayDump::fill_hits(PHCompositeNode* topNode, PHG4TruthInfoContainer* truth)
{
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, m_cemcHitNode);
  if (!hits)
  {
    if (m_missingHitWarnings < 5)
    {
      std::cout << "Pi0EventDisplayDump::process_event - missing optional hit node " << m_cemcHitNode << std::endl;
    }
    ++m_missingHitWarnings;
    return;
  }

  PHG4HitContainer::ConstRange range = hits->getHits();
  for (PHG4HitContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    PHG4Hit* hit = iter->second;
    if (!hit)
    {
      continue;
    }

    b_trkid = hit->get_trkid();
    b_ancestor_pi0 = truth ? find_ancestor(truth, b_trkid, 111, nullptr) : kInvalidInt;
    b_ancestor_gamma = truth ? find_ancestor(truth, b_trkid, 22, nullptr) : kInvalidInt;
    b_x0 = finite_or_invalid(hit->get_x(0));
    b_y0 = finite_or_invalid(hit->get_y(0));
    b_z0 = finite_or_invalid(hit->get_z(0));
    b_x1 = finite_or_invalid(hit->get_x(1));
    b_y1 = finite_or_invalid(hit->get_y(1));
    b_z1 = finite_or_invalid(hit->get_z(1));
    b_r0 = radius(b_x0, b_y0);
    b_r1 = radius(b_x1, b_y1);
    b_edep = finite_or_invalid(hit->get_edep());
    b_eion = finite_or_invalid(hit->get_eion());
    b_light_yield = finite_or_invalid(hit->get_light_yield());

    if (m_cemcHitsTree)
    {
      m_cemcHitsTree->Fill();
    }
    ++b_n_cemc_hits;
  }
}

int Pi0EventDisplayDump::find_ancestor(PHG4TruthInfoContainer* truth, int track_id, int pid, int* generation) const
{
  if (generation)
  {
    *generation = kInvalidInt;
  }
  if (!truth)
  {
    return kInvalidInt;
  }

  int current_id = track_id;
  int depth = 0;
  std::set<int> visited;
  while (current_id != 0 && current_id != kInvalidInt && visited.insert(current_id).second)
  {
    PHG4Particle* particle = truth->GetParticle(current_id);
    if (!particle)
    {
      break;
    }
    if (particle->get_pid() == pid)
    {
      if (generation)
      {
        *generation = depth;
      }
      return particle->get_track_id();
    }
    current_id = particle->get_parent_id();
    ++depth;
  }

  return kInvalidInt;
}

bool Pi0EventDisplayDump::has_direct_parent_pid(PHG4TruthInfoContainer* truth, const PHG4Particle* particle, int pid) const
{
  if (!truth || !particle)
  {
    return false;
  }

  PHG4Particle* parent = truth->GetParticle(particle->get_parent_id());
  return parent && parent->get_pid() == pid;
}

bool Pi0EventDisplayDump::get_vertex(PHG4TruthInfoContainer* truth, int vtx_id, double& x, double& y, double& z) const {
  if (!truth) {
    return false;
  }
  PHG4VtxPoint* vertex = truth->GetVtx(vtx_id);
  if (!vertex) {
    return false;
  }
  x = vertex->get_x();
  y = vertex->get_y();
  z = vertex->get_z();
  return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

bool Pi0EventDisplayDump::project_to_radius(double x0, double y0, double z0, double px, double py, double pz, double radius_target, double& x1, double& y1, double& z1) const
{
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (!std::isfinite(p) || p <= 0.0 || radius_target <= 0.0)
  {
    return false;
  }

  const double ux = px / p;
  const double uy = py / p;
  const double uz = pz / p;
  const double a = ux * ux + uy * uy;
  const double b = 2.0 * (x0 * ux + y0 * uy);
  const double c = x0 * x0 + y0 * y0 - radius_target * radius_target;
  if (a <= std::numeric_limits<double>::epsilon())
  {
    return false;
  }

  const double discriminant = b * b - 4.0 * a * c;
  if (discriminant < 0.0)
  {
    return false;
  }

  const double root = std::sqrt(discriminant);
  const double s1 = (-b + root) / (2.0 * a);
  const double s2 = (-b - root) / (2.0 * a);
  double s = std::numeric_limits<double>::max();
  if (s1 > 0.0)
  {
    s = std::min(s, s1);
  }
  if (s2 > 0.0)
  {
    s = std::min(s, s2);
  }
  if (s == std::numeric_limits<double>::max())
  {
    return false;
  }

  x1 = x0 + s * ux;
  y1 = y0 + s * uy;
  z1 = z0 + s * uz;
  return true;
}

bool Pi0EventDisplayDump::find_first_daughter_vertex(PHG4TruthInfoContainer* truth, int track_id, double& x, double& y, double& z) const
{
  if (!truth)
  {
    return false;
  }

  PHG4TruthInfoContainer::ConstRange range = truth->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    PHG4Particle* particle = iter->second;
    if (!particle || particle->get_parent_id() != track_id)
    {
      continue;
    }
    if (get_vertex(truth, particle->get_vtx_id(), x, y, z))
    {
      return true;
    }
  }

  return false;
}

double Pi0EventDisplayDump::get_display_radius(RawTowerGeomContainer* geom) const {
  if (!geom) {
    return kDefaultDisplayRadius;
  }

  const RawTowerGeom* sample = nullptr;
  for (int ieta = 0; ieta < 96 && !sample; ++ieta) {
    for (int iphi = 0; iphi < 256 && !sample; ++iphi) {
      const unsigned int key = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));
      sample = geom->get_tower_geometry(key);
    }
  }

  if (!sample) {
    return kDefaultDisplayRadius;
  }

  const double r = radius(sample->get_center_x(), sample->get_center_y());
  return std::isfinite(r) && r > 0.0 ? r : kDefaultDisplayRadius;
}

RawTowerGeom* Pi0EventDisplayDump::get_tower_geom(RawTowerGeomContainer* geom, unsigned int raw_key) const
{
  return geom ? geom->get_tower_geometry(raw_key) : nullptr;
}

double Pi0EventDisplayDump::radius(double x, double y)
{
  if (!std::isfinite(x) || !std::isfinite(y) || x == kInvalidDouble || y == kInvalidDouble) {
    return kInvalidDouble;
  }
  return std::sqrt(x * x + y * y);
}

double Pi0EventDisplayDump::eta_from_xyz(double x, double y, double z) {
  const double r = radius(x, y);
  if (r <= 0.0 || r == kInvalidDouble || !std::isfinite(z) || z == kInvalidDouble) {
    return kInvalidDouble;
  }
  return std::asinh(z / r);
}

double Pi0EventDisplayDump::phi_from_xy(double x, double y) {
  if (!std::isfinite(x) || !std::isfinite(y) || x == kInvalidDouble || y == kInvalidDouble)
  {
    return kInvalidDouble;
  }
  return std::atan2(y, x);
}

double Pi0EventDisplayDump::delta_phi(double a, double b) {
  double dphi = a - b;
  while (dphi > M_PI)
  {
    dphi -= 2.0 * M_PI;
  }
  while (dphi <= -M_PI)
  {
    dphi += 2.0 * M_PI;
  }
  return dphi;
}

double Pi0EventDisplayDump::finite_or_invalid(double value) {
  return std::isfinite(value) ? value : kInvalidDouble;
}

int Pi0EventDisplayDump::finite_or_invalid_int(double value) {
  return std::isfinite(value) ? static_cast<int>(value) : kInvalidInt;
}

void Pi0EventDisplayDump::create_output_directory() const {
  const std::string::size_type slash_position = m_outputFileName.find_last_of('/');
  if (slash_position == std::string::npos) {
    return;
  }

  const std::string directory = m_outputFileName.substr(0, slash_position);
  if (!directory.empty()) {
    gSystem->mkdir(directory.c_str(), true);
  }
}

void Pi0EventDisplayDump::create_trees() {
  m_eventsTree = new TTree("events", "Event display event summary");
  m_eventsTree->Branch("event", &b_event);
  m_eventsTree->Branch("n_truth_particles", &b_n_truth_particles);
  m_eventsTree->Branch("n_truth_pi0", &b_n_truth_pi0);
  m_eventsTree->Branch("n_truth_gamma", &b_n_truth_gamma);
  m_eventsTree->Branch("n_clusters", &b_n_clusters);
  m_eventsTree->Branch("n_towers_above_threshold", &b_n_towers_above_threshold);
  m_eventsTree->Branch("n_cemc_hits", &b_n_cemc_hits);

  m_truthParticlesTree = new TTree("truth_particles", "Truth particles");
  m_truthParticlesTree->Branch("event", &b_event);
  m_truthParticlesTree->Branch("track_id", &b_track_id);
  m_truthParticlesTree->Branch("pid", &b_pid);
  m_truthParticlesTree->Branch("parent_id", &b_parent_id);
  m_truthParticlesTree->Branch("primary_id", &b_primary_id);
  m_truthParticlesTree->Branch("vtx_id", &b_vtx_id);
  m_truthParticlesTree->Branch("px", &b_px);
  m_truthParticlesTree->Branch("py", &b_py);
  m_truthParticlesTree->Branch("pz", &b_pz);
  m_truthParticlesTree->Branch("e", &b_e);
  m_truthParticlesTree->Branch("vx", &b_vx);
  m_truthParticlesTree->Branch("vy", &b_vy);
  m_truthParticlesTree->Branch("vz", &b_vz);
  m_truthParticlesTree->Branch("pt", &b_pt);
  m_truthParticlesTree->Branch("p", &b_p);
  m_truthParticlesTree->Branch("eta", &b_eta);
  m_truthParticlesTree->Branch("phi", &b_phi);
  m_truthParticlesTree->Branch("is_primary", &b_is_primary);
  m_truthParticlesTree->Branch("is_pi0", &b_is_pi0);
  m_truthParticlesTree->Branch("is_gamma", &b_is_gamma);
  m_truthParticlesTree->Branch("is_pi0_daughter", &b_is_pi0_daughter);
  m_truthParticlesTree->Branch("ancestor_pi0", &b_ancestor_pi0);
  m_truthParticlesTree->Branch("ancestor_gamma", &b_ancestor_gamma);
  m_truthParticlesTree->Branch("generation", &b_generation);

  m_truthSegmentsTree = new TTree("truth_segments", "Truth display line segments");
  m_truthSegmentsTree->Branch("event", &b_event);
  m_truthSegmentsTree->Branch("track_id", &b_track_id);
  m_truthSegmentsTree->Branch("pid", &b_pid);
  m_truthSegmentsTree->Branch("parent_id", &b_parent_id);
  m_truthSegmentsTree->Branch("ancestor_pi0", &b_ancestor_pi0);
  m_truthSegmentsTree->Branch("ancestor_gamma", &b_ancestor_gamma);
  m_truthSegmentsTree->Branch("x0", &b_x0);
  m_truthSegmentsTree->Branch("y0", &b_y0);
  m_truthSegmentsTree->Branch("z0", &b_z0);
  m_truthSegmentsTree->Branch("x1", &b_x1);
  m_truthSegmentsTree->Branch("y1", &b_y1);
  m_truthSegmentsTree->Branch("z1", &b_z1);
  m_truthSegmentsTree->Branch("r0", &b_r0);
  m_truthSegmentsTree->Branch("r1", &b_r1);
  m_truthSegmentsTree->Branch("eta0", &b_eta0);
  m_truthSegmentsTree->Branch("phi0", &b_phi0);
  m_truthSegmentsTree->Branch("eta1", &b_eta1);
  m_truthSegmentsTree->Branch("phi1", &b_phi1);
  m_truthSegmentsTree->Branch("segment_type", &b_segment_type);

  m_cemcClustersTree = new TTree("cemc_clusters", "CEMC clusters");
  m_cemcClustersTree->Branch("event", &b_event);
  m_cemcClustersTree->Branch("cluster_node", &b_cluster_node);
  m_cemcClustersTree->Branch("cluster_id", &b_cluster_id);
  m_cemcClustersTree->Branch("energy", &b_energy);
  m_cemcClustersTree->Branch("ecore", &b_ecore);
  m_cemcClustersTree->Branch("chi2", &b_chi2);
  m_cemcClustersTree->Branch("prob", &b_prob);
  m_cemcClustersTree->Branch("merged_cluster_prob", &b_merged_cluster_prob);
  m_cemcClustersTree->Branch("x", &b_x);
  m_cemcClustersTree->Branch("y", &b_y);
  m_cemcClustersTree->Branch("z", &b_z);
  m_cemcClustersTree->Branch("r", &b_r);
  m_cemcClustersTree->Branch("phi", &b_phi);
  m_cemcClustersTree->Branch("eta_vtx0", &b_eta);
  m_cemcClustersTree->Branch("ntowers", &b_ntowers);
  m_cemcClustersTree->Branch("lead_tower_key", &b_lead_tower_key);
  m_cemcClustersTree->Branch("lead_tower_ieta", &b_lead_tower_ieta);
  m_cemcClustersTree->Branch("lead_tower_iphi", &b_lead_tower_iphi);
  m_cemcClustersTree->Branch("lead_tower_energy", &b_lead_tower_energy);
  m_cemcClustersTree->Branch("nearest_truth_gamma_track_id", &b_nearest_truth_gamma_track_id);
  m_cemcClustersTree->Branch("nearest_truth_gamma_delta_eta", &b_nearest_truth_gamma_delta_eta);
  m_cemcClustersTree->Branch("nearest_truth_gamma_delta_phi", &b_nearest_truth_gamma_delta_phi);
  m_cemcClustersTree->Branch("nearest_truth_gamma_delta_r", &b_nearest_truth_gamma_delta_r);
  m_cemcClustersTree->Branch("nearest_truth_gamma_delta_tower", &b_nearest_truth_gamma_delta_tower);

  m_cemcClusterTowersTree = new TTree("cemc_cluster_towers", "CEMC cluster tower contributions");
  m_cemcClusterTowersTree->Branch("event", &b_event);
  m_cemcClusterTowersTree->Branch("cluster_node", &b_cluster_node);
  m_cemcClusterTowersTree->Branch("cluster_id", &b_cluster_id);
  m_cemcClusterTowersTree->Branch("tower_key", &b_tower_key);
  m_cemcClusterTowersTree->Branch("ieta", &b_ieta);
  m_cemcClusterTowersTree->Branch("iphi", &b_iphi);
  m_cemcClusterTowersTree->Branch("eta", &b_eta);
  m_cemcClusterTowersTree->Branch("phi", &b_phi);
  m_cemcClusterTowersTree->Branch("cluster_tower_value", &b_cluster_tower_value);
  m_cemcClusterTowersTree->Branch("tower_energy", &b_tower_energy);

  m_cemcTowersTree = new TTree("cemc_towers", "CEMC towers");
  m_cemcTowersTree->Branch("event", &b_event);
  m_cemcTowersTree->Branch("tower_node", &b_tower_node);
  m_cemcTowersTree->Branch("channel", &b_channel);
  m_cemcTowersTree->Branch("key", &b_tower_key);
  m_cemcTowersTree->Branch("ieta", &b_ieta);
  m_cemcTowersTree->Branch("iphi", &b_iphi);
  m_cemcTowersTree->Branch("eta", &b_eta);
  m_cemcTowersTree->Branch("phi", &b_phi);
  m_cemcTowersTree->Branch("x", &b_x);
  m_cemcTowersTree->Branch("y", &b_y);
  m_cemcTowersTree->Branch("z", &b_z);
  m_cemcTowersTree->Branch("energy", &b_energy);
  m_cemcTowersTree->Branch("time", &b_time);
  m_cemcTowersTree->Branch("is_good", &b_is_good);
  m_cemcTowersTree->Branch("status", &b_status);

  m_cemcHitsTree = new TTree("cemc_hits", "CEMC G4 hits");
  m_cemcHitsTree->Branch("event", &b_event);
  m_cemcHitsTree->Branch("trkid", &b_trkid);
  m_cemcHitsTree->Branch("ancestor_pi0", &b_ancestor_pi0);
  m_cemcHitsTree->Branch("ancestor_gamma", &b_ancestor_gamma);
  m_cemcHitsTree->Branch("x0", &b_x0);
  m_cemcHitsTree->Branch("y0", &b_y0);
  m_cemcHitsTree->Branch("z0", &b_z0);
  m_cemcHitsTree->Branch("x1", &b_x1);
  m_cemcHitsTree->Branch("y1", &b_y1);
  m_cemcHitsTree->Branch("z1", &b_z1);
  m_cemcHitsTree->Branch("r0", &b_r0);
  m_cemcHitsTree->Branch("r1", &b_r1);
  m_cemcHitsTree->Branch("edep", &b_edep);
  m_cemcHitsTree->Branch("eion", &b_eion);
  m_cemcHitsTree->Branch("light_yield", &b_light_yield);
}

void Pi0EventDisplayDump::close_output_file()
{
  if (!m_outputFile)
  {
    return;
  }

  m_outputFile->cd();
  if (m_eventsTree)
  {
    m_eventsTree->Write();
  }
  if (m_truthParticlesTree)
  {
    m_truthParticlesTree->Write();
  }
  if (m_truthSegmentsTree)
  {
    m_truthSegmentsTree->Write();
  }
  if (m_cemcClustersTree)
  {
    m_cemcClustersTree->Write();
  }
  if (m_cemcClusterTowersTree)
  {
    m_cemcClusterTowersTree->Write();
  }
  if (m_cemcTowersTree)
  {
    m_cemcTowersTree->Write();
  }
  if (m_cemcHitsTree)
  {
    m_cemcHitsTree->Write();
  }
  m_outputFile->Close();
  delete m_outputFile;
  m_outputFile = nullptr;
}
