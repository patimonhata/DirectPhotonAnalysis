// Read-only, bounded DST diagnostic using the frozen production evaluator.
// Outputs are CSVs in a NEW directory; no production setting or library is changed.
#include <Pi0AnchorTopologyEvaluator.h>
#include <Pi0ClusterTruthMatcher.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>
#include <g4calo/RawTowerBuilder.h>
#include <g4detectors/PHG4FullProjSpacalCellReco.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/getClass.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libPhotonAnalysisTree.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4calo.so)

namespace common_matching_diagnostic
{
using namespace photon_tree;
struct Below { const RawCluster* cluster; double eta, phi; std::vector<Pi0ClusterTruthMatch> matches; };
class Audit : public SubsysReco
{
 public:
  explicit Audit(const std::string& output) : SubsysReco("CommonMatchingAudit"), candidates(output + "/candidates.csv"),
      failures(output + "/first_failures.csv"), anchors(output + "/anchors.csv"), events(output + "/events.csv")
  {
    if (!candidates || !failures || !anchors || !events) throw std::runtime_error("Cannot create diagnostic CSVs");
    candidates << "event_in_file,hepmc_event_number,pdg,candidate_index,parent_source,parent_id,photon0_track,photon1_track,global_usable,above_count,above_failures,below_count,below_failures,projection_ok,interaction_ok,anchor_count,all_anchors_usable,representatives_found,representatives_usable\n";
    failures << "event_in_file,hepmc_event_number,pdg,candidate_index,stage,cluster_id,cluster_energy,cluster_eta,cluster_phi,min_photon_delta_r,is_anchor,is_representative,status,failure,ancestry_complete,coverage,tower_count,matched_towers,failure_ieta,failure_iphi,gamma0_edep,gamma1_edep,total_edep\n";
    anchors << "event_in_file,hepmc_event_number,pdg,candidate_index,cluster_id,et,main_fraction,topology,applicability,topology_reason,anchor_match_usable,anchor_match_status,anchor_match_coverage\n";
    events << "event_in_file,hepmc_event_number,pdg,status,candidates,anchors,above_clusters,below_clusters\n";
    for (auto* stream : {&candidates, &failures, &anchors, &events}) *stream << std::setprecision(12);
  }
  int process_event(PHCompositeNode* top) override
  {
    const unsigned int index = event_index++;
    auto* event_map = findNode::getClass<PHHepMCGenEventMap>(top, "PHHepMCGenEventMap");
    const auto* sub = event_map ? event_map->get(1) : nullptr;
    const int hep_event = sub && sub->getEvent() ? sub->getEvent()->event_number() : -1;
    auto* clusters = findNode::getClass<RawClusterContainer>(top, "CLUSTERINFO_CEMC");
    auto* towers = findNode::getClass<TowerInfoContainer>(top, "TOWERINFO_CALIB_CEMC");
    auto* raw = findNode::getClass<RawTowerContainer>(top, "TOWER_SIM_CEMC");
    auto* cells = findNode::getClass<PHG4CellContainer>(top, "G4CELL_CEMC");
    auto* hits = findNode::getClass<PHG4HitContainer>(top, "G4HIT_CEMC");
    auto* truth = findNode::getClass<PHG4TruthInfoContainer>(top, "G4TruthInfo");
    for (int pdg : {111, 221})
    {
      Pi0AnchorTopologyConfig config;
      config.common_composition = true; config.eta_pseudo_primary = true; config.parent_pdg_id = pdg;
      config.require_hepmc_two_photon_decay = pdg == 221;
      config.min_cluster_energy = std::nextafter(0.12, std::numeric_limits<double>::infinity());
      config.tagging_partner_min_cluster_energy = 0.5; config.missing_energy_min = 0.12; config.missing_energy_max = 0.5;
      config.tagging_pi0_mass_min = pdg == 111 ? 0.10 : 0.45; config.tagging_pi0_mass_max = pdg == 111 ? 0.20 : 0.65;
      config.max_abs_vertex_z = 60.; config.enable_missing_diagnostics = true;
      Pi0AnchorTopologyEvaluator evaluator; evaluator.configure(config);
      auto result = evaluator.evaluate(top);
      if (result.status != Pi0TopologyEventStatus::accepted)
      { events << index << ',' << hep_event << ',' << pdg << ',' << int(result.status) << ",0,0,0,0\n"; continue; }
      std::vector<std::size_t> indices;
      std::vector<std::array<int, 2>> photon_ids;
      for (std::size_t i = 0; i < result.candidates.size(); ++i) if (result.candidates[i].topology_evaluated)
      { indices.push_back(i); photon_ids.push_back(result.candidates[i].photon_track_ids); }
      std::vector<Below> below;
      const auto range = clusters->getClusters();
      for (auto it = range.first; it != range.second; ++it)
      {
        const auto* c = it->second;
        if (!c || !std::isfinite(c->get_energy()) || !(c->get_energy() > 0) || c->get_energy() >= config.min_cluster_energy ||
            !std::isfinite(c->get_x()) || !std::isfinite(c->get_y()) || !std::isfinite(c->get_z())) continue;
        const double dx = c->get_x() - result.collision_vertex[0], dy = c->get_y() - result.collision_vertex[1];
        const double dz = c->get_z() - result.collision_vertex[2], transverse = std::hypot(dx, dy), distance = std::hypot(transverse, dz);
        if (!(distance > 0) || !(transverse > 0)) continue;
        const double eta = std::asinh(dz / transverse), phi = std::atan2(dy, dx), et = c->get_energy() * transverse / distance;
        if (!std::isfinite(et) || !(et > 0) || !std::isfinite(eta) || !std::isfinite(phi)) continue;
        auto matches = Pi0ClusterTruthMatcher().match_many(c, towers, raw, cells, hits, truth, photon_ids, true);
        if (matches.size() != indices.size()) throw std::runtime_error("Diagnostic match size mismatch");
        for (auto& m : matches)
        {
          if (m.status == Pi0ClusterTruthMatchStatus::partial && m.cluster_member_energy_coverage >= config.min_direct_match_cluster_energy_coverage) m.usable = true;
          if (!m.ancestry_complete) m.usable = false;
        }
        below.push_back({c, eta, phi, std::move(matches)});
      }
      events << index << ',' << hep_event << ',' << pdg << ',' << int(result.status) << ',' << result.candidates.size() << ','
             << result.anchors.size() << ',' << result.clusters.size() << ',' << below.size() << '\n';
      for (std::size_t pos = 0; pos < indices.size(); ++pos)
      {
        const std::size_t ci = indices[pos]; const auto& c = result.candidates[ci];
        std::map<unsigned int, const Pi0TopologyAnchorRecord*> owned;
        bool anchors_usable = true; int representatives_found = 0; bool representatives_usable = true;
        for (const auto& a : result.anchors) if (a.candidate_index == ci)
        {
          const auto& cl = result.clusters[a.cluster_index]; const auto& m = c.cluster_matches[a.cluster_index];
          owned[cl.cluster_id] = &a; anchors_usable &= m.usable;
          anchors << index << ',' << hep_event << ',' << pdg << ',' << ci << ',' << cl.cluster_id << ',' << cl.et << ',' << a.main_fraction << ','
                  << int(a.topology) << ',' << int(cl.composition.applicability) << ',' << cl.composition.topology_reason << ','
                  << m.usable << ',' << int(m.status) << ',' << m.cluster_member_energy_coverage << '\n';
        }
        for (const auto& p : c.truth_partner_clusters)
        { representatives_found += p.found; representatives_usable &= p.found; }
        // A found representative is selected only from usable matches by the production evaluator.
        const auto failure = [&](const char* stage, unsigned int id, double energy, double eta, double phi, const Pi0ClusterTruthMatch& m) {
          double dr = std::numeric_limits<double>::infinity();
          for (int i = 0; i < 2; ++i)
          {
            const double ref_eta = c.photon_projection_valid[i] ? c.photon_projection_eta[i] : c.photon_eta[i];
            const double ref_phi = c.photon_projection_valid[i] ? c.photon_projection_phi[i] : c.photon_phi[i];
            const double dphi = std::atan2(std::sin(phi - ref_phi), std::cos(phi - ref_phi));
            dr = std::min(dr, std::hypot(eta - ref_eta, dphi));
          }
          const bool representative = (c.truth_partner_clusters[0].found && c.truth_partner_clusters[0].cluster_id == id) ||
                                      (c.truth_partner_clusters[1].found && c.truth_partner_clusters[1].cluster_id == id);
          failures << index << ',' << hep_event << ',' << pdg << ',' << ci << ',' << stage << ',' << id << ',' << energy << ',' << eta << ',' << phi << ','
                   << dr << ',' << owned.count(id) << ',' << representative << ',' << int(m.status) << ',' << pi0_cluster_truth_match_failure_name(m.failure) << ','
                   << m.ancestry_complete << ',' << m.cluster_member_energy_coverage << ',' << m.tower_count << ',' << m.matched_tower_count << ','
                   << m.failure_ieta << ',' << m.failure_iphi << ',' << m.gamma_edep[0] << ',' << m.gamma_edep[1] << ',' << m.total_edep << '\n';
        };
        int above_bad = 0, below_bad = 0;
        for (std::size_t i = 0; i < result.clusters.size(); ++i) if (!c.cluster_matches[i].usable)
        {
          const auto& cl = result.clusters[i];
          if (above_bad++ == 0) failure("above", cl.cluster_id, cl.energy, cl.eta, cl.phi, c.cluster_matches[i]);
        }
        for (const auto& cl : below) if (!cl.matches[pos].usable)
        { if (below_bad++ == 0) failure("below", cl.cluster->get_id(), cl.cluster->get_energy(), cl.eta, cl.phi, cl.matches[pos]); }
        if (c.matching_usable != (above_bad == 0 && below_bad == 0)) throw std::runtime_error("Failed to reproduce production matching gate");
        const bool projection = (c.photon_projection_valid[0] || c.photon_projection_no_intersection[0]) &&
                                (c.photon_projection_valid[1] || c.photon_projection_no_intersection[1]);
        candidates << index << ',' << hep_event << ',' << pdg << ',' << ci << ',' << c.parent_identity.source << ',' << c.parent_identity.id << ','
                   << c.photon_track_ids[0] << ',' << c.photon_track_ids[1] << ',' << c.matching_usable << ',' << result.clusters.size() << ',' << above_bad << ','
                   << below.size() << ',' << below_bad << ',' << projection << ',' << c.interaction_info_valid << ',' << owned.size() << ',' << anchors_usable << ','
                   << representatives_found << ',' << representatives_usable << '\n';
      }
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }
 private:
  unsigned int event_index = 0;
  std::ofstream candidates, failures, anchors, events;
};
}

int DiagnoseCommonMatching(const std::string manifest_path, const int manifest_row, const std::string output, const int n_events = 30)
{
  if (n_events <= 0 || manifest_row < 0 || !gSystem->AccessPathName(output.c_str())) throw std::runtime_error("Require positive event limit and new output directory");
  std::ifstream manifest(manifest_path); std::string suffix;
  for (int i = 0; i <= manifest_row; ++i) if (!std::getline(manifest, suffix)) throw std::runtime_error("Cannot read manifest row");
  if (suffix.empty() || suffix.find('/') != std::string::npos || suffix.find('\r') != std::string::npos) throw std::runtime_error("Invalid suffix");
  const char* calibration = std::getenv("CALIBRATIONROOT");
  if (!calibration) throw std::runtime_error("CALIBRATIONROOT missing");
  if (gSystem->mkdir(output.c_str(), true) != 0) throw std::runtime_error("Cannot create output directory");
  std::ofstream provenance(output + "/input.txt");
  provenance << "manifest=" << manifest_path << "\nrow=" << manifest_row << "\nsuffix=" << suffix << "\nrequested_events=" << n_events
             << "\nselection=all accepted vertex events; no stitching or kinematic selection in this diagnostic\n";
  auto* server = Fun4AllServer::instance(); server->Verbosity(0);
  for (const std::string& prefix : {"DST_CALO_CLUSTER", "DST_MBD_EPD", "DST_TRUTH_JET", "G4Hits"})
  {
    auto* input = new Fun4AllDstInputManager(prefix); input->AddFile(prefix + "_" + suffix); server->registerInputManager(input);
  }
  auto* cells = new PHG4FullProjSpacalCellReco("CEMCCYLCELLRECO"); cells->Detector("CEMC");
  cells->get_light_collection_model().load_data_file(std::string(calibration) + "/CEMC/LightCollection/Prototype3Module.xml", "data_grid_light_guide_efficiency", "data_grid_fiber_trans");
  server->registerSubsystem(cells);
  auto* towers = new RawTowerBuilder("CEMCTruthRawTowerBuilder"); towers->Detector("CEMC");
  towers->set_sim_tower_node_prefix("SIM"); towers->set_towerinfo(RawTowerBuilder::ProcessTowerType::kRawTowerOnly); server->registerSubsystem(towers);
  server->registerSubsystem(new common_matching_diagnostic::Audit(output));
  const int run_status = server->run(n_events), end_status = server->End(); delete server;
  if (run_status != Fun4AllReturnCodes::EVENT_OK || end_status != Fun4AllReturnCodes::EVENT_OK) throw std::runtime_error("DST diagnostic did not complete normally");
  return 0;
}
