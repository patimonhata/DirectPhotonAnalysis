#include "EtaCascadeTruth.h"
#include "Pi0ClusterTruthMatcher.h"
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4MCProcessDefs.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <algorithm>
#include <cmath>
#include <set>

namespace
{
int decay_mode(std::vector<int> pdgs)
{
  std::sort(pdgs.begin(), pdgs.end());
  if (pdgs == std::vector<int>{22, 22}) return 1;
  if (pdgs == std::vector<int>{111, 111, 111}) return 2;
  if (pdgs == std::vector<int>{-211, 111, 211}) return 3;
  return pdgs.empty() ? 0 : 4;
}
const HepMC::GenParticle* unique_parent(const HepMC::GenParticle* particle)
{
  const auto* vertex = particle ? particle->production_vertex() : nullptr;
  return vertex && vertex->particles_in_size() == 1 ? *vertex->particles_in_const_begin() : nullptr;
}
int hepmc_mode(const HepMC::GenParticle* eta)
{
  const auto* vertex = eta ? eta->end_vertex() : nullptr;
  std::vector<int> pdgs;
  if (vertex) for (auto it = vertex->particles_out_const_begin(); it != vertex->particles_out_const_end(); ++it) pdgs.push_back((*it)->pdg_id());
  return decay_mode(pdgs);
}
}
namespace photon_tree
{
EtaCascadeTruth::EtaCascadeTruth(PHG4TruthInfoContainer* truth, const PHHepMCGenEventMap* events) : truth_(truth), events_(events)
{
  if (!truth_) return;
  const auto range = truth_->GetParticleRange();
  for (auto it = range.first; it != range.second; ++it) if (it->second) children_[it->second->get_parent_id()].push_back(it->second);
}
int EtaCascadeTruth::g4_mode(int eta_track) const
{
  const auto found = children_.find(eta_track);
  if (found == children_.end()) return 0;
  std::vector<int> pdgs;
  const PHG4VtxPoint* first = nullptr;
  for (const auto* child : found->second)
  {
    const auto* vertex = truth_->GetVtx(child->get_vtx_id());
    // Require an explicit G4 decay, rather than reinterpreting a material interaction.
    if (!vertex || vertex->get_process() != PHG4MCProcess::kPDecay) return 0;
    if (first && (first->get_x() != vertex->get_x() || first->get_y() != vertex->get_y() ||
                  first->get_z() != vertex->get_z() || first->get_t() != vertex->get_t())) return 0;
    first = vertex;
    pdgs.push_back(child->get_pid());
  }
  return decay_mode(pdgs);
}
bool EtaCascadeTruth::is_eta_daughter_pi0(const PHG4Particle* particle) const
{
  if (!particle || particle->get_pid() != 111) return false;
  const auto* parent = truth_->GetParticle(particle->get_parent_id());
  if (!parent || parent->get_pid() != 221 || parent->get_parent_id() != 0) return false;
  const int mode = g4_mode(parent->get_track_id());
  return mode == 2 || mode == 3;
}
void EtaCascadeTruth::annotate(TruthContributor& contributor) const
{
  const auto* current = truth_ ? truth_->GetParticle(contributor.g4_track_id) : nullptr;
  const PHG4Particle* root = current;
  std::set<int> visited;
  while (current && visited.insert(current->get_track_id()).second)
  {
    if (current->get_pid() == 111 && contributor.pi0_parent_track_id == 0) contributor.pi0_parent_track_id = current->get_track_id();
    if (current->get_pid() == 221)
    {
      contributor.eta_parent_track_id = current->get_track_id();
      contributor.eta_decay_mode = g4_mode(current->get_track_id());
      // A primary eta has a HepMC identity; secondary barcodes need not be unique.
      if (current->get_parent_id() == 0)
      {
        const auto* subevent = events_ ? events_->get(contributor.embedding_id) : nullptr;
        const auto* event = subevent ? subevent->getEvent() : nullptr;
        const auto* eta = event ? event->barcode_to_particle(current->get_barcode()) : nullptr;
        if (eta && eta->pdg_id() == 221) contributor.eta_parent_barcode = eta->barcode();
      }
      return;
    }
    root = current;
    current = truth_->GetParticle(current->get_parent_id());
  }
  const auto* subevent = events_ ? events_->get(contributor.embedding_id) : nullptr;
  const auto* event = subevent ? subevent->getEvent() : nullptr;
  const auto* particle = event && root ? event->barcode_to_particle(root->get_barcode()) : nullptr;
  std::set<const HepMC::GenParticle*> seen;
  bool saw_pi0 = contributor.pi0_parent_track_id != 0;
  while (particle && seen.insert(particle).second)
  {
    const int pdg = particle->pdg_id();
    if (pdg == 111 && contributor.pi0_parent_barcode == 0) { contributor.pi0_parent_barcode = particle->barcode(); saw_pi0 = true; }
    if (pdg == 221)
    {
      contributor.eta_parent_barcode = particle->barcode();
      contributor.eta_decay_mode = hepmc_mode(particle);
      return;
    }
    // Follow only photon bookkeeping and the relevant pion/eta decay ancestry.
    if (pdg != 22 && pdg != 111 && std::abs(pdg) != 211) break;
    if (pdg == 22 && saw_pi0) break;
    particle = unique_parent(particle);
  }
}
void EtaCascadeTruth::refine(ClusterTruthMatch& match, const RawCluster* cluster, TowerInfoContainer* towers, RawTowerContainer* raw,
                             PHG4CellContainer* cells, PHG4HitContainer* hits) const
{
  std::vector<TruthContributor> refined;
  for (auto contributor : match.contributors)
  {
    annotate(contributor);
    const int mode = contributor.eta_decay_mode;
    if (contributor.g4_pdg_id != 221 || (mode != 2 && mode != 3)) { refined.push_back(contributor); continue; }
    const auto found = children_.find(contributor.g4_track_id);
    std::vector<std::array<int, 2>> tracks;
    for (const auto* child : found->second) tracks.push_back({child->get_track_id(), 0});
    const auto matches = Pi0ClusterTruthMatcher().match_many(cluster, towers, raw, cells, hits, truth_, tracks, true);
    const auto eta_match = Pi0ClusterTruthMatcher().match(cluster, towers, raw, cells, hits, truth_, {contributor.g4_track_id, 0}, true);
    double children_edep = 0.0;
    bool complete = matches.size() == tracks.size() && eta_match.status == Pi0ClusterTruthMatchStatus::complete;
    for (const auto& item : matches) complete = complete && item.status == Pi0ClusterTruthMatchStatus::complete;
    if (complete) for (std::size_t i = 0; i < matches.size(); ++i) children_edep += matches[i].gamma_edep[0];
    const double eta_edep = complete ? eta_match.gamma_edep[0] : 0.0;
    // Preserve the original shower energy and denominator. Only its measured daughter shares change.
    if (!complete || !(eta_edep > 0.0) || std::abs(children_edep - eta_edep) > 1e-5 * eta_edep)
    {
      contributor.eta_split_status = -1;
      refined.push_back(contributor);
      continue;
    }
    for (std::size_t i = 0; i < found->second.size(); ++i)
    {
      const auto* child = found->second[i];
      const double share = matches[i].gamma_edep[0] / children_edep;
      if (!(share > 0.0)) continue;
      auto daughter = contributor;
      daughter.g4_track_id = child->get_track_id(); daughter.g4_pdg_id = child->get_pid();
      daughter.hepmc_barcode = 0; daughter.hepmc_valid = 0; daughter.hepmc_pdg_id = -999;
      daughter.photon_category = daughter.photon_source = daughter.immediate_parent_pdg = daughter.classification_parent_pdg = -999;
      daughter.edep = contributor.edep * share; daughter.fraction = contributor.fraction * share; daughter.eta_split_status = 1;
      if (child->get_pid() == 111) daughter.pi0_parent_track_id = child->get_track_id();
      refined.push_back(daughter);
    }
  }
  std::sort(refined.begin(), refined.end(), [](const auto& a, const auto& b) { return a.edep == b.edep ? a.g4_track_id < b.g4_track_id : a.edep > b.edep; });
  match.contributors = std::move(refined);
}
}
