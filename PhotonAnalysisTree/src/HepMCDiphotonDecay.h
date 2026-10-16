#ifndef PHOTON_TREE_HEPMC_DIPHOTON_DECAY_H
#define PHOTON_TREE_HEPMC_DIPHOTON_DECAY_H

#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <set>

namespace photon_tree
{
inline bool has_hepmc_two_photon_decay(const HepMC::GenParticle* parent)
{
  const auto* vertex = parent ? parent->end_vertex() : nullptr;
  if (!vertex) return false;
  unsigned int count = 0;
  for (auto iterator = vertex->particles_out_const_begin(); iterator != vertex->particles_out_const_end(); ++iterator)
  {
    if (!*iterator || (*iterator)->pdg_id() != 22) return false;
    ++count;
  }
  return count == 2U;
}

// Follow unique photon bookkeeping copies, stopping at the first physical parent.
inline const HepMC::GenParticle* hepmc_photon_parent(const HepMC::GenParticle* photon, int parent_pdg)
{
  if (!photon || photon->pdg_id() != 22) return nullptr;
  const HepMC::GenParticle* current = photon;
  std::set<int> visited;
  if (current->barcode() != 0) visited.insert(current->barcode());
  while (current)
  {
    const auto* vertex = current->production_vertex();
    if (!vertex) return nullptr;
    const HepMC::GenParticle* parent = nullptr;
    for (auto iterator = vertex->particles_in_const_begin(); iterator != vertex->particles_in_const_end(); ++iterator)
    {
      if (!*iterator) continue;
      if (parent) return nullptr;
      parent = *iterator;
    }
    if (!parent) return nullptr;
    if (parent->pdg_id() != 22) return parent->pdg_id() == parent_pdg ? parent : nullptr;
    current = parent;
    if (current->barcode() != 0 && !visited.insert(current->barcode()).second) return nullptr;
  }
  return nullptr;
}
}
#endif
