#include "HepMCPhotonClassifier.h"

#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace {
  std::vector<const HepMC::GenParticle*> incoming(const HepMC::GenVertex* vertex) {
    std::vector<const HepMC::GenParticle*> result;
    if (!vertex) {
      return result;
    }
    for (auto iter = vertex->particles_in_const_begin(); iter != vertex->particles_in_const_end(); ++iter) {
      if (*iter) {
        result.push_back(*iter);
      }
    }
    return result;
  }

  std::vector<const HepMC::GenParticle*> outgoing(const HepMC::GenVertex* vertex) {
    std::vector<const HepMC::GenParticle*> result;
    if (!vertex) {
      return result;
    }
    for (auto iter = vertex->particles_out_const_begin(); iter != vertex->particles_out_const_end(); ++iter) {
      if (*iter) {
        result.push_back(*iter);
      }
    }
    return result;
  }
}

namespace photon_tree
{
  HepMCPhotonClassification HepMCPhotonClassifier::classify(const HepMC::GenParticle* particle) const {
    HepMCPhotonClassification result;
    if (!particle) {
      return result;
    }

    result.valid = true;
    if (particle->pdg_id() != 22) {
      result.category = -1;
      result.source = -1;
      result.immediate_parent_pdg = 0;
      result.classification_parent_pdg = 0;
      return result;
    }

    result.category = 0;
    result.source = 0;
    const auto immediate_parents = incoming(particle->production_vertex());
    if (!immediate_parents.empty()) {
      result.immediate_parent_pdg = immediate_parents.front()->pdg_id();
    } else {
      result.immediate_parent_pdg = 0;
    }

    for (const HepMC::GenParticle* parent : immediate_parents) {
      if (parent->pdg_id() == 22) {
        result.source = 1;
        break;
      }
      if (parent->pdg_id() == 111) {
        result.source = 2;
      } else if (parent->pdg_id() == 221 && result.source == 0) {
        result.source = 3;
      }
    }

    // Pythia can store photon->photon bookkeeping copies. Follow the unique
    // photon parent backwards before classifying the physical production vertex.
    const HepMC::GenParticle* current = particle;
    std::set<int> visited_barcodes;
    if (current->barcode() != 0) {
      visited_barcodes.insert(current->barcode());
    }
    std::vector<const HepMC::GenParticle*> parents;
    const HepMC::GenVertex* vertex = nullptr;
    while (current) {
      vertex = current->production_vertex();
      parents = incoming(vertex);
      if (parents.size() != 1 || parents.front()->pdg_id() != 22) {
        break;
      }
      const HepMC::GenParticle* parent = parents.front();
      if (parent->barcode() != 0 && !visited_barcodes.insert(parent->barcode()).second) {
        result.valid = false;
        result.category = -999;
        return result;
      }
      current = parent;
    }

    if (!vertex) {
      return result;
    }
    const auto children = outgoing(vertex);
    result.classification_parent_pdg = parents.empty() ? 0 : parents.front()->pdg_id();
    const bool has_photon = std::any_of(children.begin(), children.end(), [](const auto* child) {
      return child->pdg_id() == 22;
    });
    if (!has_photon) {
      return result;
    }

    // This intentionally reproduces the analysis-specific rules used by
    // PhotonAna::photon_type; it is not an sPHENIX framework definition.
    if (parents.size() == 2 && children.size() == 2 &&
        std::abs(parents[0]->pdg_id()) <= 22 && std::abs(parents[1]->pdg_id()) <= 22 &&
        std::abs(children[0]->pdg_id()) <= 22 && std::abs(children[1]->pdg_id()) <= 22)
    {
      result.category = 1;
    } 
    else if (parents.size() == 1) {
      const int parent_pdg = parents.front()->pdg_id();
      if (std::abs(parent_pdg) <= 11 && children.size() == 2) {
        const bool retains_parent = std::any_of(children.begin(), children.end(), [parent_pdg](const auto* child) {
          return child->pdg_id() == parent_pdg;
        });
        if (retains_parent) {
          result.category = 2;
        }
      }
      if (std::abs(parent_pdg) > 37) {
        result.category = 3;
      }
    }
    return result;
  }
}
