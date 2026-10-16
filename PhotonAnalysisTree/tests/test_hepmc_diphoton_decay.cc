#include "HepMCDiphotonDecay.h"
#include <HepMC/GenEvent.h>
#include <cassert>
#include <iostream>

int main()
{
  HepMC::GenEvent event;
  auto* vertex = new HepMC::GenVertex();
  auto* eta = new HepMC::GenParticle(HepMC::FourVector(0, 0, 0, 1), 221, 2);
  auto* gamma0 = new HepMC::GenParticle(HepMC::FourVector(0, 0, .5, .5), 22, 1);
  auto* gamma1 = new HepMC::GenParticle(HepMC::FourVector(0, 0, -.5, .5), 22, 1);
  vertex->add_particle_in(eta); vertex->add_particle_out(gamma0); vertex->add_particle_out(gamma1); event.add_vertex(vertex);
  assert(photon_tree::has_hepmc_two_photon_decay(eta));
  assert(photon_tree::hepmc_photon_parent(gamma0, 221) == eta);
  assert(photon_tree::hepmc_photon_parent(gamma1, 111) == nullptr);
  auto* copies = new HepMC::GenVertex();
  auto* copy = new HepMC::GenParticle(gamma0->momentum(), 22, 1);
  copies->add_particle_in(gamma0); copies->add_particle_out(copy); event.add_vertex(copies);
  assert(photon_tree::hepmc_photon_parent(copy, 221) == eta);
  assert(photon_tree::has_hepmc_two_photon_decay(eta));
  // Two transported photons do not imply an exclusive two-photon HepMC decay.
  auto* extra = new HepMC::GenParticle(HepMC::FourVector(), 111, 2);
  vertex->add_particle_out(extra);
  assert(!photon_tree::has_hepmc_two_photon_decay(eta));
  auto* pion_decay = new HepMC::GenVertex();
  auto* pion_gamma0 = new HepMC::GenParticle(HepMC::FourVector(), 22, 1);
  auto* pion_gamma1 = new HepMC::GenParticle(HepMC::FourVector(), 22, 1);
  pion_decay->add_particle_in(extra); pion_decay->add_particle_out(pion_gamma0); pion_decay->add_particle_out(pion_gamma1); event.add_vertex(pion_decay);
  assert(photon_tree::hepmc_photon_parent(pion_gamma0, 111) == extra);
  assert(photon_tree::hepmc_photon_parent(pion_gamma0, 221) == nullptr);
  assert(!photon_tree::has_hepmc_two_photon_decay(nullptr));
  assert(photon_tree::hepmc_photon_parent(eta, 221) == nullptr);
  std::cout << "HepMC direct eta decay, photon copies, extra daughters and pi0 cascade: passed" << std::endl;
}
