#ifndef RYOTARO_HEPMCPHOTONCLASSIFIER_H_20260810
#define RYOTARO_HEPMCPHOTONCLASSIFIER_H_20260810

namespace HepMC
{
class GenParticle;
}

namespace photon_tree
{
struct HepMCPhotonClassification
{
  bool valid = false;
  int category = -999;
  int source = -999;
  int immediate_parent_pdg = -999;
  int classification_parent_pdg = -999;
};

class HepMCPhotonClassifier
{
 public:
  static constexpr int kAlgorithmVersion = 1;

  // category: -1 non-photon, 0 unknown, 1 direct, 2 fragmentation, 3 decay.
  // source:   -1 non-photon, 0 other, 1 photon, 2 pi0, 3 eta parent.
  HepMCPhotonClassification classify(const HepMC::GenParticle* particle) const;
};
}

#endif
