#ifndef RYOTARO_SHOWER_SHAPE_CALCULATOR_H_20260720
#define RYOTARO_SHOWER_SHAPE_CALCULATOR_H_20260720

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

class RawCluster;
class TowerInfoContainer;

class ShowerShapeCalculator
{
 public:
  static constexpr int kAlgorithmVersion = 1;
  static constexpr int kCemcEtaBins = 96;
  static constexpr int kCemcPhiBins = 256;
  static constexpr int kPatchSide = 7;
  static constexpr std::size_t kPatchSize = kPatchSide * kPatchSide;

  struct Config
  {
    float min_tower_energy = 0.07F;
  };

  struct Result
  {
    bool valid = false;
    bool full_containment = false;
    bool edge_padded = false;
    bool tower_data_complete = false;

    float cog_ieta = std::numeric_limits<float>::quiet_NaN();
    float cog_iphi = std::numeric_limits<float>::quiet_NaN();
    float cluster_energy_above_threshold = 0.0F;
    float owned_patch_energy = 0.0F;

    float w_eta_cogx = 0.0F;
    float w_phi_cogx = 0.0F;

    float e11 = 0.0F;
    float e33 = 0.0F;
    float e32 = 0.0F;
    float e35 = 0.0F;
    float e11_over_e33 = 0.0F;
    float e32_over_e35 = 0.0F;

    // These are the E_t1...E_t4 variables in the photon-ID table. The
    // "t" labels normalized 2x2 combinations; these are not transverse energies.
    float et1 = 0.0F;
    float et2 = 0.0F;
    float et3 = 0.0F;
    float et4 = 0.0F;

    // Raw calibrated energy is retained even for bad/below-threshold towers.
    // The good and ownership masks make the calculation reproducible later.
    std::array<float, kPatchSize> patch_energy = {};
    std::array<std::uint8_t, kPatchSize> patch_good = {};
    std::array<std::uint8_t, kPatchSize> patch_owned = {};
  };

  ShowerShapeCalculator() = default;
  explicit ShowerShapeCalculator(const Config &config);

  Result calculate(const RawCluster &cluster, TowerInfoContainer &cemc_towers) const;

 private:
  Config config_;
};

#endif  // RYOTARO_SHOWER_SHAPE_CALCULATOR_H_20260720
