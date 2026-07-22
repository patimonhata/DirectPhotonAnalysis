#ifndef PI0GAMMAONNX_H
#define PI0GAMMAONNX_H

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace pi0gamma
{
  class Pi0GammaOnnx
  {
   public:
    using GlobalFeatures = std::array<float, 3>;
    using TowerFeatures = std::array<float, 4>;

    explicit Pi0GammaOnnx(const std::string& model_path, int intra_op_threads = 1);
    ~Pi0GammaOnnx();

    Pi0GammaOnnx(Pi0GammaOnnx&&) noexcept;
    Pi0GammaOnnx& operator=(Pi0GammaOnnx&&) noexcept;
    Pi0GammaOnnx(const Pi0GammaOnnx&) = delete;
    Pi0GammaOnnx& operator=(const Pi0GammaOnnx&) = delete;

    float predict(const GlobalFeatures& global_features,
                  const std::vector<TowerFeatures>& tower_features) const;

   private:
    class Impl;
    std::unique_ptr<Impl> impl_;
  };
}

#endif
