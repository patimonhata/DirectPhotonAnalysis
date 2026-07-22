#include "Pi0GammaOnnx.h"

#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>

namespace
{
  constexpr const char* input_names[] = {
      "global_features", "point_features", "point_mask"};
  constexpr const char* output_names[] = {"p_gamma"};

  template <std::size_t Size>
  void require_finite(const std::array<float, Size>& values, const char* description)
  {
    if (!std::all_of(values.begin(), values.end(), [](float value) { return std::isfinite(value); }))
    {
      throw std::invalid_argument(std::string(description) + " contains a non-finite value");
    }
  }
}

namespace pi0gamma
{
  class Pi0GammaOnnx::Impl
  {
   public:
    Impl(const std::string& model_path, int intra_op_threads)
      : environment_(ORT_LOGGING_LEVEL_WARNING, "Pi0GammaOnnx")
      , session_options_()
      , session_(nullptr)
      , memory_info_(Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault))
    {
      if (intra_op_threads <= 0)
      {
        throw std::invalid_argument("intra_op_threads must be positive");
      }
      session_options_.SetIntraOpNumThreads(intra_op_threads);
      session_options_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);
      session_ = Ort::Session(environment_, model_path.c_str(), session_options_);
      if (session_.GetInputCount() != 3 || session_.GetOutputCount() != 1)
      {
        throw std::runtime_error("Unexpected Pi0Gamma ONNX input/output count");
      }
    }

    float predict(const GlobalFeatures& global_features,
                  const std::vector<TowerFeatures>& tower_features) const
    {
      require_finite(global_features, "global_features");
      if (tower_features.empty())
      {
        throw std::invalid_argument("At least one tower is required");
      }

      std::vector<float> flat_tower_features;
      flat_tower_features.reserve(tower_features.size() * TowerFeatures{}.size());
      for (const TowerFeatures& tower : tower_features)
      {
        require_finite(tower, "tower_features");
        flat_tower_features.insert(flat_tower_features.end(), tower.begin(), tower.end());
      }

      const int64_t number_of_towers = static_cast<int64_t>(tower_features.size());
      const std::array<int64_t, 2> global_shape = {1, 3};
      const std::array<int64_t, 3> point_shape = {1, number_of_towers, 4};
      const std::array<int64_t, 2> mask_shape = {1, number_of_towers};
      GlobalFeatures global_buffer = global_features;
      std::unique_ptr<bool[]> point_mask = std::make_unique<bool[]>(tower_features.size());
      std::fill_n(point_mask.get(), tower_features.size(), true);

      std::array<Ort::Value, 3> inputs = {
          Ort::Value::CreateTensor<float>(memory_info_, global_buffer.data(), global_buffer.size(),
                                          global_shape.data(), global_shape.size()),
          Ort::Value::CreateTensor<float>(memory_info_, flat_tower_features.data(), flat_tower_features.size(),
                                          point_shape.data(), point_shape.size()),
          Ort::Value::CreateTensor<bool>(memory_info_, point_mask.get(), tower_features.size(),
                                         mask_shape.data(), mask_shape.size())};

      std::vector<Ort::Value> outputs = session_.Run(
          Ort::RunOptions{nullptr}, input_names, inputs.data(), inputs.size(),
          output_names, std::size(output_names));
      if (outputs.size() != 1 || !outputs.front().IsTensor() ||
          outputs.front().GetTensorTypeAndShapeInfo().GetElementCount() != 1)
      {
        throw std::runtime_error("Pi0Gamma ONNX returned an invalid output");
      }

      const float probability = outputs.front().GetTensorData<float>()[0];
      if (!std::isfinite(probability) || probability < 0.0F || probability > 1.0F)
      {
        throw std::runtime_error("Pi0Gamma ONNX returned an invalid probability");
      }
      return probability;
    }

   private:
    Ort::Env environment_;
    Ort::SessionOptions session_options_;
    mutable Ort::Session session_;
    Ort::MemoryInfo memory_info_;
  };

  Pi0GammaOnnx::Pi0GammaOnnx(const std::string& model_path, int intra_op_threads)
    : impl_(std::make_unique<Impl>(model_path, intra_op_threads))
  {
  }

  Pi0GammaOnnx::~Pi0GammaOnnx() = default;
  Pi0GammaOnnx::Pi0GammaOnnx(Pi0GammaOnnx&&) noexcept = default;
  Pi0GammaOnnx& Pi0GammaOnnx::operator=(Pi0GammaOnnx&&) noexcept = default;

  float Pi0GammaOnnx::predict(const GlobalFeatures& global_features,
                              const std::vector<TowerFeatures>& tower_features) const
  {
    return impl_->predict(global_features, tower_features);
  }
}
