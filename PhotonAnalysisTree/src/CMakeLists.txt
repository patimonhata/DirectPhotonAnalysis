cmake_minimum_required(VERSION 3.16)
project(PhotonAnalysisTree LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED ENV{OFFLINE_MAIN})
  message(FATAL_ERROR "OFFLINE_MAIN is not set; source the sPHENIX environment first")
endif()

set(SPHENIX_PREFIX "$ENV{OFFLINE_MAIN}")
if(DEFINED PHOTON_CONFIGURED_PREFIX AND NOT PHOTON_CONFIGURED_PREFIX STREQUAL SPHENIX_PREFIX)
  message(FATAL_ERROR "This build directory belongs to another sPHENIX release; configure a new build directory")
endif()
set(PHOTON_CONFIGURED_PREFIX "${SPHENIX_PREFIX}" CACHE INTERNAL "Configured sPHENIX release")
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist EG TMVA TMVAUtils)

find_path(
  PI0RECONSTRUCTION_INCLUDE_DIR
  ShowerShapeCalculator.h
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../../Pi0Reconstruction/install/include/pi0reconstruction"
)
find_library(
  PI0RECONSTRUCTION_LIBRARY
  NAMES Pi0Reconstruction
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../../Pi0Reconstruction/install/lib"
)
if(NOT PI0RECONSTRUCTION_INCLUDE_DIR OR NOT PI0RECONSTRUCTION_LIBRARY)
  message(FATAL_ERROR "Build and install ../../Pi0Reconstruction before PhotonAnalysisTree")
endif()

foreach(library phool SubsysReco fun4all calo_io g4dst g4eval phhepmc HepMC)
  string(TOUPPER "${library}" library_upper)
  find_library(
    ${library_upper}_LIBRARY
    NAMES ${library}
    HINTS "${SPHENIX_PREFIX}/lib" "${SPHENIX_PREFIX}/lib64"
    REQUIRED
  )
  list(APPEND SPHENIX_LIBRARIES "${${library_upper}_LIBRARY}")
endforeach()

add_library(
  PhotonAnalysisTree SHARED
  DecayOriginDiagnostics.cc
  EtaCascadeTruth.cc
  HepMCPhotonClassifier.cc
  Pi0AnchorTopologyEvaluator.cc
  Pi0ClusterTruthMatcher.cc
  PhotonAnalysisTree.cc
  PythiaClusterEtSpectrum.cc
  PythiaPi0ClusterMultiplicity.cc
  PhotonTreeCommon.cc
  PythiaClusterTruthMatcher.cc
  PythiaPhotonCandidateTree.cc
  PythiaPhotonAnalysisTree.cc
  PythiaTruthPtSpectrum.cc
  PythiaTruthSpectrumTree.cc
)
target_include_directories(
  PhotonAnalysisTree
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    "${SPHENIX_PREFIX}/include"
    "${PI0RECONSTRUCTION_INCLUDE_DIR}"
)
target_link_libraries(
  PhotonAnalysisTree
  PRIVATE
    ${SPHENIX_LIBRARIES}
    "${PI0RECONSTRUCTION_LIBRARY}"
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::EG
    ROOT::TMVA
    ROOT::TMVAUtils
)
target_compile_options(PhotonAnalysisTree PRIVATE -Wall -Wextra -Wpedantic)

find_path(
  ONNXRUNTIME_INCLUDE_DIR
  onnxruntime_cxx_api.h
  HINTS "${SPHENIX_PREFIX}/include" "${SPHENIX_PREFIX}/include/onnxruntime"
)
find_library(
  ONNXRUNTIME_LIBRARY
  NAMES onnxruntime
  HINTS "${SPHENIX_PREFIX}/lib" "${SPHENIX_PREFIX}/lib64"
)
if(NOT ONNXRUNTIME_INCLUDE_DIR OR NOT ONNXRUNTIME_LIBRARY)
  message(FATAL_ERROR "ONNX Runtime was not found in OFFLINE_MAIN")
endif()

add_library(Pi0GammaOnnx SHARED Pi0GammaOnnx.cc)
target_include_directories(
  Pi0GammaOnnx
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE "${ONNXRUNTIME_INCLUDE_DIR}"
)
target_link_libraries(Pi0GammaOnnx PRIVATE "${ONNXRUNTIME_LIBRARY}")
target_compile_options(Pi0GammaOnnx PRIVATE -Wall -Wextra -Wpedantic)

set_target_properties(PhotonAnalysisTree Pi0GammaOnnx PROPERTIES INSTALL_RPATH "$ORIGIN/../dependencies/Pi0Reconstruction/lib")
install(TARGETS PhotonAnalysisTree Pi0GammaOnnx)
install(
  FILES
    DecayOriginDiagnostics.h
    EtaCascadeTruth.h
    EtaOriginFractions.h
    HepMCPhotonClassifier.h
    Pi0AnchorTopologyEvaluator.h
    Pi0ClusterTruthMatcher.h
    PhotonAnalysisTree.h
    PhotonTreeCommon.h
    Pi0GammaOnnx.h
    PythiaClusterTruthMatcher.h
    PythiaClusterEtSpectrum.h
    PythiaPi0ClusterMultiplicity.h
    PythiaPhotonAnalysisTree.h
    PythiaPhotonCandidateTree.h
    PythiaTruthPtSpectrum.h
    PythiaTruthSpectrumTree.h
    "${PI0RECONSTRUCTION_INCLUDE_DIR}/ShowerShapeCalculator.h"
  DESTINATION include
)
install(FILES
  ../models/best_model.onnx
  ../models/onnx_metadata.json
  ../models/model_ppg15_nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_single_tmva.root
  DESTINATION ${CMAKE_INSTALL_DATADIR}/PhotonAnalysisTree/models)
