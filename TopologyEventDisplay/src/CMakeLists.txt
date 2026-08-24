cmake_minimum_required(VERSION 3.16)
project(TopologyEventDisplay LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED ENV{OFFLINE_MAIN})
  message(FATAL_ERROR "OFFLINE_MAIN is not set; source the sPHENIX environment first")
endif()
set(SPHENIX_PREFIX "$ENV{OFFLINE_MAIN}")
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

find_path(
  PHOTON_ANALYSIS_TREE_INCLUDE_DIR
  Pi0AnchorTopologyEvaluator.h
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../../PhotonAnalysisTree/install/include"
  REQUIRED
)
find_library(
  PHOTON_ANALYSIS_TREE_LIBRARY
  NAMES PhotonAnalysisTree
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../../PhotonAnalysisTree/install/lib64"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../PhotonAnalysisTree/install/lib"
  REQUIRED
)

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

add_library(TopologyEventDisplay SHARED TopologyEventDisplayDump.cc)
target_include_directories(
  TopologyEventDisplay
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    "${SPHENIX_PREFIX}/include"
    "${PHOTON_ANALYSIS_TREE_INCLUDE_DIR}"
)
target_link_libraries(
  TopologyEventDisplay
  PRIVATE
    "${PHOTON_ANALYSIS_TREE_LIBRARY}"
    ${SPHENIX_LIBRARIES}
    ROOT::Core ROOT::RIO ROOT::Tree
)
target_compile_options(TopologyEventDisplay PRIVATE -Wall -Wextra -Wpedantic)

install(TARGETS TopologyEventDisplay)
install(FILES TopologyEventDisplayDump.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
