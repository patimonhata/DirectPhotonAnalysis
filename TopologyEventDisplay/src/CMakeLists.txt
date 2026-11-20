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

set(PHOTON_ANALYSIS_TREE_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../../PhotonAnalysisTree/install/dev/ana.565"
    CACHE PATH "PhotonAnalysisTree install to use for this build")
unset(PHOTON_ANALYSIS_TREE_INCLUDE_DIR CACHE)
unset(PHOTON_ANALYSIS_TREE_LIBRARY CACHE)
find_path(PHOTON_ANALYSIS_TREE_INCLUDE_DIR PhotonCandidateTopology.h
          PATHS "${PHOTON_ANALYSIS_TREE_PREFIX}/include" NO_DEFAULT_PATH REQUIRED)
find_library(PHOTON_ANALYSIS_TREE_LIBRARY NAMES PhotonAnalysisTree
             PATHS "${PHOTON_ANALYSIS_TREE_PREFIX}/lib64" "${PHOTON_ANALYSIS_TREE_PREFIX}/lib" NO_DEFAULT_PATH REQUIRED)
get_filename_component(PHOTON_ANALYSIS_TREE_LIB_DIR "${PHOTON_ANALYSIS_TREE_LIBRARY}" DIRECTORY)

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

set_target_properties(TopologyEventDisplay PROPERTIES INSTALL_RPATH "${PHOTON_ANALYSIS_TREE_LIB_DIR}")
install(TARGETS TopologyEventDisplay)
install(FILES TopologyEventDisplayDump.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
