cmake_minimum_required(VERSION 3.16)
project(TopoClusterHCalStudy LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED ENV{OFFLINE_MAIN})
  message(FATAL_ERROR "OFFLINE_MAIN is not set; source the sPHENIX environment first")
endif()

set(SPHENIX_PREFIX "$ENV{OFFLINE_MAIN}")
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

foreach(library phool SubsysReco fun4all calo_io)
  string(TOUPPER "${library}" library_upper)
  find_library(
    ${library_upper}_LIBRARY
    NAMES ${library}
    HINTS "${SPHENIX_PREFIX}/lib" "${SPHENIX_PREFIX}/lib64"
    REQUIRED
  )
  list(APPEND SPHENIX_LIBRARIES "${${library_upper}_LIBRARY}")
endforeach()

add_library(TopoClusterHCalStudy SHARED TopoClusterHCalTree.cc)
target_include_directories(
  TopoClusterHCalStudy
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE "${SPHENIX_PREFIX}/include"
)
target_link_libraries(
  TopoClusterHCalStudy
  PRIVATE
    ${SPHENIX_LIBRARIES}
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
)
target_compile_options(TopoClusterHCalStudy PRIVATE -Wall -Wextra -Wpedantic)

install(TARGETS TopoClusterHCalStudy)
install(FILES TopoClusterHCalTree.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
