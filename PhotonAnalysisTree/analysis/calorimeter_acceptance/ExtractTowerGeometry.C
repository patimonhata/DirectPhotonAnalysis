#include <G4_CEmc_Spacal.C>
#include <G4_HcalIn_ref.C>
#include <G4_HcalOut_ref.C>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom_Spacalv3.h>
#include <g4main/PHG4HitContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/PHIODataNode.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <TFile.h>
#include <TNamed.h>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>

namespace calo_acceptance
{
void require(bool condition, const std::string& message)
{
  if (!condition) throw std::runtime_error(message);
}

void extract(const std::string& output)
{
  const std::string release = gSystem->Getenv("OFFLINE_MAIN") ? gSystem->Getenv("OFFLINE_MAIN") : "";
  require(release.substr(release.find_last_of('/') + 1) == "ana.565", "Set up ana.565 before extracting geometry");
  require(gSystem->mkdir(output.c_str(), true) == 0 || !gSystem->AccessPathName(output.c_str()), "Cannot create output directory");
  auto top = std::make_unique<PHCompositeNode>("TOP");
  auto* run = new PHCompositeNode("RUN");
  auto* dst = new PHCompositeNode("DST");
  top->addNode(run); top->addNode(dst); top->addNode(new PHCompositeNode("PAR"));

  // Run the installed standard setup macros only to register configured subsystems.
  // Do not initialize Geant4, load a CDB payload, or transport any particles.
  PHG4Reco setup;
  CEmc(&setup, 0., 0);
  HCalInner(&setup, 0., 0);
  HCalOuter(&setup, 0., 0);
  auto subsystem = [&](const std::string& name) {
    auto* result = dynamic_cast<PHG4DetectorSubsystem*>(setup.getSubsystem(name + "_0"));
    require(result != nullptr, "Missing configured subsystem " + name);
    return result;
  };
  auto* cemc = subsystem("CEMC");
  require(cemc->ReadParamsFromFile("CEMC", PHG4DetectorSubsystem::xml, 1) == 0, "Cannot read CEMC geometry XML");
  cemc->UpdateParametersWithMacro();
  auto* cylinder = new PHG4CylinderGeom_Spacalv3();
  cylinder->ImportParameters(*cemc->GetParams());
  cylinder->set_layer(0);
  auto* cylinders = new PHG4CylinderGeomContainer();
  cylinders->AddLayerGeom(0, cylinder);
  run->addNode(new PHIODataNode<PHObject>(cylinders, "CYLINDERGEOM_CEMC", "PHObject"));
  dst->addNode(new PHIODataNode<PHObject>(new PHG4HitContainer("G4HIT_CEMC"), "G4HIT_CEMC", "PHObject"));
  PHG4FullProjSpacalCellReco cells;
  cells.Detector("CEMC");
  require(cells.InitRun(top.get()) == 0, "CEMC cell geometry initialization failed");
  RawTowerBuilder towers;
  towers.Detector("CEMC");
  towers.set_towerinfo(RawTowerBuilder::kRawTowerOnly);
  require(towers.InitRun(top.get()) == 0, "CEMC tower geometry initialization failed");

  for (const std::string name : {"HCALIN", "HCALOUT"})
  {
    auto* detector = subsystem(name);
    detector->UpdateParametersWithMacro();
    detector->GetParams()->SaveToNodeTree(run, "G4GEOPARAM_" + name, 0);
    HcalRawTowerBuilder builder;
    builder.Detector(name);
    builder.set_double_param("phistart", detector->GetParams()->get_double_param("phistart"));
    builder.set_towerinfo(HcalRawTowerBuilder::kRawTowerOnly);
    require(builder.InitRun(top.get()) == 0, name + " tower geometry initialization failed");
  }

  TFile snapshot((output + "/tower_geometry.root").c_str(), "RECREATE");
  require(!snapshot.IsZombie(), "Cannot write tower geometry snapshot");
  TNamed("release", release.c_str()).Write();
  TNamed("definition", "Full-phi outer envelope of tower eta bins at the common tower-center radius; x_v=y_v=0; internal gaps ignored").Write();
  std::ofstream summary(output + "/geometry.csv"), bins(output + "/eta_bins.csv");
  require(summary.good() && bins.good(), "Cannot write geometry CSV files");
  summary << std::setprecision(17) << "detector,radius_cm,eta_min,eta_max,z_min_cm,z_max_cm,eta_bins,phi_bins,container_radius_cm,container_thickness_cm\n";
  bins << std::setprecision(17) << "detector,bin,eta_low,eta_high\n";
  for (const std::string name : {"CEMC", "HCALIN", "HCALOUT"})
  {
    auto* geom = findNode::getClass<RawTowerGeomContainer>(top.get(), "TOWERGEOM_" + name);
    require(geom && geom->get_etabins() > 0 && geom->get_phibins() > 0, "Missing/empty tower geometry " + name);
    double radius = -1., etaMin = 1e9, etaMax = -1e9;
    std::size_t count = 0;
    const auto range = geom->get_tower_geometries();
    for (auto it = range.first; it != range.second; ++it)
    {
      auto* tower = it->second;
      const double r = std::hypot(tower->get_center_x(), tower->get_center_y());
      if (radius < 0.) radius = r;
      require(std::isfinite(r) && r > 0. && std::abs(r - radius) < 1e-7, name + " is not a common-radius tower surface");
      require(std::abs(std::asinh(tower->get_center_z() / r) - geom->get_etacenter(tower->get_bineta())) < 1e-7,
              name + " tower-center eta disagrees with eta-bin geometry");
      ++count;
    }
    require(count == static_cast<std::size_t>(geom->get_etabins() * geom->get_phibins()), name + " incomplete tower grid");
    for (int bin = 0; bin < geom->get_etabins(); ++bin)
    {
      const auto edge = geom->get_etabounds(bin);
      require(std::isfinite(edge.first) && std::isfinite(edge.second) && edge.first < edge.second, name + " invalid eta bounds");
      etaMin = std::min(etaMin, edge.first); etaMax = std::max(etaMax, edge.second);
      bins << name << ',' << bin << ',' << edge.first << ',' << edge.second << '\n';
    }
    summary << name << ',' << radius << ',' << etaMin << ',' << etaMax << ',' << radius * std::sinh(etaMin) << ','
            << radius * std::sinh(etaMax) << ',' << geom->get_etabins() << ',' << geom->get_phibins() << ','
            << geom->get_radius() << ',' << geom->get_thickness() << '\n';
    geom->Write(("TOWERGEOM_" + name).c_str());
    std::cout << name << ": R=" << radius << " cm, eta=[" << etaMin << ", " << etaMax << "]" << std::endl;
  }
  snapshot.Close();
}
}  // namespace calo_acceptance
