#include <sys/stat.h>
#include <filesystem>
#include <system_error>
#include <iostream>
#include <cstdlib> // for EXIT_FAILURE

#include <G4Setup_sPHENIX.C> // for G4Init(), G4Setup()

#include <DisplayOn.C>
#include <G4Setup_sPHENIX.C>
#include <G4_Mbd.C>
#include <G4_CaloTrigger.C>
#include <G4_Centrality.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_KFParticle.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>
// #include <G4_Tracking.C>
#include <G4_User.C>
#include <QA.C>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <RawClusterBuilderTemplate.h>
#include "Calo_Calib_ryotaro.C"

#include "/sphenix/u/ryotaro/DirectPhotonAnalysis/EmcalEtViewer/install/include/displaylegoplot/DisplayLegoPlot.h"


void ensure_dir(const std::string& path);

int Fun4All_SingleParticlePi0(int processID=0, int nEvents=5, bool save_tree=false) {
  // const int nEvents = 5;
  bool runTruth = false;
  std::string particle_name = "pi0";

  std::ostringstream pid;
  pid << std::setw(6) << std::setfill('0') << processID;
  std::string pid_str = pid.str();

  // char cwd[PATH_MAX];
  // char* cwd_result = getcwd(cwd, sizeof(cwd));
  // if ( cwd_result == nullptr ) {
	//   std::cerr << "Failed to get current working directory!" << std::endl;
	//   return EXIT_FAILURE;
  // }

  std::error_code ec;
  std::filesystem::path cwd = std::filesystem::current_path(ec);
  if (ec) {
    std::cerr << "Failed to get current working directory: " << ec.message() << std::endl;
    return EXIT_FAILURE;
  }

  // std::string baseDir(cwd);
  std::string baseDir("/sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output");
  std::string outDir  = baseDir + "/DST_" + particle_name;
  std::string outDir2 = baseDir + "/ana_" + particle_name;
  std::string outDir3 = baseDir + "/jobtime_" + particle_name;
  ensure_dir(outDir);
  ensure_dir(outDir2);
  ensure_dir(outDir3);
  ensure_dir(outDir + "/qa");
  
  Fun4AllServer *fun4allServer = Fun4AllServer::instance();

  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  PHRandomSeed::Verbosity(1);

  recoConsts *rc = recoConsts::instance();
  Input::VERBOSITY = 10;
  Input::SIMPLE = true;
  InputInit();  // This creates the input generator(s)

  INPUTGENERATOR::SimpleEventGenerator[0]->add_particles(particle_name, 1);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus,
                                                                            PHG4SimpleEventGenerator::Gaus,
                                                                            PHG4SimpleEventGenerator::Gaus);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(0., 0.);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(0.0, 0.0);
  // INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(9., 9.);

  // register all input generators with Fun4All
  InputRegister();

  // Simulation setup
  Enable::MBDFAKE = true;
  Enable::PIPE = true;
  Enable::PIPE_ABSORBER = true;
  Enable::MVTX = true;
  Enable::INTT = true;
  Enable::TPC = true;
  Enable::MICROMEGAS = true;

  Enable::MAGNET = true;
  Enable::MAGNET_ABSORBER = true;

  Enable::PLUGDOOR_ABSORBER = true;

  Enable::CEMC = true;
  Enable::CEMC_ABSORBER = true;
  Enable::CEMC_CELL = true;
  Enable::CEMC_TOWER = true;
  // Enable::CEMC_CLUSTER = true;
  // Enable::CEMC_EVAL =  false;
  // Enable::CEMC_QA =  false;

  Enable::HCALIN = true;
  Enable::HCALIN_ABSORBER = true;
  Enable::HCALIN_CELL =  true;
  Enable::HCALIN_TOWER =  true;
  // Enable::HCALIN_CLUSTER =  true;
  // Enable::HCALIN_EVAL =  true;
  // Enable::HCALIN_QA =  true;

  Enable::HCALOUT = true;
  Enable::HCALOUT_ABSORBER = true;
  Enable::HCALOUT_CELL =  true;
  Enable::HCALOUT_TOWER =  true;
  // Enable::HCALOUT_CLUSTER = false;
  // Enable::HCALOUT_EVAL =  false;
  // Enable::HCALOUT_QA =  false;

  Enable::GLOBAL_RECO = true;

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", CDB::global_tag);
  rc->set_uint64Flag("TIMESTAMP", CDB::timestamp);

  MagnetInit();

  G4Init(); // Initialize the selected subsystems
  if (!Input::READHITS)
  {
    G4Setup();
  }

  if (Enable::MBD || Enable::MBDFAKE) Mbd_Reco();
  if (Enable::MVTX_CELL) Mvtx_Cells();
  if (Enable::INTT_CELL) Intt_Cells();
  if (Enable::TPC_CELL) TPC_Cells();
  if (Enable::MICROMEGAS_CELL) Micromegas_Cells();
  if (Enable::CEMC_CELL) CEMC_Cells();
  if (Enable::HCALIN_CELL) HCALInner_Cells();
  if (Enable::HCALOUT_CELL) HCALOuter_Cells();
  if (Enable::CEMC_TOWER) CEMC_Towers();
  if (Enable::CEMC_CLUSTER) CEMC_Clusters();
  if (Enable::HCALIN_TOWER) HCALInner_Towers();
  if (Enable::HCALIN_CLUSTER) HCALInner_Clusters();
  if (Enable::HCALOUT_TOWER) HCALOuter_Towers();
  if (Enable::HCALOUT_CLUSTER) HCALOuter_Clusters();
  if (Enable::TOPOCLUSTER) TopoClusterReco();

  Process_Calo_Calib_ryotaro();

  // if (Enable::TRACKING_TRACK) TrackingInit();
  // if (Enable::MVTX_CLUSTER) Mvtx_Clustering();
  // if (Enable::INTT_CLUSTER) Intt_Clustering();
  // if (Enable::TPC_CLUSTER) TPC_Clustering();
  // if (Enable::MICROMEGAS_CLUSTER) Micromegas_Clustering();
  // if (Enable::TRACKING_TRACK) Tracking_Reco();
  if (Enable::GLOBAL_RECO && Enable::GLOBAL_FASTSIM) {
    cout << "You can only enable Enable::GLOBAL_RECO or Enable::GLOBAL_FASTSIM, not both" << endl;
    gSystem->Exit(1);
  }
  if (Enable::GLOBAL_RECO) {
    Global_Reco();
  } else if (Enable::GLOBAL_FASTSIM) {
    Global_FastSim();
  }
  if (Enable::CENTRALITY) Centrality();
  if (Enable::CALOTRIGGER) CaloTrigger_Sim();
  if (Enable::JETS) Jet_Reco();
  if (Enable::HIJETS) HIJetReco();
  if (Enable::PARTICLEFLOW) ParticleFlow();

  std::string outputName = outDir + "/DST_single_" + particle_name + "_";
  if (runTruth) {
  	outputName += "truth";
  } else {
    outputName += "reconstructed";
  }
  outputName += "Info_" + pid_str + ".root";
  std::string outputName2 = outDir2 + "/ana_" + pid_str + "ClusterBuilder_5events.root";

  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outputName);
  fun4allServer->registerOutputManager(out);

  RawClusterBuilderTemplate* rawClusterBuilder = new RawClusterBuilderTemplate("myClusterBuilder");
  rawClusterBuilder->Detector("CEMC");
  rawClusterBuilder->setSubclusterSplitting(false);
  rawClusterBuilder->set_UseTowerInfo(1);
  rawClusterBuilder->set_threshold_energy(0.070);
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  rawClusterBuilder->LoadProfile(emc_prof);
  rawClusterBuilder->setSubclusterSplitting(false);
  rawClusterBuilder->setOutputClusterNodeName("CLUSTERINFO_CEMC_NO_SPLIT");
  rawClusterBuilder->set_UseTowerInfo(1); // to use towerinfo objects rather than old RawTower
  fun4allServer->registerSubsystem(rawClusterBuilder);  

  // int run=0;
  // DisplayLegoPlot* module_display_lego_plot = new DisplayLegoPlot("myDisplayLego", run, pid_str, save_tree);
  // module_display_lego_plot->set_output_file( Form("/sphenix/u/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/%s.root", pid_str.c_str()));
  // fun4allServer->registerSubsystem(module_display_lego_plot);  

  fun4allServer->run(nEvents);

  fun4allServer->End();
  std::cout << "All done" << std::endl;
  delete fun4allServer;
  if (Enable::PRODUCTION) {
	  Production_MoveOutput();
  }

  gSystem->Exit(0);
  return EXIT_SUCCESS;
}

void ensure_dir(const std::string& path) {
  struct stat info;

  if (stat(path.c_str(), &info) != 0)
  {
    std::cout << "Directory " << path << " does not exist. Creating..." << std::endl;
    if (mkdir(path.c_str(), 0777) != 0)
    {
      std::cerr << "Failed to create directory " << path << std::endl;
      exit(1);
    }
  }
  else if (!(info.st_mode & S_IFDIR))
  {
    std::cerr << "Path " << path << " exists but is not a directory!" << std::endl;
    exit(1);
  }
};

