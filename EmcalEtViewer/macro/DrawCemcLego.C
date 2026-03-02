// ROOT macro: interactive CEMC E_T lego plot
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iostream>
#include <string>
#include <vector>

void DrawCemcLego(const char *file_name,
                  const char *tree_name = "cemc_et",
                  const char *mode = "vtx")
{
  std::cout << "Came to 1" << std::endl;

  TFile *file = TFile::Open(file_name, "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Failed to open file: " << file_name << std::endl;
    return;
  }
  std::cout << "Came to 2" << std::endl;

  TTree *tree = dynamic_cast<TTree *>(file->Get(tree_name));
  if (!tree) {
    std::cout << "Failed to find TTree: " << tree_name << std::endl;
    return;
  }
  std::cout << "Came to 3" << std::endl;

  int event = 0;
  float vtx_x = 0.0f;
  float vtx_y = 0.0f;
  float vtx_z = 0.0f;

  std::vector<float> *eta0 = nullptr;
  std::vector<float> *etavtx = nullptr;
  std::vector<float> *phi0 = nullptr;
  std::vector<float> *phivtx = nullptr;
  std::vector<float> *et0 = nullptr;
  std::vector<float> *etvtx = nullptr;

  tree->SetBranchStatus("*", 0);  // いったん全部OFF

  tree->SetBranchStatus("event", 1);
  tree->SetBranchStatus("vtx_x", 1);
  tree->SetBranchStatus("vtx_y", 1);
  tree->SetBranchStatus("vtx_z", 1);
  // // tree->SetBranchStatus("eta0", 1);
  // tree->SetBranchStatus("etavtx", 1);
  // // tree->SetBranchStatus("phi0", 1);
  // tree->SetBranchStatus("phivtx", 1);
  // // tree->SetBranchStatus("et0", 1);
  // tree->SetBranchStatus("etvtx", 1);
  
  const bool use_vtx = (std::string(mode) == "vtx");

  if (use_vtx) {
    tree->SetBranchStatus("etavtx", 1);
    tree->SetBranchStatus("phivtx", 1);
    tree->SetBranchStatus("etvtx", 1);
  } else {
    tree->SetBranchStatus("eta0", 1);
    tree->SetBranchStatus("phi0", 1);
    tree->SetBranchStatus("et0", 1);
  }

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("vtx_x", &vtx_x);
  tree->SetBranchAddress("vtx_y", &vtx_y);
  tree->SetBranchAddress("vtx_z", &vtx_z);

  if (use_vtx) {
    tree->SetBranchAddress("etavtx", &etavtx);
    tree->SetBranchAddress("phivtx", &phivtx);
    tree->SetBranchAddress("etvtx", &etvtx);
  } else {  
    tree->SetBranchAddress("eta0", &eta0);
    tree->SetBranchAddress("phi0", &phi0);
    tree->SetBranchAddress("et0", &et0);
  }

  std::cout << "Came to 4" << std::endl;

  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("cemc_lego", "CEMC E_{T} lego", 1200, 800);
  TH2F *h = new TH2F("h_et", ";#eta;#phi;E_{T} [GeV]",
                     96, -1.2, 1.2,
                     256, -TMath::Pi(), TMath::Pi());

  const Long64_t n_entries = tree->GetEntries();
  std::cout << "Came to 5" << std::endl;
  for (Long64_t ievt = 0; ievt < n_entries; ++ievt) {
    tree->GetEntry(ievt);
    h->Reset("ICES");

    const std::vector<float> &eta = use_vtx ? *etavtx : *eta0;
    const std::vector<float> &phi = use_vtx ? *phivtx : *phi0;
    const std::vector<float> &et = use_vtx ? *etvtx : *et0;

    for (size_t i = 0; i < et.size(); ++i) {
      h->Fill(eta[i], phi[i], et[i]);
    }

    h->SetTitle(Form("Event %d  vtx=(%.2f, %.2f, %.2f) cm  mode=%s", event, vtx_x, vtx_y, vtx_z, mode));
    h->Draw("lego2");
    canvas->Update();
    // canvas->SaveAs("event.pdf");
    canvas->SaveAs("event.png");

    /* Update the canvas upon some primitive action on the gPad */
    // std::cout << "Event " << event << " displayed. Click the canvas or press any key inside it to advance." << std::endl;
    // gPad->WaitPrimitive();
    
    /* or update upon a keyboard input */
    std::cout << "Event " << event << " displayed. Press Enter to advance (q + Enter to quit): " << std::endl;
    std::string line;
    std::getline(std::cin, line);
    if (line == "q" || line == "Q") break;
    
    /* or update in a constant seconds */
    // gSystem->Sleep(1000); // 1000 ms = 1秒
    
    gSystem->ProcessEvents();
  }
}
