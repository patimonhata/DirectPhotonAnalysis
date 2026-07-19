#ifndef RYOTARO_Pi0Reconstruction_H_20260210
#define RYOTARO_Pi0Reconstruction_H_20260210

#include <fun4all/SubsysReco.h>

#include <TH1D.h>

#include <array>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;

class Pi0Reconstruction : public SubsysReco {
  public:
    enum class VertexMode {
      Origin,
      GlobalVertexMap
    };

    explicit Pi0Reconstruction(const std::string &name = "Pi0Reconstruction");
    ~Pi0Reconstruction() override;

    void set_output_file_name(const std::string &output_file_name);
    void set_cluster_node_name(const std::string &cluster_node_name);
    void set_process_id(unsigned int process_id);
    void set_vertex_node_name(const std::string &vertex_node_name);
    void set_vertex_mode(VertexMode vertex_mode);
    void set_abort_on_missing_cluster_node(bool abort_on_missing_cluster_node);
    void set_abort_on_missing_vertex_node(bool abort_on_missing_vertex_node);
    void set_min_cluster_energy(double min_cluster_energy);
    void set_mass_histogram_bins(int nbins, double min, double max);

    int Init(PHCompositeNode *topNode) override;
    int InitRun(PHCompositeNode *topNode) override;
    int process_event(PHCompositeNode *topNode) override;
    int ResetEvent(PHCompositeNode *topNode) override;
    int Reset(PHCompositeNode *topNode) override;
    int EndRun(const int runnumber) override;
    int End(PHCompositeNode *topNode) override;

  private:
    struct PhotonCandidate {
      double energy = 0.0;
      std::array<double, 3> momentum = {0.0, 0.0, 0.0};
    };

    bool get_event_vertex(PHCompositeNode *topNode, std::array<double, 3> &vertex);
    bool build_photon_candidate(double energy, double x, double y, double z, const std::array<double, 3> &vertex, PhotonCandidate &candidate) const;
    void reset_tree_variables();
    void create_tree_branches();
    void create_output_directory() const;

    std::string output_file_name_ = "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/pi0_reconstruction.root";
    std::string cluster_node_name_ = "CLUSTER_CEMC";
    std::string vertex_node_name_ = "GlobalVertexMap";
    VertexMode vertex_mode_ = VertexMode::Origin;
    unsigned int process_id_ = 0;

    bool abort_on_missing_cluster_node_ = true;
    bool abort_on_missing_vertex_node_ = false;
    double min_cluster_energy_ = 0.0;

    int mass_histogram_nbins_ = 100;
    double mass_histogram_min_ = 0.0;
    double mass_histogram_max_ = 1.0;

    TFile *output_file_ = nullptr;
    TTree *event_tree_ = nullptr;
    TH1D *h_m_gg_ = nullptr;
    TH1D *h_ncluster_ = nullptr;
    TH1D *h_cluster_e_ = nullptr;
    TH1D *h_pair_e_asym_ = nullptr;

    unsigned int tree_process_id_ = 0;
    unsigned int tree_event_ = 0;
    unsigned long long tree_event_uid_ = 0;
    unsigned int tree_ncluster_ = 0;
    unsigned int tree_ncluster_all_ = 0;
    double tree_min_cluster_energy_ = 0.0;
    double tree_vertex_x_ = 0.0;
    double tree_vertex_y_ = 0.0;
    double tree_vertex_z_ = 0.0;

    std::vector<double> tree_cluster_e_;
    std::vector<double> tree_cluster_x_;
    std::vector<double> tree_cluster_y_;
    std::vector<double> tree_cluster_z_;
    std::vector<double> tree_cluster_px_;
    std::vector<double> tree_cluster_py_;
    std::vector<double> tree_cluster_pz_;
    std::vector<unsigned int> tree_pair_cluster_i_;
    std::vector<unsigned int> tree_pair_cluster_j_;
    std::vector<double> tree_pair_m_gg_;
    std::vector<double> tree_pair_e_asym_;

    unsigned int event_counter_ = 0;
    unsigned int missing_cluster_node_warnings_ = 0;
    unsigned int missing_vertex_node_warnings_ = 0;
};

#endif // RYOTARO_Pi0Reconstruction_H_20260210
