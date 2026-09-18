#include "WorkingPointScan.h"

int ReduceWorkingPointScan(const std::string family, const std::string map_root, const std::string output_base,
    const bool require_complete, const std::string sample, const int shard_index = 0, const int shard_count = 1, const Long64_t max_events = 0)
{
  using namespace tagging_mass;
  const auto definitions = samples(family);
  const auto definition = std::find_if(definitions.begin(), definitions.end(), [&](const auto& s) { return s.name == sample; });
  if (definition == definitions.end() || output_base.empty() || shard_count <= 0 || shard_index < 0 || shard_index >= shard_count || max_events < 0) return 1;
  const std::string output_path = output_base + "/working_point_scan.root";
  if (!gSystem->AccessPathName(output_path.c_str())) return 2;
  Metadata meta; meta.family = family; meta.sample = sample; meta.require_complete = require_complete;
  meta.shard_index = shard_index; meta.shard_count = shard_count; meta.max_events = max_events;
  Input input; double tag_threshold = -1;
  if (!inspect(*definition, map_root + "/" + sample + "/map_*.root", require_complete, input, meta.release, meta.model, tag_threshold, meta.selection_settings, meta.region_settings)) return 3;
  const auto& s = meta.selection_settings;
  if (s[0] > .12 || !close(s[1], .15) || !close(s[2], .5) || !close(s[3], .08) || !close(s[4], .18) || !close(s[5], .45) || !close(s[6], .65)) return 3;
  meta.map_count = input.maps; meta.sumw = input.sumw;
  TChain tree("event_tree"); if (tree.Add(input.pattern.c_str()) <= 0) return 3;
  meta.total_entries = tree.GetEntries(); meta.entry_begin = meta.total_entries * shard_index / shard_count;
  meta.entry_end = meta.total_entries * (shard_index + 1) / shard_count;
  if (max_events) meta.entry_end = std::min(meta.entry_end, meta.entry_begin + max_events);
  if (meta.entry_end <= meta.entry_begin) return 3;
  unsigned char valid=0, stitch_valid=0, stitch_pass=0;
  unsigned int region_count=0;
  double numerator=0;
  std::vector<unsigned int>* id=nullptr;
  std::vector<double> *e=nullptr, *et=nullptr, *eta=nullptr, *phi=nullptr;
  std::vector<unsigned char> *region=nullptr, *prompt=nullptr, *ptag=nullptr, *etag=nullptr, *final=nullptr;
  std::vector<float> *fraction=nullptr, *pi0_fraction=nullptr, *contributor_fraction=nullptr;
  std::vector<unsigned char>* pi0_valid=nullptr;
  std::vector<unsigned int>* offset=nullptr;
  std::vector<int> *pdg=nullptr, *embedding=nullptr, *source=nullptr, *parent=nullptr;
  tree.SetBranchStatus("*", false); tree.SetCacheSize(32LL*1024*1024);
#define B(name, var) tree.SetBranchStatus(name, true); if (!bind_branch(tree, name, &var)) return 4
  B("event_weight_valid",valid); B("sample_stitching_valid",stitch_valid); B("sample_stitching_pass",stitch_pass);
  B("region_a_count",region_count); B("weight_numerator_pb",numerator);
  B("split_cluster_id",id); B("split_cluster_e",e); B("split_cluster_et",et); B("split_cluster_eta",eta); B("split_cluster_phi",phi);
  B("split_cluster_pass_region_a",region); B("split_cluster_truth_prompt_cluster",prompt); B("split_cluster_truth_dominant_fraction",fraction);
  B("split_cluster_pi0_tag",ptag); B("split_cluster_eta_tag",etag); B("split_cluster_pass_final_photon",final);
  B("split_cluster_pi0_anchor_valid",pi0_valid); B("split_cluster_pi0_anchor_main_fraction",pi0_fraction);
  B("split_cluster_truth_contributor_offset",offset); B("split_cluster_truth_contributor_g4_pdg_id",pdg);
  B("split_cluster_truth_contributor_embedding_id",embedding); B("split_cluster_truth_contributor_fraction",contributor_fraction);
  B("split_cluster_truth_contributor_photon_source",source); B("split_cluster_truth_contributor_classification_parent_pdg",parent);
#undef B
  const std::vector<const char*> gate = {"event_weight_valid","sample_stitching_valid","sample_stitching_pass","region_a_count"};
  const std::vector<const char*> detail = {"weight_numerator_pb","split_cluster_id","split_cluster_e","split_cluster_et","split_cluster_eta","split_cluster_phi",
      "split_cluster_pass_region_a","split_cluster_truth_prompt_cluster","split_cluster_truth_dominant_fraction","split_cluster_pi0_tag","split_cluster_eta_tag","split_cluster_pass_final_photon"};
  for (auto name : gate) tree.AddBranchToCache(name,true);
  for (auto name : detail) tree.AddBranchToCache(name,true);
  for (auto name : {"split_cluster_pi0_anchor_valid", "split_cluster_pi0_anchor_main_fraction", "split_cluster_truth_contributor_offset",
      "split_cluster_truth_contributor_g4_pdg_id", "split_cluster_truth_contributor_embedding_id", "split_cluster_truth_contributor_fraction",
      "split_cluster_truth_contributor_photon_source", "split_cluster_truth_contributor_classification_parent_pdg"}) tree.AddBranchToCache(name,true);
  tree.StopCacheLearningPhase();
  working_point_scan::Histograms hist;
  unsigned long long selected = 0;
  for (Long64_t entry=meta.entry_begin; entry<meta.entry_end; ++entry)
  {
    if(tree.GetEntry(entry)<=0)return 5;
    if(!valid || !stitch_valid || !stitch_pass || region_count==0) continue;
    if(!id || !e || !et || !eta || !phi || !region || !prompt || !fraction || !ptag || !etag || !final || !pi0_valid || !pi0_fraction || !offset || !pdg || !embedding || !contributor_fraction || !source || !parent) return 5;
    const std::size_t n=id->size();
    for(auto size:{e->size(),et->size(),eta->size(),phi->size(),region->size(),prompt->size(),fraction->size(),ptag->size(),etag->size(),final->size(),pi0_valid->size(),pi0_fraction->size()}) if(size!=n)return 5;
    if (offset->size()!=n+1 || offset->front()!=0 || offset->back()!=contributor_fraction->size() ||
        !std::is_sorted(offset->begin(),offset->end())) return 5;
    for (auto size : {pdg->size(),embedding->size(),source->size(),parent->size()}) if (size!=contributor_fraction->size()) return 5;
    const double weight=numerator/input.sumw; if(!std::isfinite(weight))return 5;
    for(std::size_t i=0;i<n;++i)
    {
      if(!(*region)[i])continue;
      ++selected;
      std::array<bool,working_point_scan::npoints> tags{}; bool eta_tag=false;
      for(std::size_t j=0;j<n;++j)
      {
        if(i==j || (*id)[i]==(*id)[j] || (*e)[j]<=.12)continue;
        const double mass=pair_mass((*e)[i],(*eta)[i],(*phi)[i],(*e)[j],(*eta)[j],(*phi)[j]);
        if(!std::isfinite(mass) || mass<0)return 5;
        for(int w=0;w<working_point_scan::npoints;++w)if((*e)[j]>working_point_scan::cuts[w] && mass>working_point_scan::lower[w] && mass<working_point_scan::upper[w])tags[w]=true;
        if((*e)[j]>.5 && mass>.45 && mass<.65)eta_tag=true;
      }
      if(tags[0]!=bool((*ptag)[i]) || eta_tag!=bool((*etag)[i]) || bool((*final)[i])!=(!tags[0]&&!eta_tag))
      {std::cerr<<"Baseline veto mismatch: "<<entry<<" cluster "<<(*id)[i]<<std::endl;return 6;}
      unsigned mask=0;for(int w=0;w<working_point_scan::npoints;++w)if(!tags[w]&&!eta_tag)mask|=1U<<w;
      double eta_fraction=0;
      for (unsigned j=(*offset)[i]; j<(*offset)[i+1]; ++j)
        if ((*embedding)[j]==1 && (std::abs((*pdg)[j])==221 || ((*pdg)[j]==22 && ((*parent)[j]==221 || (*source)[j]==3))))
          eta_fraction+=(*contributor_fraction)[j];
      if (!std::isfinite(eta_fraction) || eta_fraction < -1e-6 || eta_fraction > 1+1e-5) return 5;
      const bool is_prompt=(*prompt)[i] && (*fraction)[i]>.5F, is_pi0=(*pi0_valid)[i] && (*pi0_fraction)[i]>.5F;
      const bool exclusive=int(is_prompt)+int(is_pi0)+int(eta_fraction>.5)<=1;
      hist.fill(is_prompt && exclusive,is_pi0 && exclusive,(*et)[i],mask,weight);
    }
  }
  if(gSystem->mkdir(output_base.c_str(),true)!=0 && gSystem->AccessPathName(output_base.c_str()))return 7;
  const std::string tmp=output_path+".tmp."+std::to_string(gSystem->GetPid());
  TFile out(tmp.c_str(),"CREATE");if(out.IsZombie())return 7;
  hist.write();meta.write();out.Close();
  if(out.TestBit(TFile::kWriteError) || !gSystem->AccessPathName(output_path.c_str()) || gSystem->Rename(tmp.c_str(),output_path.c_str())!=0)return 7;
  std::cout<<"Completed "<<sample<<" shard "<<shard_index<<"; events="<<meta.entry_end-meta.entry_begin<<" candidates="<<selected<<std::endl;
  return 0;
}
