"""Synthetic schema-8 map -> reduce -> merge, including eta cascade modes, pi0 priority and distinct-parent rejection."""
import os
from pathlib import Path
import re
import tempfile
import ROOT

repo = Path(__file__).resolve().parents[2]
os.chdir(repo)
ROOT.gROOT.SetBatch(True)
assert ROOT.gROOT.LoadMacro('PhotonAnalysisTree/workflows/photon_candidate_selection/MergePythiaPhotonCandidateSelection.C') == 0
ROOT.gInterpreter.Declare(r'''
#include <numeric>
void make_eta_test_map(const std::string& path, const std::string& sample_name) {
  const auto samples = sample_definitions("photonjet");
  const auto sample = *std::find_if(samples.begin(), samples.end(), [&](const auto& s) { return sample_name == s.name; });
  TFile file(path.c_str(), "RECREATE");
  MapMetadata map;
  map.schema_version = 8; map.decay_origin_algorithm_version = 1; map.decay_origin_max_proper_length_cm = 1e-4; map.decay_origin_max_flight_cm = 1e-2; map.pi0_topology_algorithm_version = 12; map.eta_topology_algorithm_version = 2;
  map.input_manifest = "synthetic"; map.sample_name = sample.name; map.analysis_release = "ana.565"; map.model_sha256 = "synthetic";
  map.manifest_begin = 0; map.manifest_end = 1; map.input_file_count = 1;
  map.cross_section_pb = sample.cross_section_pb; map.window_min = sample.window_min; map.window_max = sample.window_max;
  map.upper_unbounded = sample.upper_unbounded; map.sum_generator_weight_processed = 23;
  map.events_processed = map.events_written = 23;
  map.min_cluster_energy = .2; map.partner_diagnostic_min_cluster_energy = 0;
  map.pi0_partner_min_energy = map.eta_partner_min_energy = .5;
  map.pi0_mass_min = .1; map.pi0_mass_max = .2; map.eta_mass_min = .45; map.eta_mass_max = .65;
  map.missing_energy_min = .2; map.missing_energy_max = .5; map.min_photon_energy_recovery = 0;
  TTree metadata("metadata", "Synthetic map");
#define M(name, member) metadata.Branch(name, &map.member)
  M("decay_origin_algorithm_version", decay_origin_algorithm_version); M("decay_origin_max_proper_length_cm", decay_origin_max_proper_length_cm);
  M("decay_origin_max_flight_cm", decay_origin_max_flight_cm); M("schema_version", schema_version); M("input_manifest", input_manifest); M("sample_name", sample_name);
  M("analysis_release", analysis_release); M("model_sha256", model_sha256); M("manifest_begin", manifest_begin); M("manifest_end", manifest_end);
  M("input_file_count", input_file_count); M("map_chunk_id", map_chunk_id); M("sample_cross_section_pb", cross_section_pb);
  M("sample_window_min", window_min); M("sample_window_max", window_max); M("sample_upper_unbounded", upper_unbounded);
  M("sum_generator_weight_processed", sum_generator_weight_processed); M("n_events_processed", events_processed);
  M("n_events_written", events_written); M("n_events_vertex_rejected", events_vertex_rejected); M("n_events_invalid", events_invalid);
  M("min_cluster_energy", min_cluster_energy); M("partner_diagnostic_min_cluster_energy", partner_diagnostic_min_cluster_energy);
  M("pi0_partner_min_energy", pi0_partner_min_energy); M("eta_partner_min_energy", eta_partner_min_energy);
  M("pi0_mass_min", pi0_mass_min); M("pi0_mass_max", pi0_mass_max); M("eta_mass_min", eta_mass_min); M("eta_mass_max", eta_mass_max);
  M("missing_energy_min", missing_energy_min); M("missing_energy_max", missing_energy_max); M("min_photon_energy_recovery", min_photon_energy_recovery);
  M("pi0_topology_algorithm_version", pi0_topology_algorithm_version); M("eta_topology_algorithm_version", eta_topology_algorithm_version);
#undef M
  metadata.Fill(); metadata.Write();
  TTree events("event_tree", "Synthetic event");
  unsigned int ncluster = 1;
  unsigned char enabled = 1;
  double weight = sample.cross_section_pb;
  std::vector<double> et = {7};
  std::vector<float> score = {.9}, dominant = {.8}, pi0_fraction = {-.999}, eta_fraction = {-.999}, fraction = {.8};
  std::vector<unsigned char> yes = {1}, no = {0}, isolated = {1}, nonisolated = {0}, region = {1}, final = {1};
  std::vector<unsigned char> prompt = {0}, pi0_valid = {0}, eta_valid = {0}, pi0_tag = {0}, eta_tag = {0};
  std::vector<int> pi0_topology = {-999}, eta_topology = {-999}, pi0_missing = {-999}, eta_missing = {-999};
  std::vector<unsigned int> offsets = {0, 1};
  std::vector<int> pdg = {22}, embedding = {1}, source = {0}, parent = {0}, photon_category = {0};
#define E(name, value) events.Branch(name, &value)
  E("split_ncluster", ncluster); E("event_weight_valid", enabled); E("sample_stitching_valid", enabled); E("sample_stitching_pass", enabled);
  E("weight_numerator_pb", weight); E("split_cluster_et", et); E("split_cluster_bdt_score", score); E("split_cluster_bdt_valid", yes);
  E("split_cluster_pass_kinematics", yes); E("split_cluster_pass_preselection", yes); E("split_cluster_pass_tight", yes);
  E("split_cluster_pass_isolated", isolated); E("split_cluster_pass_nonisolated", nonisolated); E("split_cluster_pass_region_a", region);
  E("split_cluster_pass_final_photon", final); E("split_cluster_pi0_tag", pi0_tag); E("split_cluster_eta_tag", eta_tag);
  E("split_cluster_truth_valid", yes); E("split_cluster_truth_prompt_cluster", prompt); E("split_cluster_truth_dominant_fraction", dominant);
  E("split_cluster_pi0_anchor_valid", pi0_valid); E("split_cluster_pi0_anchor_main_fraction", pi0_fraction);
  E("split_cluster_pi0_anchor_topology", pi0_topology); E("split_cluster_pi0_anchor_missing_category", pi0_missing);
  E("split_cluster_eta_anchor_valid", eta_valid); E("split_cluster_eta_anchor_main_fraction", eta_fraction);
  E("split_cluster_eta_anchor_topology", eta_topology); E("split_cluster_eta_anchor_missing_category", eta_missing);
  E("split_cluster_truth_contributor_offset", offsets); E("split_cluster_truth_contributor_g4_pdg_id", pdg);
  E("split_cluster_truth_contributor_embedding_id", embedding); E("split_cluster_truth_contributor_fraction", fraction);
  E("split_cluster_truth_contributor_photon_category", photon_category);
  E("split_cluster_truth_contributor_photon_source", source); E("split_cluster_truth_contributor_classification_parent_pdg", parent);
  std::vector<int> eta_parent_barcode = {0}; E("split_cluster_truth_contributor_eta_parent_barcode", eta_parent_barcode);
  std::vector<int> eta_parent_track_id = {0}; E("split_cluster_truth_contributor_eta_parent_track_id", eta_parent_track_id);
  std::vector<int> eta_decay_mode = {0}; E("split_cluster_truth_contributor_eta_decay_mode", eta_decay_mode);
  std::vector<int> eta_split_status = {0}; E("split_cluster_truth_contributor_eta_split_status", eta_split_status);
  std::vector<int> pi0_parent_barcode = {0}; E("split_cluster_truth_contributor_pi0_parent_barcode", pi0_parent_barcode);
  std::vector<int> pi0_parent_track_id = {0}; E("split_cluster_truth_contributor_pi0_parent_track_id", pi0_parent_track_id);
#undef E
  constexpr int values[] = {1, 2, 4, 3, 0};
  std::vector<unsigned char> d_valid = {1}; std::vector<unsigned int> d_offset;
  events.Branch("split_cluster_decay_origin_valid", &d_valid); events.Branch("split_cluster_decay_origin_offset", &d_offset);
  std::vector<int> d_pdg; events.Branch("split_cluster_decay_origin_pdg", &d_pdg);
  std::vector<int> d_embedding; events.Branch("split_cluster_decay_origin_embedding", &d_embedding);
  std::vector<int> d_source; events.Branch("split_cluster_decay_origin_source", &d_source);
  std::vector<int> d_id; events.Branch("split_cluster_decay_origin_id", &d_id);
  std::vector<int> d_root_track_id; events.Branch("split_cluster_decay_origin_root_track_id", &d_root_track_id);
  std::vector<int> d_parent_pdg; events.Branch("split_cluster_decay_origin_parent_pdg", &d_parent_pdg);
  std::vector<int> d_parent_source; events.Branch("split_cluster_decay_origin_parent_source", &d_parent_source);
  std::vector<int> d_parent_id; events.Branch("split_cluster_decay_origin_parent_id", &d_parent_id);
  std::vector<int> d_status; events.Branch("split_cluster_decay_origin_status", &d_status);
  std::vector<int> d_pathway; events.Branch("split_cluster_decay_origin_pathway", &d_pathway);
  std::vector<int> d_diphoton; events.Branch("split_cluster_decay_origin_diphoton", &d_diphoton);
  std::vector<int> d_long_lived; events.Branch("split_cluster_decay_origin_long_lived", &d_long_lived);
  std::vector<float> d_fraction; events.Branch("split_cluster_decay_origin_fraction", &d_fraction);
  std::vector<float> d_gamma0_fraction; events.Branch("split_cluster_decay_origin_gamma0_fraction", &d_gamma0_fraction);
  std::vector<float> d_gamma1_fraction; events.Branch("split_cluster_decay_origin_gamma1_fraction", &d_gamma1_fraction);
  std::vector<float> d_flight_cm; events.Branch("split_cluster_decay_origin_flight_cm", &d_flight_cm);
  std::vector<float> d_parent_flight_cm; events.Branch("split_cluster_decay_origin_parent_flight_cm", &d_parent_flight_cm);
  std::vector<float> d_proper_length_cm; events.Branch("split_cluster_decay_origin_proper_length_cm", &d_proper_length_cm);
  std::vector<int> d_daughter_count, d_daughter_pi0_count, d_daughter_gamma_count, d_daughter_charged_count;
  events.Branch("split_cluster_decay_origin_daughter_count", &d_daughter_count);
  events.Branch("split_cluster_decay_origin_daughter_pi0_count", &d_daughter_pi0_count);
  events.Branch("split_cluster_decay_origin_daughter_gamma_count", &d_daughter_gamma_count);
  events.Branch("split_cluster_decay_origin_daughter_charged_count", &d_daughter_charged_count);
  for (int row = 0; row < 23; ++row) {
    offsets = {0, 1}; pdg = {22}; embedding = {1}; fraction = {.8}; source = {0}; parent = {0};
    prompt[0] = pi0_valid[0] = eta_valid[0] = pi0_tag[0] = eta_tag[0] = 0;
    pi0_fraction[0] = eta_fraction[0] = -999; pi0_topology[0] = eta_topology[0] = pi0_missing[0] = eta_missing[0] = -999;
    isolated[0] = row != 8 && row != 9; nonisolated[0] = row == 8; region[0] = isolated[0];
    if (row == 0) prompt[0] = 1;
    else if (row <= 5) { pi0_valid[0] = 1; pi0_fraction[0] = .8; pi0_topology[0] = values[row-1]; pi0_missing[0] = row == 4 ? 1 : 0; parent[0] = 111; }
    else if (row <= 15) {
      eta_valid[0] = 1; eta_fraction[0] = .8; eta_topology[0] = row <= 10 ? values[row-6] : 3;
      eta_missing[0] = row == 9 ? 1 : row >= 11 ? row - 9 : 0; parent[0] = 221;
      eta_tag[0] = row == 6; pi0_tag[0] = row == 7;
    } else if (row == 16) { offsets = {0, 2}; pdg = {22, 22}; embedding = {1, 1}; fraction = {.3, .3}; source = {3, 3}; parent = {221, 221}; }
    else if (row == 17) { eta_valid[0] = 1; eta_fraction[0] = .5; fraction[0] = .5; eta_topology[0] = 0; eta_missing[0] = 0; parent[0] = 221; }
    else if (row == 18) { pdg[0] = 221; }
    else { pdg[0] = 22; parent[0] = 223; }
    eta_parent_barcode.assign(fraction.size(), 0); eta_parent_track_id.assign(fraction.size(), 0);
    eta_decay_mode.assign(fraction.size(), 0); eta_split_status.assign(fraction.size(), 0);
    pi0_parent_barcode.assign(fraction.size(), 0); pi0_parent_track_id.assign(fraction.size(), 0);
    for (std::size_t i = 0; i < fraction.size(); ++i) if (parent[i] == 221 || pdg[i] == 221)
    { eta_parent_barcode[i] = row == 16 ? i + 100 : 100; eta_decay_mode[i] = pdg[i] == 221 ? 2 : 1; }
    if (row >= 20)
    {
      offsets = {0, 3}; pdg = {211, -211, 111}; embedding = {1, 1, 1}; fraction = {.3, .3, .4};
      source = {0, 0, 0}; parent = {0, 0, 0}; eta_parent_barcode = {100, 100, 100}; eta_parent_track_id = {0, 0, 0};
      eta_decay_mode = {3, 3, 3}; eta_split_status = {0, 0, 0}; pi0_parent_barcode = {0, 0, 111}; pi0_parent_track_id = {0, 0, 0};
      if (row == 21)
      {
        pdg = {111, 111, 111}; fraction = {.2, .2, .6}; eta_decay_mode = {2, 2, 2}; pi0_parent_barcode = {111, 112, 113};
        pi0_valid[0] = 1; pi0_fraction[0] = .6; pi0_topology[0] = 2; pi0_missing[0] = 0;
      }
      if (row == 22) { fraction = {.6, .2, .2}; eta_parent_barcode = {0, 0, 0}; eta_decay_mode = {0, 0, 0}; }
    }
    photon_category.clear();
    for (std::size_t i = 0; i < pdg.size(); ++i) photon_category.push_back(pdg[i] == 22 ? (std::abs(parent[i]) > 37 ? 3 : 0) : -1);
    d_offset = offsets; d_pdg = pdg; d_embedding = embedding; d_fraction = fraction; d_parent_pdg = parent;
    const auto size = fraction.size(); d_source.assign(size, 2); d_id.resize(size); std::iota(d_id.begin(), d_id.end(), 1);
    d_root_track_id = d_id; d_parent_source.assign(size, 0); d_parent_id.assign(size, 0);
    d_daughter_count.assign(size, 0); d_daughter_pi0_count.assign(size, 0); d_daughter_gamma_count.assign(size, 0); d_daughter_charged_count.assign(size, 0);
    d_status.assign(size, 0); d_pathway.assign(size, 0); d_diphoton.assign(size, 0); d_long_lived.assign(size, 0);
    d_gamma0_fraction.assign(size, -1); d_gamma1_fraction.assign(size, -1); d_flight_cm.assign(size, -1);
    d_parent_flight_cm.assign(size, -1); d_proper_length_cm.assign(size, -1);
    if (row == 19) { d_pdg[0] = 111; d_diphoton[0] = 1; d_gamma0_fraction[0] = .5; d_gamma1_fraction[0] = .3; d_pathway[0] = 2; d_status[0] = 1; }
    final[0] = region[0] && !pi0_tag[0] && !eta_tag[0]; events.Fill();
  }
  events.Write();
}
''')

with tempfile.TemporaryDirectory(prefix='eta_pipeline_') as temporary:
    root = Path(temporary)
    for sample in ('photonjet3', 'photonjet5', 'photonjet10', 'photonjet20'):
        maps = root / 'maps' / sample
        maps.mkdir(parents=True)
        ROOT.make_eta_test_map(str(maps / 'map_000000.root'), sample)
        partial = root / 'partial' / sample / 'shard_0' / 'photon_candidate_selection'
        assert ROOT.ReducePythiaPhotonCandidateSelection('photonjet', str(root / 'maps'), str(partial), False, 40, 40., sample, 0) == 0
        f = ROOT.TFile.Open(str(partial) + '.root')
        assert f.Get('composition/kinematic/h_candidate_eta_et_count').Integral() == 12
        assert f.Get('eta_anchor_topology/kinematic/h_region_a_eta_anchor_et_count').Integral() == 10
        assert f.Get('eta_anchor_topology/region_a_tagging_veto/h_region_a_eta_anchor_et_count').Integral() == 6
        assert f.Get('anchor_topology/kinematic/h_region_a_pi0_anchor_et_count').Integral() == 6
        for name in ('separated', 'merged', 'single_contaminated', 'missing', 'other'):
            assert f.Get(f'anchor_topology/kinematic/h_region_a_{name}_et_count').Integral() == (2 if name == 'merged' else 1)
        assert f.Get('eta_decay/kinematic/h_candidate_eta_three_pi0_et_count').Integral() == 1
        assert f.Get('eta_decay/kinematic/h_candidate_eta_charged_pions_pi0_et_count').Integral() == 1
        assert f.Get('composition/kinematic/h_candidate_other_et_count').Integral() == 4
        assert f.Get('other_origin/kinematic/h_other_decay_photon_et_count').Integral() == 1
        assert f.Get('other_origin/kinematic/h_other_charged_hadron_et_count').Integral() == 1
        assert f.Get('other_origin/kinematic/h_other_mixed_et_count').Integral() == 2
        f.Close()
    output = root / 'merged'
    assert ROOT.MergePythiaPhotonCandidateSelection('photonjet', str(root / 'partial'), str(output)) == 0
    f = ROOT.TFile.Open(str(output / 'candidate_composition_eta_topology.root'))
    assert f.Get('kinematic/detailed/h_candidate_eta_residual_et_count').Integral() == 0
    assert f.Get('kinematic/detailed/h_candidate_eta_three_pi0_et_count').Integral() == 4
    assert f.Get('kinematic/detailed/h_candidate_eta_charged_pions_pi0_et_count').Integral() == 4
    f.Close()
    import csv
    census = list(csv.DictReader((output / 'eta_cascade_composition_summary.csv').open()))
    other = next(row for row in census if row['selection'] == 'kinematic' and row['bin'] == '-1' and row['category'] == 'h_candidate_other_et_count')
    assert float(other['raw_count']) == 16 and float(other['total_raw_count']) == 92
    detailed = ROOT.TFile.Open(str(output / 'other_origin.root'))
    assert detailed.Get('kinematic/h_other_decay_photon_et_count').Integral() == 4
    assert detailed.Get('kinematic/h_other_charged_hadron_et_count').Integral() == 4
    assert detailed.Get('kinematic/h_other_mixed_et_count').Integral() == 8
    assert detailed.Get('kinematic/h_other_decay_photon_main_fraction_vs_et_pb').GetEntries() == 4
    detailed.Close()
    other_rows = list(csv.DictReader((output / 'other_origin_summary.csv').open()))
    row = next(x for x in other_rows if x['selection'] == 'kinematic' and x['bin'] == '-1' and x['category'] == 'decay_photon')
    assert float(row['raw_count']) == 4 and float(row['other_raw_count']) == 16 and float(row['all_raw_count']) == 92
    assert (output / 'preselection/photon_candidate_other_origin_category_fraction_stack.pdf').is_file()
    assert (output / 'preselection/photon_candidate_composition_other_origin_category_fraction_stack_detailed.pdf').is_file()
    diagnostic = ROOT.TFile.Open(str(output / 'decay_origin.root'))
    assert diagnostic.Get('kinematic/h_decay_origin_pi0_gg_et_count').Integral() == 4
    assert diagnostic.Get('kinematic/h_decay_origin_charged_et_count').Integral() == 4
    assert diagnostic.Get('kinematic/h_decay_origin_mixed_et_count').Integral() == 8
    diagnostic.Close()
    diagnostic_rows = list(csv.DictReader((output / 'decay_origin_summary.csv').open()))
    row = next(x for x in diagnostic_rows if x['selection'] == 'kinematic' and x['bin'] == '-1' and x['category'] == 'pi0_gg')
    assert float(row['raw_count']) == 4 and float(row['other_raw_count']) == 16
    assert (output / 'region_a_tagging_veto/photon_candidate_other_decay_origin_category_fraction_stack.pdf').is_file()
    f = ROOT.TFile.Open(str(output / 'pi0_anchor_topology.root'))
    assert f.Get('kinematic/h_region_a_pi0_anchor_et_count').Integral() == 24
    f.Close()
    for directory in (output / 'kinematic', output / 'region_a_tagging_veto'):
        for parent in ('pi0', 'eta'):
            assert (directory / f'region_a_{parent}_anchor_topology.pdf').is_file()
            for suffix in ('', '_isolation'):
                pdf = directory / f'{parent}_bdt_score_by_topology{suffix}.pdf'
                assert len(re.findall(rb'/Type\s*/Page\b', pdf.read_bytes())) == 6
    print('Schema-8 map -> reduce -> merge with supplemental decay origins, eta residuals, veto survival, both plot families: passed')
