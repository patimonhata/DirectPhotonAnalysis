# Photon-candidate definition notes

Development format: map schema 10, partial schema 13, common composition v1, pi0 topology v14, eta topology v4, direct hit matcher v5. Saved schema-8 productions retain their original definitions and results. A new production requires a new clean-source software snapshot. Map-affecting changes require regenerated maps; reduce-only corrections still follow production provenance checks. Section headings below retain the version in which a feature was introduced. Use [the workflow README](../workflows/photon_candidate_selection/README.md) for current paths and execution; historical command blocks were removed.

## BDT model and input order

The development map uses `models/model_ppg15_nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_single_tmva.root`, with SHA-256
`5494946f40213a1035dc9e883bc71ddfdcebf86015d279300c8f44774e846248`.
`TMVA::Experimental::RBDT` loads key `myBDT`. `run_map.sh` passes the model from
`PHOTON_SOURCE_ROOT/models/` explicitly; the direct map macro defaults to the same
model.

For the retained `newsignal_v1` production, `models/model_tmva_newsignaldefinition.root` has SHA-256
`4cf8a46e3f25d9d68e499f2a2bf395beb053d0dbc613070edf763c4d18f70e86` and is byte-identical to the completed **nominal** acceptance model at
`/sphenix/tg/tg01/coldqcd/jaein213/photon/bdt_v11_acceptance_20261102/models/nominal/model_tmva.root`.
Its SHA-256 also matches that model's `COMPLETE.json` artifact record.
The training/export source for that model is
`/gpfs/mnt/gpfs02/sphenix/user/jaein213/photon/BDT/PPG15PhotonAN/bdt_training/studies/v11_acceptance_20261102/`:
`vendor/v11_photonlike_fullsample_20261011/common.py::FEATURES`, `export_wide.py`,
`train.py`, and the output `resolved_config.yaml` confirm the 11 positions below.
The export uses `SaveXGBoost(..., 'myBDT', ..., num_inputs=11)` without a feature scaler.
The `model_sha256` metadata value is a fixed code value, not a runtime file hash.
Keep the model file, configured path and recorded SHA-256 together when changing it.

`PythiaPhotonCandidateTree.cc` passes these 11 inputs in this exact order:

| Index | Training feature | Map input |
|---|---|---|
| 0 | `cluster_Et` | cluster ET |
| 1 | `cluster_weta_cogx` | `shower_w_eta_cogx` |
| 2 | `cluster_wphi_cogx` | `shower_w_phi_cogx` |
| 3 | `vertexz` (reconstructed MBD z) | signal HepMC collision z (`b_vertex_z_`) |
| 4 | `cluster_Eta` | cluster eta |
| 5 | `e11_over_e33` | `shower_e11_over_e33` |
| 6 | `cluster_et1` | `shower_et1` |
| 7 | `cluster_et2` | `shower_et2` |
| 8 | `cluster_et3` | `shower_et3` |
| 9 | `cluster_et4` | `shower_et4` |
| 10 | `e32_over_e35` | `shower_e32_over_e35` |

The input positions agree with the `newsignal_v1` training/export order. That training uses
reconstructed MBD z, while `fill_event_truth()` assigns the map's
`b_vertex_z_` from `PHHepMCGenEvent::get_collision_vertex()`. The current map also
uses this signal vertex for cluster kinematics. The model-only comparison retains
these existing map definitions and the tight/non-tight formulas in Selection
definitions below. A model switch changes
stored scores and score-dependent flags, so regenerate maps from DST in a new
production before reduce/merge. Replotting stored histograms cannot apply a new model.

## Boundary independent composition (map schema 9)

This section supersedes the historical origin/population rules below. The composition change itself did not alter reconstructed selections, scores, models, training labels or the reconstructed tag veto. The current model is specified separately above.

`CommonDecayTruth` resolves an event-local graph across HepMC/G4 using `(event_uid, embedding, source, barcode/track ID)`. The shared physical particle uses its unambiguous matching HepMC identity when available. G4 daughter vertices must be recorded decays; generator decays have one incoming hadron. Beam particles, partons, hadronization and interactions never define parent sums. Photon bookkeeping copies are collapsed. Shower electrons/photons remain under the original independent photon; they are not extra primary photons.

The exclusive units are individual pi0 (including non-gg decays), direct eta-to-gg, independent photons outside those units, charged particles and residual particles. Contributions of distinct particles are never combined merely because their PDGs agree. First select a unit with deposit **strictly greater than half** the original truth deposit. This includes charged particles and direct decay photons. Only if no unit has a majority, sum confirmed descendants of each individual hadron and choose the downstream-most majority parent. Without any majority, use mixed when the required information is known; otherwise preserve unknown and its reason. Signal embedding is 1. A non-signal majority is separate from unknown. Prompt requires the existing HepMC category 1/2; an unclassified photon is not implicitly prompt.

The denominator is the original sum of shower deposits from `PythiaClusterTruthMatcher`, not reconstructed cluster energy. `original_cluster_truth_*` retains the original primary contributions. Common expansion uses complete hit ancestry with the existing split-tower allocation and measured conditional daughter ratios. A split is accepted only when child deposits conserve the measured parent deposit; otherwise retain the known parent with unresolved daughter shares. Parent summary rows overlap the unit table and must never be added into the deposit total. Legacy eta refinement feeds only the retained legacy contributor table and migration reference, never the common contribution input.

There is **no lifetime or flight-distance expansion cut** in either representation. Proper lengths, known production/decay vertices, flight lengths, and long-lived ancestry are diagnostics. Unknown flight alone does not invalidate known composition. The `long_lived_ancestor` diagnostic uses c-tau >= 1 micrometre independently of classification; this threshold cannot change main categories.

Known parent origin with unresolved daughters (`parent_unresolved`) is separate from confirmed parent majority and mixed. A known pi0 majority remains pi0 when other deposits are unresolved. Individual photon deposits, largest/other independent-photon deposits, charged deposits and unresolved deposits are retained for future purity work. These fields do not change BDT training or determine clean-photon labels.

### Saved common tables and quality

`split_cluster_common_*` stores main, selected identity/PDG/fraction, original total, attributed total, unresolved deposit, independent-photon purity components, quality flags, topology applicability/reason, matching status/minimum coverage, and legacy overview origin (0 Prompt, 1 pi0, 2 eta, 3 Other). Main codes follow `composition::Main` in `CommonTruthComposition.h`: unknown=0, prompt=1, pi0=2, eta_gg=3, decay_photon=4, charged=5, particle=6, parent=7, parent_unresolved=8, mixed=9, non_signal=10, photon_unknown=11.

`offset` indexes particle rows per cluster. `table` is 0 units, 1 independent photons, 2 overlapping parent totals. Rows retain particle IDs, deposit, mode, prompt category, unresolved flag, lifetime/flight and vertices. Known independent photon IDs are retained even when their shares are unresolved: table-1 deposit -1 is an unavailable value, never a negative physical contribution or zero deposit. Maximum/other-photon deposits sum only resolved rows and therefore are lower bounds when unresolved photons exist. The cluster unresolved deposit includes both unresolved classification units and unresolved photon shares within known pi0/eta units; row-level flags distinguish these scopes. `ancestor_offset` indexes each row's full confirmed ancestry; ancestor IDs retain embedding/source as well as PDG. Empty vertex arrays and negative flight/proper lengths mean unavailable. Quality bits are missing particle=1, ambiguous identity=2, incomplete split=4, interaction=8, unknown species=16, ancestry cycle=32, invalid/unavailable truth deposit or denominator=64, unknown decay=128.

### Composition and topology populations

Pi0 composition equals topology-evaluable + topology-not-applicable + topology-unknown. A confirmed direct gg decay with identifiable photons and adequate matching of that anchor is evaluable. Confirmed non-gg pi0 is not applicable; unknown decay mode, identity, necessary projection, or insufficient matching is unknown. Eta overview includes direct eta-to-gg, eta parent majority and eta origin with unresolved daughter shares. Direct eta topology and the remaining eta population are recorded separately.

`applicability`: 0 unknown, 1 evaluable, 2 not applicable. `topology_reason`: 0 available, 1 confirmed non-gg/parent sum, 2 unknown decay, 3 unavailable candidate/photon identity, 4 insufficient matching, 5 unavailable projection, 6 unavailable first-daughter interaction information. Mapping statuses remain invalid=0, partial=1, complete=2. Topology accepts complete or the existing partial matching with member-energy coverage >= 0.5; daughter splitting requires complete conserved attribution. Coverage is distinct from photon fraction and from reconstructed energy purity. Missing hit ancestry makes matching unusable. Noise energy is never assigned to a fabricated particle. Candidate-wide minimum coverage and matching status remain saved, with ET and reason in `matching_count`/`matching_weighted`, as diagnostics only. Starting with map schema 10, one failed searched cluster or sibling anchor no longer vetoes all anchors of that meson. Applicability requires the current anchor's existing `match.usable`, both photons' valid projection or known non-intersection, and first-daughter interaction information. Partner selection and `partner_direct_match_incomplete` diagnostics are unchanged; no partner representative is required merely to evaluate a missing topology. The selected maximum is therefore a maximum over available matching information, not proof that incomplete clusters contain no larger deposit.

`split_cluster_common_anchor_matching_status`, `anchor_matching_coverage`, and `anchor_matching_usable` separately save the current anchor's match. Before anchor evaluation they are -1, -1, and false. Evaluated status uses 0/1/2 as above. `anchor_matching_count`/`anchor_matching_weighted` preserve ET, local status, local coverage and topology reason, only for evaluated anchors. Candidate aggregates retain their original meaning even when the anchor is evaluable. Map schema 9 and partial schema 11 must not be mixed with map schema 10 / partial schema 12 or 13; regenerate maps from DST for the new applicability policy.

Topology fractions and topology BDT distributions use only the evaluable population. `other` still means an evaluable anchor outside the named topology cases; it does not absorb missing information. The single-contaminated definition remains exactly one photon with first daughter-vertex radius < 90 cm; both photons interacting falls under merged. Partner selection maximizes absolute deposit over positive-energy clusters before applying the actual configured partner threshold. A known ray with no forward intersection with the CEMC reference cylinder is saved separately as `photon*_projection_no_intersection`; it is geometric non-arrival (missing/acceptance when unrecovered), not an unavailable-vertex unknown. Ties choose the smallest cluster ID. A failing representative is never replaced by a lower-deposit candidate. Reconstructed veto flags remain independent.

### Outputs, denominators and migration

All six selections reuse the same stored main class. Fractions count each selected cluster once with its event weight and reconstructed ET. Summary composition remains Prompt/pi0/eta/Other. Detailed common composition includes pi0 and direct-eta topology, unavailable/not-applicable topology, eta/other parent majority versus unresolved daughters, decay photons, charged/residual particles, mixed, unknown and non-signal.

`common_composition.root`, `common_composition_summary.csv`, and `common_composition_migration.csv` cover **all** origins, selections, ET bins and flows. Count, weighted yield and Sumw2 close for exclusive details, legacy-to-common transitions and the new topology population partitions. Sparse tables preserve ancestry species and matching quality versus ET. Migration-reason bits record possible boundary generalization (1), charged/direct-photon priority (2), long-lived ancestry (4), explicit topology quality (8), and missing information (16); these overlapping diagnostic flags are not a causal decomposition. Exact event-local particle identities remain in maps.

Composition fractions divide by all selected clusters; topology fractions divide by evaluable mesons. Survival divides by the same new class at kinematic or the specified predecessor: kinematic for preselection; preselection for tight/isolation/Region A; Region A for final veto. BDT shapes normalize each ET/topology/isolation distribution independently. ET pages remain [5,6), [6,8), [8,10), [10,18), [18,30), and the full selected range; 30–35 appears only in the full-range page.

Old Other-origin and supplemental decay-origin outputs retain the **legacy Other population**, with their old lifetime/flight policy. They are historical diagnostics, not tests of the new formal classification. No `workinprogress` directory or duplicate WIP plots are created. Existing WIP outputs are retained.

## Region A/B/C/D content

Every stored split cluster has these independent flag branches:

- `split_cluster_pass_region_a`: isolated and tight
- `split_cluster_pass_region_b`: non-isolated and tight
- `split_cluster_pass_region_c`: isolated and non-tight
- `split_cluster_pass_region_d`: non-isolated and non-tight
- `split_cluster_pass_final_photon`: Region A after the pi0-or-eta tag veto

The event tree also stores `region_a_count`, `region_b_count`, `region_c_count`, `region_d_count`, and `final_photon_count`. The one-entry metadata tree stores their map-file totals as `n_clusters_region_a/b/c/d` and `n_clusters_final_photon`. Clusters in the isolation gap or BDT gap have no A/B/C/D flag, which is intentional.

Thus the future ABCD purity calculation can be performed entirely from these intermediate TTrees. No Region-B/C/D cluster is discarded at the map stage.


## Selection definitions

- The signal HepMC collision vertex must satisfy `|z| < 60 cm`; missing or rejected vertices are not written.
- `TOPOCLUSTER_ALLCALO` is reconstructed in the job with the requested EMCal+HCal topological-clustering configuration.
- Only split clusters from `CLUSTERINFO_CEMC` with `E > min_cluster_energy` are stored; the threshold is strict and configured per map production.
- Shower shapes use the existing 7x7 calculator and a 70 MeV per-tower threshold. Tower patches, constituent towers, and the all-pairs table are not populated.
- Raw isolation is `sum(TopoCluster ET, deltaR < 0.4) - candidate ET`.
- Corrected simulation isolation is `1.2 * raw + 0.1 GeV`.
- Isolated means `corrected_iso < 0.490 + 0.037 * candidate_ET`.
- Non-isolated means `corrected_iso > 0.490 + 0.037 * candidate_ET + 0.8 GeV`.
- The common ABCD candidate requirement is `5 < ET < 35 GeV`, `|eta| < 0.7`, and the shower-shape preselection.
- Tight means `score > 0.8156 - 0.00156 ET`.
- Non-tight is strictly `0.7333 - 0.01333 ET < score < 0.6844 + 0.00156 ET`.

Meson tagging considers every other valid reconstructed split cluster in `CLUSTERINFO_CEMC`, without a partner eta cut. Pi0 and eta each require their own strict `E > threshold` and `mass_min < mass < mass_max`. Any qualifying pair sets the tag; among multiple pairs, the saved veto partner is the pair closest to the nominal meson mass. Truth matching is not used in this veto. `split_cluster_pass_final_photon` remains `region_a && !pi0_tag && !eta_tag`.


## Partner and missing definitions

The anchor selection keeps `min_cluster_energy`. Partner lookup first selects the usable truth match with the largest absolute direct daughter deposit across all positive-energy clusters. This single representative is used for both topology and saved truth-pair diagnostics. It must then satisfy strict `E > pi0_partner_min_energy` to count as recovered: a distinct partner gives separated, while an anchor-identical partner enters merged/single-contaminated classification. If it fails the energy requirement, the anchor enters missing; no lower-deposit alternative cluster is substituted. Matching-quality conditions and the disabled minimum recovered-energy fraction are unchanged. The calibrated photon energy estimate and recovery ratio remain diagnostics.

For missing diagnostics, all valid positive-energy clusters are searched, including those at or below 0.1 GeV. The representative truth partner maximizes direct daughter deposit among usable truth matches, independently of the tagging threshold. Its full reconstructed cluster energy (not the truth-photon energy or the attributed fraction) determines the missing energy category. The reconstruction's own clustering thresholds remain unchanged.

The six exclusive missing categories are:

| Code | Key | Definition |
|---|---|---|
| 1 | `energy_band_taggable` | `L < E <= U`, representative pair inside pi0 mass window |
| 4 | `energy_band_not_taggable` | `L < E <= U`, finite representative pair mass outside pi0 window |
| 5 | `low_energy` | `E <= L` |
| 6 | `unclustered_or_no_cemc_deposit` | No associated cluster found, or no CEMC deposit |
| 2 | `acceptance` | Partner projection outside CEMC acceptance |
| 3 | `other` | Remaining cases, including unresolved matching, invalid mass, same cluster as anchor, and `E > U` |

`L=0.2` and `U=0.5 GeV` by default. The mass-only taggability in this table asks whether lowering the partner energy threshold could tag the pair; it does not impose the production tagging threshold. Invalid projection goes to other; acceptance takes priority over energy classification. Incomplete matching is not labelled unclustered when it prevents identifying a partner. Displacement no longer defines a missing category. Legacy detailed reason codes remain diagnostic information, not additional histogram categories.

The independent diagnostic axes remain:

- `partner_alignment`: 0 not applicable, 1 near, 2 displaced (`deltaR > 0.15`), 3 projection invalid, 4 cluster unavailable;
- `truth_partner_tag_status`: 0 not applicable, 1 taggable, 2 unavailable, 3 same as anchor, 4 below the production pi0 tagging threshold, 5 mass outside pi0 window, 6 invalid mass;
- `tag_result`: 0 not applicable, 1 survived, 2 truth pair taggable veto, 3 combinatorial-only veto.

These axes use the deposit-selected truth partner, which need not be the reconstructed partner selected by the actual veto.


## Configurable parameters

Map submit files define `pi0_partner_min_energy`, `eta_partner_min_energy`, `pi0_mass_min`, `pi0_mass_max`, `eta_mass_min`, `eta_mass_max`, `missing_energy_min`, and `missing_energy_max` (GeV). Initial pi0/eta thresholds preserve each submit file's old common threshold; the mass windows default to 0.10/0.20 and 0.45/0.65. Missing boundaries default to 0.2/0.5. Bounds must be finite, non-negative, and ordered.

The optional tail of `run_map.sh` after `MIN_CLUSTER_ENERGY_GEV` is:

`run_small_sample.sh` and `check_sample_map_outputs.sh` accept the same tail after their stored-cluster threshold / existing pi0 threshold position. Omitted pi0 energy defaults to the stored-cluster threshold, omitted eta energy to pi0, and omitted windows/bounds to the values above. All settings are preserved in ROOT metadata and compared when reducing and merging. Use a separate output directory for every configuration.


## Stitching and normalization

PhotonJet windows use the leading terminal prompt HepMC photon pT (classifier category 1 or 2). Jet windows use the leading `AntiKt_Truth_r04` jet pT. Events are retained regardless of the window and carry `sample_stitching_valid` and `sample_stitching_pass`; downstream aggregation must require both.

For each complete sample, the reduce stage must calculate

where `weight_numerator_pb = sample_cross_section_pb * generator_weight`. The denominator includes events rejected later by the vertex cut. Do not normalize map files independently. Cross sections and half-open stitching windows are stored in metadata; `jet40` has no upper bound.


## Pi0 BDT distributions and common population

Topology spectra and BDT diagnostics use the topology-evaluable subset of common pi0 composition. Pi0 topology-not-applicable and unknown remain in composition and are excluded from topology/BDT denominators. Reduce and merge check count, weighted yield, Sumw2 and BDT totals under this inclusion relation. New maps must be regenerated from DST; reduce/merge alone cannot migrate schema-8 results. Use a new software/production manifest and output directory.

Every selection directory (`kinematic`, `preselection`, `preselection_isolation`, `preselection_tight`, `region_a`, `region_a_tagging_veto`) receives:

- `pi0_bdt_score_by_topology.pdf`: six pages, two columns (weighted yield with log axis; independently unit-normalized topology shapes).
- `pi0_bdt_score_by_topology_isolation.pdf`: six pages, three columns (all, isolated, nonisolated) and two rows (weighted yield; normalized shapes), including selections that already require isolation.
- `pi0_bdt_score.root`: count/weighted score histograms, weighted score-versus-ET histograms, selection/invalid-score counters, and population/selection/model metadata.
- `pi0_bdt_score_summary.csv`: selected and invalid-score counts and weights, valid-score weighted integrals, and score underflow/overflow counts for each topology/isolation group.

The five BDT categories are separated, merged, single contaminated, missing and other. No missing subdivisions are added to BDT diagnostics. Existing missing-detail plots remain available. All includes the isolation gap; isolated and nonisolated use stored independent flags. Gap is saved for closure checks. Invalid scores are counted separately. Score uses 100 bins over [0,1]; exactly 1 is included in the last visible bin. True underflow/overflow is retained and reported; unit-area plots normalize visible bins only. Empty distributions are not normalized. Isolation-required stages retain empty nonisolated panels with an explanatory label. Tight selections intentionally retain their BDT cut.

Run `python3 PhotonAnalysisTree/tests/test_pi0_bdt.py` in the ROOT environment for synthetic classification, partition, merge and output checks.


### ET pages and redraw from existing BDT ROOT files

Both BDT PDFs use pages 1--5 for reconstructed anchor ET in [5,6), [6,8), [8,10), [10,18), [18,30) GeV, respectively, and page 6 for the original full selection (5 < ET < 35 GeV). The existing kinematic cut excludes exactly 5 GeV. The 30--35 GeV population appears only on page 6. Each slice projects the stored 0.5-GeV-bin score-versus-ET histogram, preserving weighted errors. Shape normalization is independent for each ET/isolation/topology combination. Empty pages/panels are retained.

The saved ROOT histograms and CSV continue to contain the original full-ET counters and 2D weighted distributions. Since unweighted and invalid-BDT counters were not stored by ET, slice pages show the ET range instead of those counters; the full-ET page retains the original counter annotations. No reduce rerun or schema change is required.

To update only the two BDT PDFs in all six selection directories of an existing merged output, load `ReplotPi0BdtScores.C` and call `ReplotPi0BdtScores("/absolute/path/to/merged/jet")` in ROOT. This reads each `pi0_bdt_score.root` without modifying it or the CSV. Future ordinary merges automatically produce the six-page PDFs as well.



## Direct eta-to-two-photon topology extension

Schema-6 maps add `eta_topology_*` candidate/daughter/anchor records and
`split_cluster_eta_anchor_*` cluster links. The shared evaluator is run separately
for PDG 221. Generator candidates must have exactly two photon daughters in the
HepMC decay vertex and exactly two corresponding G4-primary photons. Unique photon
bookkeeping copies are followed. G4-primary eta candidates must have exactly two
direct G4 photon daughters. Eta-to-pi0 cascades and all other eta decays are excluded.
When this extension was introduced in map schema 6, the pi0 candidate selection and algorithm version 11 were unchanged. Current algorithm versions are listed at the start of this document.

Eta topology uses the existing `eta_partner_min_energy` and `eta_mass_min/max`
settings, including for missing energy-band mass classification and truth-pair
taggability. Anchor energy, eta acceptance, missing-energy bands, deposit matching,
and pre-CEMC interaction definitions are shared with pi0. Recovery remains disabled.

Legacy composition is unchanged: `eta_fraction()` still sums compatible eta-origin
contributors without grouping by parent. New eta topology and BDT histograms use a
subset of that existing composition eta category: a valid direct eta-to-two-photon
anchor with the **same-parent** contribution strictly greater than 0.5. An exact
0.5 is excluded. Multiple eta parents whose sum exceeds 0.5, unsupported decays,
and other eta composition entries without such an anchor remain in the residual
`Other eta contributions` category in the new composition plots. That residual is
formed by subtracting the anchor subset from legacy eta; its Sumw2 is subtracted
as well, respecting subset covariance. The five eta topology categories partition
all selected eta anchors, and the six missing categories partition missing.

All six existing selection directories gain these separate outputs:

- `region_a_eta_anchor_topology*.pdf`: eight spectrum, fraction, stack and
  kinematic-relative survival plots with the same layout and normalization as pi0.
- `eta_bdt_score_by_topology.pdf` and `eta_bdt_score_by_topology_isolation.pdf`:
  the same six ET pages, isolation groups and invalid-score accounting as pi0.
- `eta_bdt_score.root` and `eta_bdt_score_summary.csv`.
- `photon_candidate_composition_eta_topology_category_fraction_stack_detailed.pdf`:
  prompt, aggregate pi0, five eta topology categories, residual eta and other origin.
- `photon_candidate_composition_eta_topology_category_fraction_stack_superdetailed.pdf`:
  the same partition with eta missing expanded into its six categories.
- `photon_candidate_composition_eta_topology_survival_fraction_relative_to_kinematic_detailed.pdf`:
  the survival of the detailed components relative to their own kinematic population.

The additional spectra, fractions, uncertainties and source metadata are saved in
`eta_anchor_topology.root` and `candidate_composition_eta_topology.root` (schema 1).
Existing pi0 and composition file names, histogram categories and plot layouts are
retained. That eta extension did not retrain or change the BDT; the current model
is specified above. `ReplotEtaBdtScores.C` redraws only
eta BDT PDFs from their stored ROOT histograms.

Regenerate maps from DST, then reduce and merge. Old schema-5 maps and schema-6
partials cannot supply the new eta records. Current submit files default to
`production_directory = UNCONFIGURED`; initialize a new production manifest and
pass its configuration directory explicitly. The small-sample QA wrapper defaults
to `output/photon_candidate_selection/qa/common_composition_v1/` with threshold,
sample and segment-count subdirectories. No script automatically submits jobs.
The pi0-only spectrum workflow and event display are outside this extension.

Map production selects the headers, libraries, model and source snapshot bound to
the production manifest. Commit changes, build a development prefix, freeze a new
software ID and initialize a new production; existing frozen software remains intact.
See [software/build provenance](software.md). The tests include HepMC direct
eta decay and photon-copy ancestry, rejection of extra daughters and pi0 cascades,
eta mass-window boundaries, ROOT branch roundtrips, category/BDT closure, residual
covariance, and a synthetic four-sample map-to-reduce-to-merge run that verifies both
pi0 and eta plot families.


## Eta decay flight diagnostic (independent of origin classification)

`CheckEtaDecayFlight.C` reads existing G4Hits DSTs directly. It does not load the
analysis library, regenerate maps, or change composition/topology/origin selection.
It records one row per G4-primary eta in the requested embedding (default 1),
regardless of reconstruction, cluster multiplicity, collision-z or stitching cuts.

After sourcing the sPHENIX environment, run from the repository root:

Arguments, in order:

1. Input G4Hits list: one local path or FROG-resolvable DST name per line; blank
   lines and lines beginning with `#` are ignored. Paths relative to the list are
   supported. Duplicate inputs are rejected.
2. Output base, without extension. Existing ROOT/PDF/CSV/overview files are never
   overwritten; choose a new base for another run.
3. First file index (default 0).
4. Number of files (default 1; 0 means all remaining files).
5. Maximum events per file (default 0 means all). Inclusive mode reads the first N
   entries in each selected file; use a list with dispersed file indices to avoid
   treating this as an unbiased random sample.
6. Optional request CSV (default empty), selecting specific eta tracks rather than
   an inclusive population. Required columns are `source_file_id,event_in_file,eta_track_id`;
   `g4_eta_track` is accepted as an alias. `source_file_id` is the final segment
   number in the input DST filename, not the position in the list. Optional
   `sample,eta_barcode,hepmc_event_number` columns prevent provenance mismatches.
   `sample` uses names such as `jet40`; it is matched to `pythia8_Jet40` in the
   filename. Prior `eta_residual/dst_decay_modes.csv` is accepted directly.
   Duplicate track requests count once. Unmatched selected requests fail the run
   when no event limit is set; file ranges may intentionally exclude other requests.
7. Embedding ID (default 1).
8. Positive sample normalization factor (default 1), multiplied by the first
   HepMC event weight when available. Missing/empty HepMC weights fall back to one
   and are flagged per eta. Negative or nonfinite event weights remain in raw
   counts but are excluded from weighted probability/CDF calculations.

For example, recheck the previously audited residual eta decays for Jet40:

The eta production vertex comes from `eta.get_vtx_id()`. Direct daughters are
selected by `daughter.get_parent_id() == eta.get_track_id()` and grouped by their
production vertices. Different vertex IDs are coalesced only when their position,
time and process are exactly equal; merely nearby vertices are not merged.
A single `kPDecay` vertex identifies a confirmed decay. If process information is
unavailable, a single vertex with a recognized gamma-gamma, three-pi0,
pi-plus/pi-minus/pi0 or pi-plus/pi-minus/gamma final state is marked **inferred**.
Explicit non-decay interactions are never reinterpreted as decays. Multiple decay
vertices, missing vertices/daughters, and decays accompanied by other vertices stay
in the diagnostic census and are excluded from the usable straight-flight sample.

`length_cm` is the three-dimensional **displacement from the eta production
vertex to the decay vertex**, not the transverse radius relative to the beam and
not an accumulated track path length. `transverse_cm`, the displacement projected
parallel/perpendicular to the initial eta momentum, and its alignment cosine are
also stored. Vertex times and `delta_t` use raw `PHG4VtxPoint::get_t()` units; no
lifetime/unit conversion is assumed. Negative delta-t cases are excluded; missing
times are counted separately while finite geometric displacement remains usable.

Outputs:

- `<base>.root`: `eta_tree` with eta/daughter IDs, production/decay coordinates,
  process, mode, status, displacement and weight flags; `event_tree` with eta,
  usable/unresolved and long-displacement counts; one-entry `metadata` with schema,
  scope, input paths, settings and `completed=true` after successful output.
- `<base>_summary.csv`: raw and weighted short-distance fractions, raw median,
  95th/99th percentile and maximum, plus all diagnostic status counts. Groups
  include all eta, each decay mode, and eta-pT/energy bands 0--5, 5--10, 10--20,
  20--40, 40--80 and >=80 GeV.
- `<base>_cumulative.csv`: strict `L < threshold` CDF points and subset-covariance
  errors (binomial for raw counts), including exact zeros in their denominators.
- `<base>.pdf`: distance, transverse distance, CDF and census pages per group.
- `<base>_overview.png`: the all-eta overview page.

ROOT also stores distance-versus-eta-pT/energy and alignment histograms, per-group
status histograms and both raw/weighted CDFs. Distance logarithms cover -14 to +2
in cm; positive values outside that range are retained in under/overflow and in
the raw tree/CSV statistics. Exact zero displacement and shared production/decay
vertex IDs are counted explicitly, never replaced by a made-up positive distance.
Zero in saved truth does not establish a physically zero lifetime or a measured
sub-resolution flight distribution. Groups with no usable eta have `nan` CSV
fractions/quantiles and empty ROOT CDFs. The 1 um, 10 um, 0.1 mm and 1 mm thresholds
are reporting points only and are **not origin cuts**. Unresolved eta must not be
counted as short-flight eta when interpreting event-level counts.

An inclusive run and a requested/preselected run have different populations;
compare them as separate checks. For multiple production samples, inspect raw
counts per sample and use the appropriate cross-section/full-sample normalization
for weighted comparisons. First establish validity with this diagnostic; later
pseudo-primary daughter-pion origin classification remains a separate change.

Synthetic tests can run without rebuilding or installing PhotonAnalysisTree:

They cover translation-invariant displacement, process-confirmed/inferred decays,
zero/same-ID vertices, missing/ambiguous vertices, interactions, invalid coordinates,
negative times, strict CDF boundaries and weighted errors, plus split ROOT input,
request deduplication, output protection and missing-request failure handling.


## Historical eta daughter pseudo-primary composition (map schema 7 / topology v12, eta v2)

Candidate maps now resolve G4-primary eta decays to three pi0 or pi+ pi- pi0 into
individual daughter-pion contributions. The same ancestry information is recorded
when these decays happen in HepMC, including the mixed case where HepMC produces a
pi0 which subsequently decays in Geant4. No eta-flight cut is applied.

For Geant4 decays, direct eta daughters must have an explicit `kPDecay` vertex and
the exact exclusive daughter PDG set. Cluster tower/cell/hit ancestry measures each
pion's deposit. The daughter shares must close to the eta's hit deposit, and all
cluster towers must have complete matching information. These conditional shares
partition the **existing eta shower deposit**: total cluster truth energy and all
non-eta contributions are preserved. No generator-energy weighting is used.
Missing/partial hit ancestry or failure of the daughter sum retains the original
eta contributor, flagged unresolved; it is never assumed to be a resolved mode.

A pi0 daughter with a verified two-photon G4 decay is a topology candidate, keyed
by its G4 track ID rather than a possibly nonunique secondary barcode. Its two
photons use the existing direct-photon matching and topology logic. At reduce,
a single pi0 parent with contribution strictly greater than 0.5 takes priority.
Otherwise, descendants are grouped by **one eta parent**, and that eta must also
contribute strictly greater than 0.5. Distinct eta parents are never summed.
The eta two-photon topology remains restricted to a direct two-photon decay.
Residual eta contributions are separated into resolved three-pi0, resolved charged
pions plus pi0, unresolved daughter shares, and other eta contributions. Prompt
classification and the reconstructed photon selection/tagging cuts are unchanged.

Consequently, pi0 shares 30%/30%/40% from one eta are classified as eta->3pi0, while
20%/20%/60% is assigned to the 60% pi0's topology. A charged-pion mode is treated
by exactly the same priority rule. A sum of 30% from each of two different eta
parents gives no eta majority and remains Other unless another origin qualifies.

New flat truth-contributor branches retain `eta_parent_barcode`,
`eta_parent_track_id`, `eta_decay_mode` (0 unknown, 1 gamma gamma, 2 three pi0,
3 pi+ pi- pi0, 4 other), `eta_split_status` (1 resolved, -1 unresolved, 0 not needed),
`pi0_parent_barcode`, and `pi0_parent_track_id`. IDs are scoped to each event and
embedding; HepMC eta barcodes are used when available, with G4 eta track IDs as a
fallback. Per-cluster truth branches also store `pi0_total_fraction`,
`eta_main_fraction`, `eta_main_decay_mode`, `eta_main_barcode`,
`eta_main_track_id`, and `eta_main_unresolved` for the signal embedding (1).
The pi0 total is a diagnostic and is never substituted for a single-parent majority.

Reduce saves the disjoint eta-mode subsets under `eta_decay/<selection>`. In partial
schema 13, only common-composition `parent` with selected PDG 221 can enter the
resolved cascade subsets. Contributor annotations must match that selected parent's
embedding and identity (HepMC barcode for source 1, G4 track ID for source 2), and
all matching annotations must agree on mode 2 or 3. No legacy-largest-parent or
legacy-fraction threshold is used. Missing/conflicting modes stay in the existing
eta residual. Common `parent_unresolved` with PDG 221 enters only the unresolved
subset; `eta_gg` never enters any of these subsets, even with unknown/N/A topology.
This corrects the earlier claim of disjointness: schema 12 mixed common selection
with the legacy largest eta's mode and could count a cluster twice. All partials
must be regenerated for schema 13; schema 12 partials are rejected by merge. The
map schema remains 10, and this correction requires no new map branches. Normal
production provenance checks still apply when reusing maps with corrected software.
Merge checks raw-count and Sumw2 closure against the eta composition and direct eta
anchors. Existing eta-topology composition stack/survival PDFs now resolve these
subsets. `eta_cascade_composition_summary.csv` reports raw and weighted yields and
fractions per selection and ET bin; bin -1 is the full range including flows.
This makes the remaining Other fraction available without reading a plot.

Rebuild/install the library with `src/build.sh`, then rerun **map, reduce, merge**
into a new output directory. Old maps lack the daughter shares and are rejected.
Submit files and small-sample defaults now use `topology_v12_eta_v2`, preserving
prior v11/eta-v1 output files. The synthetic cascade test covers both generation
paths, both pion modes, 30/30/40 sharing, split-tower allocation, closure, missing
provenance, interaction rejection and grouping by distinct eta parents.


## Legacy Other-origin diagnostics (introduced with map schema 7, partial schema 9)

`OtherOriginHistograms.h` adds a disjoint subdivision of **existing composition
Other only**. Prompt, pi0, eta, their topology populations, and all reconstructed
selection/tagging decisions keep their existing definitions. No model or training
labels are changed. Commit `a8121c1` supports this subdivision using existing
v12/eta-v2 maps, with reduce and merge only. The current supplemental decay-origin
extension below requires new schema-8 maps and schema-10 partials; it retains the
same Other-origin subdivision as a separate output.

Summary plots separate long-lived charged particles, other-hadron decay-photon
majority, and the remaining Other population. Detailed plots distinguish:

- Charged hadrons: one signal-embedding contributor with abs(PDG) 211, 321 or 2212
  and fraction strictly greater than 0.5 (pi+/pi-, K+/K-, proton/antiproton).
- Electrons/positrons and muons, separately, with the same single-contributor rule.
- Other-hadron decay photon: one signal G4-primary gamma with fraction >0.5 whose
  physical HepMC parent is unique (`photon_category == 3`) and is a ROOT PDG
  Meson/Baryon (including charm/beauty classes),
  excluding pi0 and eta. Photon bookkeeping copies are already skipped in the map.
- Remaining particle majority, including neutral kaons, D/B, pi0 without a diphoton
  anchor, and photons not meeting the above parent definition.
- No single-contributor majority (an exact fraction of 0.5 belongs here).
- Unresolved truth: invalid/empty truth, a non-signal dominant contributor, or an
  existing origin overlap. These are never treated as identified charged/photon
  populations.

The three charged subcategories form the summary long-lived charged category.
Contributions from different charged particles are not summed to invent a single
particle majority. This is an origin diagnostic using stored contributors, not a
new claim that a charged track reached the EMCAL surface.

A decay-photon majority is **not yet a clean single-photon training label**: an
omega daughter gamma can dominate a cluster containing another particle's deposit.
The ROOT diagnostics retain main-contributor fraction versus ET, other G4-primary
photon fraction versus ET, BDT score versus ET, leading PDG and decay-parent PDG.
The additional-photon diagnostic does not include photons hidden inside a G4-primary
pi0/hadron; the main-fraction diagnostic includes their contamination in its
complement. No purity threshold is silently imposed. A G4-primary omega/other
hadron is retained under remaining particle majority; extracting its daughter
photon shares would require a separate map diagnostic, not an assumption in reduce.

Merge writes `other_origin.root` and `other_origin_summary.csv`. Each selection
gets these separate PDFs (with `_detailed` versions):

- `photon_candidate_other_origin_category_fraction_stack.pdf`: normalized to Other.
- `photon_candidate_composition_other_origin_category_fraction_stack.pdf`:
  normalized to all selected candidates, keeping prompt/pi0/eta aggregated.
- `photon_candidate_other_origin_survival_fraction_relative_to_kinematic.pdf`:
  category-conditional survival relative to the kinematic selection.

Detailed PDFs containing `#mu` embed their fonts using ROOT's `EmbedFonts` option
(requires Ghostscript, `gs`). The default unembedded Symbol font can display mu as
a proportional sign in some viewers. This affects rendering only; labels and data
are unchanged. Existing PDFs can be repaired by embedding fonts without rerunning
map/reduce.

Axes explicitly distinguish the two fraction denominators. CSV bin -1 reports
full-range yields including under/overflow; per-bin rows also give weighted
fractions relative to both Other and all candidates. Every ET bin (including flows)
must reproduce the existing Other raw count, weighted yield and Sumw2 before output.
Histogram axes/types and algorithm metadata are checked when partials are merged.

The original Other-only production used `.../jet/partial_other_v1` and merge
configuration suffix `_v12_eta_v2_other_v1`. Current submit files derive map and
partial directories from the explicitly selected `production_directory`.
The synthetic tests cover both
particle signs, hadron versus parton/diquark photon parents, exact-half boundaries,
non-signal/invalid truth, photon contamination, negative-weight closure and ROOT
roundtrips, plus map-to-reduce-to-merge output and unchanged pi0/eta populations.



## Legacy supplemental decay origins (introduced in map schema 8, partial schema 10)

`src/DecayOriginDiagnostics.{h,cc}` adds independent per-cluster truth information.
It never modifies the legacy truth contributors, prompt flags, pi0/eta anchors,
selection decisions, model scores or existing topology classifications. Diagnostics
are saved for every stored cluster; reduce/merge plots use **existing composition
Other only**, so the new categories form a disjoint partition of that population.

The same policy applies to hadron decay chains on both sides of the HepMC/G4
boundary. ROOT PDG lifetime, or hbar*c/width when needed, must imply proper decay
length **< 1e-4 cm (1 micrometre)**, and the recorded lab displacement must be
**<= 1e-2 cm (0.1 mm)**. These are explicit diagnostic policy parameters, not
reconstructed selection cuts or newly trained labels. Both are stored and checked
in map, partial and merged diagnostic metadata. Unknown lifetimes or missing decay
vertices are not interpreted as zero. HepMC lengths are converted from the event's
length unit; G4 vertex lengths are in cm.

Eligible short-lived hadrons are expanded recursively, stopping at an individual
pi0, a direct eta->gamma gamma, a photon, a charged particle that is not an eligible
short-lived hadron, or a long-lived/displaced parent. For example, omega->pi+ pi-
pi0 is resolved into pi+, pi- and pi0; omega->pi0 gamma into pi0 and the direct
gamma. G4 vertices must explicitly be decays; interactions are never unfolded as
prompt decays. A HepMC decay must have one incoming physical hadron; beams and
partonic/hadronization ancestry are not treated as hadron decay parents.

G4-primary particles are joined to matching HepMC identities by embedding and
barcode, with PDG agreement and unique primary-barcode mapping required. Only one
representation of a decay is expanded. Photon bookkeeping copies are skipped.
For a HepMC pi0 represented by multiple G4-primary photons, contributions are
recombined by the *same parent identity*, not by PDG. The same rule applies to
eta->gamma gamma. Long-lived/unknown ancestors remain barriers even if their decay
was already performed in HepMC: a K_S->2pi0->4gamma chain is not silently promoted
to an origin-vertex pi0. K_S remains the diagnostic parent, with its flight length
and direct daughter counts saved.

G4 daughter shares are measured using the existing cluster/tower/cell/hit ancestry
matcher and the same split-tower allocation. The sum of daughter deposits must
close to the measured parent deposit. Those conditional shares partition the
original contributor fraction; the original shower normalization and total are
preserved. Missing hits/ancestry, material interactions, or failed closure retain
the parent with an unresolved status. No daughter energy is inferred from momentum
fractions or branching ratios. Original contributor objects are never mutated.

The flat table has prefix `split_cluster_decay_origin_`:

- `valid` has one entry per cluster and carries input-truth validity and fraction
  closure; daughter resolution is given by each row's `status`. `offset` has ncluster+1 entries.
- `pdg`, `embedding`, `source`, `id` identify each diagnostic contributor. Source 1
  means an event/embedding-scoped HepMC barcode; source 2 a G4 track ID. Distinct
  meson parents remain distinct. `root_track_id` retains the input track, or 0 if
  multiple input contributors were combined.
- `fraction`, `parent_pdg`, `parent_source`, `parent_id`, `flight_cm`,
  `parent_flight_cm`, `proper_length_cm` retain contribution and physical ancestry.
  Unknown lengths are -1; stable PDG entries use a finite 1e30 cm proper-length sentinel. Parent identity is scoped to the contributor embedding.
- `daughter_count`, `daughter_pi0_count`, `daughter_gamma_count`,
  `daughter_charged_count` retain the direct decay content even for stopped parents.
- `diphoton` requires exactly two direct photons from a pi0/eta. `gamma0_fraction`
  and `gamma1_fraction` use the same cluster truth denominator. They sum to the
  meson's fraction when resolved; both are -1 if photon shares cannot be measured.
  Photon indices use sorted canonical identities, not leading/subleading energy.
- `status`: 0 retained/grouped identity, 1 contains a measured daughter split,
  -1 missing/incomplete provenance or deposit closure, -2 non-decay/inconsistent G4
  interaction, -3 unavailable lifetime/flight/PDG information. When grouping rows,
  any unresolved status is retained. `pathway` is a bit mask: 1 HepMC regrouping or
  expansion, 2 G4 regrouping or expansion. `long_lived` flags a decay retained by
  the proper-length or lab-displacement policy.

`decay_origin.root` contains counts, weighted ET yields, main fraction and BDT
score versus ET, main/parent PDG, status/pathway, log10 flight distributions, and
`second_gamma_fraction_vs_et_pb` (the smaller of the two meson photon contributions).
Unknown flights/photon shares go to underflow; an exactly zero flight is recorded
in the first log-flight bin. These photon shares do **not** implement a new global
separated/missing/merged topology or establish a clean single-photon label.

The disjoint Other categories are pi0->gamma gamma majority, eta->gamma gamma
majority, other-hadron decay-gamma majority, long-lived charged majority,
long-lived/displaced parent majority, remaining particle majority, no single
particle majority, and unresolved. Majority is strictly >0.5 for one signal-embedding
particle. Charged includes pi+/pi-, K+/K-, p/antip, electrons and muons. If no
majority exists and unresolved contributors might change that decision, the result
is unresolved rather than known mixed. Raw yields, weighted yields and Sumw2 must
close to existing Other for every ET bin, including flows.

Each selection receives two extra PDFs:

- `photon_candidate_other_decay_origin_category_fraction_stack.pdf`
- `photon_candidate_other_decay_origin_survival_fraction_relative_to_kinematic.pdf`

`decay_origin_summary.csv` gives raw and weighted yields and fractions relative to
both Other and all selected candidates. Bin -1 is the full range including flows.
The existing `other_origin.*` and all prior plot families are also still produced.
