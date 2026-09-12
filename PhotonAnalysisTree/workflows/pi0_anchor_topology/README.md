# Pythia pi0-main anchor-cluster topology workflow

This workflow starts from every selected central SPLIT CEMC cluster for which a
selected pi0 is the grouped main truth contributor. Each anchor cluster is
filled exactly once into separated, merged, single-contaminated, missing, or other. Missing is partitioned into six exclusive categories: energy band with mass inside/outside the pi0 window, low energy, unclustered or no CEMC deposit, acceptance, and other. Their sum equals total missing in every ET bin, including underflow and overflow.

## Event selection

Pythia events use the truth collision vertex from the signal embedding
(`PHHepMCGenEvent::get_collision_vertex()`, in cm). The anchor-spectrum
workflow requires

    std::abs(z_vertex) < max_abs_vertex_z

with default `max_abs_vertex_z = 60.0` cm. Events at or beyond the boundary,
including exactly +60 cm and -60 cm, are rejected before cluster truth matching.
Rejected events fill no prompt, anchor, candidate, or topology histogram or
counter and are recorded separately as `events_vertex_rejected`. Invalid input
events remain recorded as `events_invalid`.

The shared `Pi0AnchorTopologyEvaluator` keeps this cut disabled by default.
The anchor-spectrum workflow enables it explicitly, so `TopologyEventDisplay`
continues to evaluate all valid vertices.

## Anchor and pi0 definition

For each cluster, contributor fractions compatible with the same selected pi0
are summed. The pi0 with the largest grouped fraction is the main contributor
when that fraction is at least anchor_pi0_fraction_min (default 0.5).
Equal-leading assignments are retained in the anchor denominator but classified
as other.

Selected pi0s are either transported G4-primary pi0 decays or generator pi0
decays represented by exactly two G4-primary photons. Detector-secondary pi0s and Dalitz decays are excluded. No truth-eta
selection is applied to the parent pi0 or either daughter photon.

The anchor requires |eta_cluster| < 0.7 and the configured anchor energy cut. Partner lookup uses its independent strict `E > pi0_partner_min_energy` cut (defaulting to the anchor threshold) and no eta cut by default. Diagnostic partner lookup searches all valid positive-energy clusters.

## Energy-deposit topology

For each direct pi0 daughter photon, the matcher finds the cluster with maximum
absolute daughter energy deposit among clusters passing the respective anchor or partner
selection. The calibrated-energy estimate is retained for diagnostics:

    Erec_gamma = Ecluster * (Edep_gamma / Edep_total)

    Erec_gamma / Etruth_gamma >= min_photon_energy_recovery

The default recovery threshold is now zero, which disables the ratio requirement entirely. A positive value can still be passed for comparison productions. The optional cluster-composition requirement

    Edep_gamma / Edep_cluster > min_energy_contribution_fraction

is applied while selecting the maximum-deposit cluster; its default is 0.0, so
any positive daughter deposit is accepted. There is no cluster/pi0 response
cut.

For one anchor cluster:

- single-contaminated: it first satisfies the merged condition, and exactly one
  direct daughter photon has its first daughter-production vertex at transverse
  radius `r < pre_cemc_interaction_radius` (default 90 cm);
- merged: it is the recovered maximum-deposit cluster of both daughter photons,
  except for the single-contaminated case above. Events in which neither or both
  direct photons first produce daughters before 90 cm remain merged. The
  boundary is strict, so a first daughter vertex at exactly 90 cm is not
  pre-CEMC;
- separated: it is the recovered maximum-deposit cluster of one daughter and
  the other daughter has a distinct recovered maximum-deposit partner cluster;
- missing: it is the recovered maximum-deposit cluster of one daughter and
  the other daughter has no usable cluster passing the partner energy cut (and any explicitly enabled recovery cut).
  Missing is split with the following exclusive priority:

  1. invalid projection -> other;
  2. outside CEMC acceptance -> acceptance;
  3. usable representative truth partner -> classify its full cluster energy: `E <= L` is low energy; `L < E <= U` is split by pi0 mass window only; `E > U`, invalid pair mass, or the same cluster as anchor -> other;
  4. unresolved matching -> other;
  5. otherwise -> unclustered or no CEMC deposit.

`L=0.2`, `U=0.5 GeV`, and pi0 window `0.10 < mass < 0.20 GeV` are configurable. The representative partner maximizes direct daughter deposit independently of the production energy threshold. Displacement and match-incomplete are no longer histogram categories. Diagnostics-disabled cases fall into other unless acceptance applies; partner topology matching still runs.

If one pi0 produces extra pi0-main fragment clusters, only daughter
maximum-deposit clusters can be merged, separated, or missing; the extra
anchors are intentionally other. This preserves the exact cluster-level
partition without duplicate fills.

## Build and production

Build and install only after jobs using the current installed library have
finished:

    PhotonAnalysisTree/src/build.sh

submit.job is configured for the minimum-bias production. Review paths,
file/job counts, and parameters before submitting manually:

    condor_submit workflows/pi0_anchor_topology/submit.job

No repository script submits jobs automatically. Each job writes
transactionally and validates its partial before publication. This change uses
partial schema 9; do not mix it with earlier partial schemas. Point the next
production at a new, empty output directory. run_partial.sh
accepts optional CEMC_ACCEPTANCE_ETA_MAX,
MIN_DIRECT_MATCH_CLUSTER_ENERGY_COVERAGE, MISSING_DIAGNOSTIC_MAX_DELTA_R, and
ENABLE_MISSING_DIAGNOSTICS, and PRE_CEMC_INTERACTION_RADIUS arguments after
MAX_ABS_VERTEX_Z; their defaults are 1.1, 0.5, 0.15, true, and 90.0 cm.

Finalize a complete production with:

    root -l -b -q 'workflows/pi0_anchor_topology/FinalizePythiaPi0AnchorClusterSpectra.C("output/intermediate_files/pi0_anchor_topology_partial/eta07_zvtx60_full_partner_fgamma0p0_recovery0p5_clusterenergy/partial_*.root","output/plots/pi0_anchor_topology/minimum_bias/eta07_zvtx60_full_partner_fgamma0p0_recovery0p5_clusterenergy",0,200000,"Pythia8 p+p MB")'

The finalizer writes the combined raw and bin-width-normalized spectra, the
aggregate missing spectrum, and detailed and summary category fractions relative
to the anchor spectrum. It produces two spectrum PDFs:

- `output_base.pdf`: summary spectrum with the prompt-photon reference, all
  anchors, and the five exclusive topology categories;
- `output_base_detailed.pdf`: detailed spectrum retaining the total missing
  spectrum and its six exclusive subcategories.

The summary category plots use five exclusive categories:
separated, merged, single contaminated, missing, and other. The detailed plots
retain the six exclusive missing subcategories. Four category PDFs are
produced:

- `output_base_category_fractions.pdf`: summary line plot;
- `output_base_category_fraction_stack.pdf`: summary stacked plot;
- `output_base_category_fractions_detailed.pdf`: detailed line plot;
- `output_base_category_fraction_stack_detailed.pdf`: detailed stacked plot.

The final ROOT schema is 10 and stores both detailed and summary fraction
histograms. Plot annotations and legends are placed outside the histogram frame.

The optional tail after `PRE_CEMC_INTERACTION_RADIUS` in `run_partial.sh` is `PI0_PARTNER_MIN_ENERGY PI0_MASS_MIN PI0_MASS_MAX MISSING_ENERGY_MIN MISSING_ENERGY_MAX`. These values are exposed in `submit.job`, saved in partial metadata, and checked before finalization. Topology algorithm 11 selects the representative truth partner before applying its energy threshold; photon-energy recovery remains disabled by default. Commit b7eb4d6 retains the previous configurable recovery behavior. A positive recovery threshold remains available for explicit comparison runs; zero bypasses the fraction check entirely.
