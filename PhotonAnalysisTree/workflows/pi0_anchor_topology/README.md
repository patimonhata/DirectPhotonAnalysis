# Pythia pi0-main anchor-cluster topology workflow

This workflow starts from every selected central SPLIT CEMC cluster for which a
selected pi0 is the grouped main truth contributor. Each anchor cluster is
filled exactly once into separated, merged, single-contaminated,
missing(energy-threshold), missing(acceptance), missing(other), or other, so
that

    N_pi0-main anchor = N_separated + N_merged + N_single-contaminated
                       + N_missing-energy-threshold + N_missing-acceptance
                       + N_missing-other + N_other

holds both globally and in every cluster-ET bin, including underflow and
overflow. An aggregate missing spectrum is also retained and satisfies

    N_missing = N_missing-energy-threshold + N_missing-acceptance + N_missing-other

in every bin.

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

The anchor requires |eta_cluster| < 0.7. Partner lookup uses the same
configurable cluster-energy cut as the anchor, but its eta cut is disabled by
default (partner_cluster_eta_max <= 0), so every cluster in the CEMC cluster
container can be considered. The default shared energy cut is
E_cluster >= 0.2 GeV.

## Energy-deposit topology

For each direct pi0 daughter photon, the matcher finds the cluster with maximum
absolute daughter energy deposit among clusters passing the shared cluster
selection. A daughter photon is considered recovered only when its
maximum-deposit cluster satisfies the calibrated-energy estimate

    Erec_gamma = Ecluster * (Edep_gamma / Edep_total)

    Erec_gamma / Etruth_gamma >= min_photon_energy_recovery

with default threshold 0.5. The optional cluster-composition requirement

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
- missing: it is the recovered maximum-deposit cluster of one daughter and the
  other daughter has no cluster passing the photon-energy recovery cut. Missing
  is split with the following exclusive priority:

  1. acceptance: the partner projection at the CEMC radius is valid and
     |eta_projection| >= cemc_acceptance_eta_max (default 1.1);
  2. energy-threshold: within missing_diagnostic_max_delta_r (default 0.15) of
     the in-acceptance partner projection, a cluster below min_cluster_energy
     has usable direct daughter deposit;
  3. missing-other: every remaining missing case, including invalid projection,
     best cluster below recovery, incomplete direct matching, or no direct
     deposit;

- other: the anchor is the maximum-deposit cluster of neither daughter, its
  main-contributor assignment is tied, or no preceding definition applies.

The acceptance test is based on the projected partner photon, not on parent or
daughter truth eta. The boundary is exclusive: exactly
|eta_projection| = cemc_acceptance_eta_max is outside. Acceptance takes
priority over the energy-threshold diagnostic.

The anchor-spectrum workflow enables low-threshold missing diagnostics by
default. Set enable_missing_diagnostics=false to avoid the additional matching
cost; then no missing event can be labeled energy-threshold and such events fall
into missing-other unless acceptance applies. Partial metadata records the
acceptance boundary, pre-CEMC interaction radius, direct-match energy-coverage
threshold, diagnostic delta-R, diagnostic enable flag, and topology algorithm
version so incompatible productions cannot be combined.

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
partial schema 7; do not mix it with earlier partial schemas. Point the next
production at a new, empty output directory. run_partial.sh
accepts optional CEMC_ACCEPTANCE_ETA_MAX,
MIN_DIRECT_MATCH_CLUSTER_ENERGY_COVERAGE, MISSING_DIAGNOSTIC_MAX_DELTA_R, and
ENABLE_MISSING_DIAGNOSTICS, and PRE_CEMC_INTERACTION_RADIUS arguments after
MAX_ABS_VERTEX_Z; their defaults are 1.1, 0.5, 0.15, true, and 90.0 cm.

Finalize a complete production with:

    root -l -b -q 'workflows/pi0_anchor_topology/FinalizePythiaPi0AnchorClusterSpectra.C("output/pi0_anchor_topology_partial/eta07_zvtx60_full_partner_fgamma0p0_recovery0p5_clusterenergy/partial_*.root","output/plots/pi0_anchor_topology/minimum_bias/eta07_zvtx60_full_partner_fgamma0p0_recovery0p5_clusterenergy",0,200000,"Pythia8 p+p MB")'

The finalizer writes the combined raw and bin-width-normalized spectra, the
aggregate missing spectrum, and seven exclusive category fractions relative to
the anchor spectrum, output_base.pdf, and
output_base_category_fractions.pdf. It also writes the stacked category-fraction
plot output_base_category_fraction_stack.pdf. Plot annotations and legends are
placed outside the histogram frame.
