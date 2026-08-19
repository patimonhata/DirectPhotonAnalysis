# Pythia pi0-main anchor-cluster topology workflow

This workflow starts from every selected central SPLIT CEMC cluster for which a
selected pi0 is the grouped main truth contributor. Each anchor cluster is
filled exactly once into separated, merged, missing, or other, so that

    N_pi0-main anchor = N_separated + N_merged + N_missing + N_other

holds both globally and in every cluster-ET bin, including underflow and
overflow.

## Anchor and pi0 definition

For each cluster, contributor fractions compatible with the same selected pi0
are summed. The pi0 with the largest grouped fraction is the main contributor
when that fraction is at least anchor_pi0_fraction_min (default 0.5).
Equal-leading assignments are retained in the anchor denominator but classified
as other.

Selected pi0s are either transported G4-primary pi0 decays or generator pi0
decays represented by exactly two G4-primary photons. Detector-secondary pi0s
and Dalitz decays are excluded. The default truth-parent acceptance is
|eta_pi0| < 0.7.

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

- merged: it is the recovered maximum-deposit cluster of both daughter photons;
- separated: it is the recovered maximum-deposit cluster of one daughter and
  the other daughter has a distinct recovered maximum-deposit partner cluster;
- missing: it is the recovered maximum-deposit cluster of one daughter and the
  other daughter has no cluster passing the photon-energy recovery cut;
- other: the anchor is the maximum-deposit cluster of neither daughter, its
  main-contributor assignment is tied, or no preceding definition applies.

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
transactionally and validates its partial before publication.

Finalize a complete production with:

    root -l -b -q 'workflows/pi0_anchor_topology/FinalizePythiaPi0AnchorClusterSpectra.C("output/pi0_anchor_topology_partial/eta07_full_partner_fgamma0p0_recovery0p5_clusterenergy/partial_*.root","output/plots/pi0_anchor_topology/minimum_bias/eta07_full_partner_fgamma0p0_recovery0p5_clusterenergy",0,200000,"Pythia8 p+p MB")'

The finalizer writes the combined raw and bin-width-normalized spectra, category
fractions relative to the anchor spectrum, output_base.pdf, and
output_base_category_fractions.pdf. It also writes the stacked category-fraction
plot output_base_category_fraction_stack.pdf. Plot annotations and legends are
placed outside the histogram frame.
