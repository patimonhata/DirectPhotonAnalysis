# Tagging-purity diagnostic

This workflow reads the existing schema-5 photon-candidate map trees. It does not alter the production selection.

For every Region-A anchor and 1 GeV ET bin it records raw and cross-section-weighted counts before tagging, pi0-only vetoes, eta-only vetoes, overlaps, and final survivors. The same flow is split into prompt and background anchors. Truth-pi0 anchors are further split by truth-pair taggability, combinatorial-only veto, selected-partner truth matching, displaced selected truth partners, and missing category.

The selected pi0 partner ID is joined to the stored split-cluster arrays. The output compares prompt combinatorial partners with selected true, ambiguous, and combinatorial-only truth-pi0 partners in partner E/ET/eta, deltaR, energy asymmetry, mass, event cluster multiplicity, ntower, shower validity, containment, tower-data completeness, shower shapes, and the diagnostic BDT score. The BDT value is descriptive only because low-energy partners can be outside its training range.

The existing 0.2 GeV maps contain the selected partners and cluster features needed for this study. The generic all-pairs and tower-detail arrays are not populated. Alternative pair cuts can be recomputed from cluster four-vectors for thresholds at or above 0.2 GeV, but detector bad/dead-tower proximity cannot be reconstructed from these maps.

## Smoke test

~~~bash
root -l -b -q \
  -e '.L PhotonAnalysisTree/workflows/tagging_purity_diagnostic/ReduceTaggingPurityDiagnostic.C' \
  -e 'gSystem->Exit(ReduceTaggingPurityDiagnostic("PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p2","PhotonAnalysisTree/output/qa/tagging_purity_diagnostic/jet8_10000events.root",0.2,"jet8",0,1,true,10000))'
~~~

## Full Jet reduction

The two submit files create 12 ordinary-sample jobs and 20 Jet12 shards, covering both production thresholds.

~~~bash
mkdir -p PhotonAnalysisTree/output/condor/tagging_purity_diagnostic/reduce
condor_submit PhotonAnalysisTree/workflows/tagging_purity_diagnostic/submit_reduce_jet_samples.job
condor_submit PhotonAnalysisTree/workflows/tagging_purity_diagnostic/submit_reduce_jet12_shards.job
~~~

No script submits jobs automatically. Existing partial outputs are never overwritten.

## Merge and tabulation

After all partials finish:

~~~bash
PhotonAnalysisTree/workflows/tagging_purity_diagnostic/run_merge.sh cluster_e_gt_0p2 \
  PhotonAnalysisTree/output/plots/tagging_purity_diagnostic/cluster_e_gt_0p2/jet
PhotonAnalysisTree/workflows/tagging_purity_diagnostic/run_merge.sh cluster_e_gt_0p5 \
  PhotonAnalysisTree/output/plots/tagging_purity_diagnostic/cluster_e_gt_0p5/jet
~~~

Each merge writes the summed ROOT file plus long-form flow and metric TSV files. Fraction errors use the weighted-subset covariance formula. The prompt definition and event weights match the production composition reducer: strict prompt contribution above 0.5 and cross section times generator weight divided by the full per-sample generator-weight sum.

## Paired pi0/eta threshold study

The paired reducer joins the 0.2 and 0.5 GeV maps by event UID and cluster ID and compares four working points on a common Region-A denominator:
0.5/0.5, 0.2/0.2, pi0=0.2 with eta=0.5, and the inverse. It also stores common survivors so comparison errors include overlap covariance.

~~~bash
mkdir -p PhotonAnalysisTree/output/condor/tagging_purity_diagnostic/hybrid_reduce
condor_submit PhotonAnalysisTree/workflows/tagging_purity_diagnostic/submit_hybrid_jet_samples.job
condor_submit PhotonAnalysisTree/workflows/tagging_purity_diagnostic/submit_hybrid_jet12_shards.job
PhotonAnalysisTree/workflows/tagging_purity_diagnostic/run_hybrid_merge.sh \
  PhotonAnalysisTree/output/plots/tagging_purity_diagnostic/hybrid_threshold/jet
~~~

The merge writes ROOT, flow TSV, metric TSV, and correlated-comparison TSV outputs. See `RESULTS.md` for the current full-statistics result.

## Exploratory selected-partner scan

~~~bash
root -l -b -q \
  -e '.L PhotonAnalysisTree/workflows/tagging_purity_diagnostic/ScanSelectedPartnerFeatures.C' \
  -e 'gSystem->Exit(ScanSelectedPartnerFeatures("INPUT.root","WORKING_POINTS.tsv",0.9,"ROC.tsv"))'
~~~

This scan is integrated over anchor ET and evaluates one binned variable at a time among already selected partners; it is not a replacement for

Current inputs use photon-candidate map schema 5 / topology algorithm 11 (representative partner selected before the tagging energy cut). The six missing categories follow photon_candidate_selection; energy-band taggability uses only the production pi0 mass window. Map settings are retained in diagnostic metadata, checked within each reduce, and compared across merged partials. The hybrid study still explicitly compares common pi0/eta thresholds 0.2 and 0.5 GeV and requires matching mass windows and diagnostic settings. Regenerate older diagnostic products before using these reducers.
