# Pythia truth-pT map-reduce workflow

This directory is the complete truth-pT analysis workflow:

- `PythiaTruthPtSpectrum`: accumulate truth-pT histograms directly from G4Hits DST nodes;
- `Fun4All_PythiaTruthPtSpectra.C`: map a manifest range of G4Hits DSTs to one partial;
- `check_pythia_truth_pt_partial.C`: validate each partial;
- `FinalizePythiaTruthPtSpectra.C`: reduce partials and write ROOT/PDF outputs;
- `run_partial.sh` and `submit.job`: execute and submit the map jobs.

## Particle selection

The photon histogram counts prompt photons. A photon must be in the final
HepMC photon collection, have valid classifier output, and have category 1
(direct 2-to-2) or category 2 (fragmentation). Category 3 decay photons are
not included. The pi0 and pi0-decay-photon definitions are documented in
`../../docs/pythia_truth_spectrum_schema.md`.

The pi0-decay-photon histogram is the union of two disjoint simulation-stage
components:

- final HepMC photons whose valid photon-copy chain reaches a production vertex
  with exactly one PDG 111 parent;
- G4 photons whose immediate parent is a signal G4-primary pi0 matched to the
  HepMC event.

The HepMC rule does not require a two-photon daughter topology, so it includes
the photon from both `pi0 -> gamma gamma` and Dalitz
`pi0 -> gamma e+ e-` decays. Component histograms and counts are retained to
check the generator/Geant4 decay handoff.

The same truth-eta selection is applied independently to each particle species.
The default partial job uses `|eta| < 0.7`, 100 uniform bins, `0 <= pT < 20 GeV`,
and unweighted particle counts.

## Map step

`Fun4All_PythiaTruthPtSpectra.C` reads a half-open range `[begin:end]` from
the suffix manifest, adds each `G4Hits_<suffix>` directly to a Fun4All input
manager, and writes:

- `h_prompt_photon_truth_pt_raw`
- `h_pi0_truth_pt_raw`
- `h_pi0_decay_photon_truth_pt_raw`, the HepMC+G4 total
- `h_hepmc_pi0_decay_photon_truth_pt_raw`
- `h_g4_pi0_decay_photon_truth_pt_raw`
- one-entry schema-v3 `metadata`

The histograms include underflow and overflow and retain `Sumw2`. They are raw:
no bin-width or event-count normalization is applied. Metadata records the exact
manifest range, first and last suffix, G4Hits prefix and node configuration,
input file and event counts, binning, eta selection, weight policy, and total
plus stage-specific particle counts. The checker requires the total histogram
and count to equal the HepMC and G4 component sums.

`run_partial.sh` writes transactionally, runs
`check_pythia_truth_pt_partial.C`, and only then renames the output to
`partial_NNNNNN.root`. Existing partials are never overwritten.

The Condor defaults are 200,000 G4Hits DSTs, 50 files per partial, and 4,000 jobs.
The schema-v3 default output directory is
`truth_pt_partial/minimum_bias/prompt_eta07_unweighted_inclusive_pi0_decay`,
keeping direct-DST partials separate from old tree-based partials.
`chunk_offset` makes it possible to extend a completed range without
resubmitting existing chunks. The required relation is:

```text
n_chunks = ceil(total_files / files_per_job) - chunk_offset
```

## Run

From the `PhotonAnalysisTree` directory, submit the full production with:

```bash
condor_submit -maxjobs 4000 workflows/truth_pt/submit.job
```

Reduce the complete contiguous range with:

```bash
root -l -b -q 'workflows/truth_pt/FinalizePythiaTruthPtSpectra.C("output/truth_pt_partial/minimum_bias/prompt_eta07_unweighted_inclusive_pi0_decay/partial_*.root","output/plots/minbias_truth_pt_prompt_eta07_inclusive_pi0_decay",0,200000)'
```

## Reduce and plot step

`FinalizePythiaTruthPtSpectra.C` reads all matching partial files and sorts them
by manifest range. Before adding histograms it requires:

- identical manifest, input directory, particle selection, binning, eta cut,
  and weight policy;
- contiguous, non-overlapping manifest ranges;
- the requested first and final manifest indices;
- valid raw histograms and zero malformed-event diagnostics.

It writes summed raw and density histograms for the three plotted species and
for the separate HepMC/G4 pi0-decay-photon components. Each density bin is
`raw bin content / bin width`, so its unit is particles/GeV (or weighted
particles/GeV). This normalization is performed only after all partials have
been added. It is not an event-count or luminosity normalization. The PDF plots
the inclusive pi0-decay-photon sum with a `counts/bin` y-axis.
