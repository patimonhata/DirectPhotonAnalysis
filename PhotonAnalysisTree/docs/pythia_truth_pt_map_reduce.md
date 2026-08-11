# Pythia truth-pT map-reduce workflow

## Particle selection

The photon histogram counts prompt photons. A photon must be in the final
HepMC photon collection, have valid classifier output, and have category 1
(direct 2-to-2) or category 2 (fragmentation). Category 3 decay photons are
not included. The pi0 and pi0-decay-photon definitions are documented in
`pythia_truth_spectrum_schema.md`.

The same truth-eta selection is applied independently to each particle species.
The default partial job uses `|eta| < 0.7`, 100 uniform bins, `0 <= pT < 20 GeV`,
and unweighted particle counts.

## Map step

`AccumulatePythiaTruthPtSpectra.C` reads a half-open range `[begin:end]` from
the suffix manifest. Every expected truth-spectrum ROOT file must exist. It
enables only the branches required by the three histograms and writes:

- `h_prompt_photon_truth_pt_raw`
- `h_pi0_truth_pt_raw`
- `h_pi0_decay_photon_truth_pt_raw`
- one-entry `metadata`

The histograms include underflow and overflow and retain `Sumw2`. They are raw:
no bin-width or event-count normalization is applied. Metadata records the exact
manifest range, first and last suffix, input file and event counts, binning,
eta selection, weight policy, and particle counts.

`run_truth_pt_partial_pythia.sh` writes transactionally, runs
`check_pythia_truth_pt_partial.C`, and only then renames the output to
`partial_NNNNNN.root`. Existing partials are never overwritten.

The Condor defaults are 40,000 files, 500 files per partial, and 80 jobs.
`chunk_offset` makes it possible to extend a completed range without
resubmitting existing chunks. The required relation is:

```text
n_chunks = ceil(total_files / files_per_job) - chunk_offset
```

## Reduce and plot step

`FinalizePythiaTruthPtSpectra.C` reads all matching partial files and sorts them
by manifest range. Before adding histograms it requires:

- identical manifest, input directory, particle selection, binning, eta cut,
  and weight policy;
- contiguous, non-overlapping manifest ranges;
- the requested first and final manifest indices;
- valid raw histograms and zero malformed-event diagnostics.

It writes both the summed raw histograms and density histograms. Each density
bin is `raw bin content / bin width`, so its unit is particles/GeV (or weighted
particles/GeV). This normalization is performed only after all partials have
been added. It is not an event-count or luminosity normalization.
