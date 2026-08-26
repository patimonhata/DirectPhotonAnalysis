# Topology event display

This package is a diagnostic view of the pi0 anchor-topology decision.  The
decision itself is not implemented in the display: both
`PythiaPi0AnchorClusterSpectrum` and this package call the shared
`photon_tree::Pi0AnchorTopologyEvaluator` in `PhotonAnalysisTree`.

## Main entry points

The four files normally used are:

- `workflows/pythia/RunPythiaTopologyEventDump.C`: read synchronized Pythia
  DST streams and write an intermediate ROOT file.
- `workflows/single_particle/RunSingleParticleTopologyEventDump.C`: do the
  same for a reconstructed single-particle DST.
- `display/BrowseTopologyEvents.C`: overwrite one PDF for the current event;
  Enter advances, `p` goes back, and `q` exits.
- `display/MakeTopologyEventBook.C`: write all selected events and anchor
  detail pages to one multipage PDF.

Code under `display/internal/` is renderer support.  Independent consistency
checks live under `validation/`; they are not analysis entry points.

## Build

From the project root, in an sPHENIX analysis environment:

```sh
PhotonAnalysisTree/src/build.sh
TopologyEventDisplay/src/build.sh
```

The first command installs the shared evaluator.  The second installs the
Fun4All dumper.

## Pythia workflow

The input manifest has the same suffix format used by
`PhotonAnalysisTree/workflows/pi0_anchor_topology`.  Each suffix resolves the
four synchronized catalog files `DST_CALO_CLUSTER_*`, `DST_MBD_EPD_*`,
`DST_TRUTH_JET_*`, and `G4Hits_*`.

Example for manifest row 0 and one event:

```sh
root -l -b -q 'TopologyEventDisplay/workflows/pythia/RunPythiaTopologyEventDump.C("PhotonAnalysisTree/input/minimum_bias/segments.list",0,1,"output/topology_pythia.root",1,true)'
```

Arguments after the output filename are `n_events`, `write_detail`, then the
topology thresholds.  A manifest range makes the extraction naturally
batchable.  For a large sample, first run ranges with `write_detail=false`
(the event/candidate/anchor catalog), select interesting source ranges and
event indices, then rerun only those ranges with `write_detail=true`.
Summary files are intentionally not merged because their local event indices
are only unique together with the manifest range stored in `metadata`.

## Single-particle workflow

A reconstructed single-particle DST normally already contains the truth-cell
and truth-tower nodes:

```sh
root -l -b -q 'TopologyEventDisplay/workflows/single_particle/RunSingleParticleTopologyEventDump.C("input.root","output/topology_gun.root",10)'
```

The optional argument after `write_detail` is `rebuild_truth_nodes`; leave it
`false` for the standard reconstructed DST and use `true` only for a stripped
DST that lacks those nodes.  The intermediate schema and renderer are the same
as for Pythia.

## Viewing

Create one multipage PDF, optionally filtering on topology and pathway:

```sh
root -l -b -q 'TopologyEventDisplay/display/MakeTopologyEventBook.C("output/topology_pythia.root","output/topology_book.pdf",-1,-1,-1)'
```

The next three arguments are `topology_filter`, `pathway_filter`, and
`max_events`.  Use `-1` for no filter.  Topology codes are 0 other,
1 separated, 2 merged, and 3 missing.  Pathway codes are 1 G4-primary pi0
decay, 2 generator-level pi0 decay, and 3 single-particle G4 decay.

Four optional range arguments follow: `vertex_z_min`, `vertex_z_max`,
`truth_pi0_pt_min`, and `truth_pi0_pt_max`.  Bounds are inclusive.  For
example, this selects up to 20 events with collision vertex z in [-10, 10] cm
and at least one selected truth pi0 candidate with pT in [5, 15] GeV/c:

```sh
root -l -b -q 'TopologyEventDisplay/display/MakeTopologyEventBook.C("output/topology_pythia.root","output/selected.pdf",-1,-1,20,-10,10,5,15)'
```

Topology, pathway and truth-pT filters apply to the same pi0 candidate.  Omit
the range arguments to leave these filters disabled.

For sequential inspection in an interactive ROOT session:

```cpp
.L TopologyEventDisplay/display/BrowseTopologyEvents.C
BrowseTopologyEvents("output/topology_pythia.root",
                     "output/current_event.pdf", -1, -1, -1,
                     -10.0, 10.0, 5.0, 15.0);
```

Every event starts with an overview page containing x-y, z-r and eta-phi
views.  One color consistently identifies one selected pi0 family: the pi0 is
solid, while its two daughter-photon branches use fine and coarse dashes.  The
open cluster ring uses the same family color, and the filled cluster marker
encodes the anchor topology.  Each anchor then gets a clearly labeled detail
page showing only that family's G4 truth segments in x-y and z-r, together with
eta-phi, the classification inputs, a local tower energy map and a
self-contained legend.  Detail pages are grouped by pi0
candidate; the stored anchor IDs and intermediate trees are unchanged.  Nested
tower outlines identify the anchor and the two best daughter-photon clusters
using the same line grammar.  Truth segments are straight display guides
between G4 vertices or to the calorimeter radius; they are not
magnetic-field-propagated charged tracks.

## Intermediate ROOT schema

- `metadata`: schema/algorithm versions, sample mode, source range,
  thresholds, detail flag and processing counters.
- `events`: pi0 population, truth/cluster population and topology counts.
- `pi0_candidates`: the two supported Pythia pathways (or gun pathway),
  daughter photons, best clusters and recovery decision.
- `anchor_decisions`: one row per anchor, including topology, reason,
  fractions and energy-deposit inputs.
- `candidate_cluster_truth`: direct daughter-photon deposit and recovery
  quantities for every evaluated candidate/cluster pair.
- `clusters`, `cluster_contributors`: reconstructed clusters and their
  truth contributors.
- `cluster_towers`: member energy, calibrated tower energy and allocation
  fraction, so shared/adjacent cluster boundaries remain visible.
- `truth_particles`, `truth_segments`: all G4 particles plus selected-pi0
  family labels and display segments.

With `write_detail=false`, the four detailed geometry/truth tables may be
empty while event, candidate and anchor tables remain usable for selection.

## Validation

The auxiliary validator is deliberately separated from the main macros:

```sh
root -l -b -q 'TopologyEventDisplay/validation/CheckTopologyEventDump.C("output/topology_pythia.root")'
```

It checks required trees and branches, category closure, anchor references,
key uniqueness and direct energy-deposit closure.
