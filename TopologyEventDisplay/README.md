# Topology event display

This package displays truth cluster origin and diphoton topology for Pythia
events. It uses the same common composition and π⁰/η topology evaluator as
`PhotonAnalysisTree/workflows/photon_candidate_selection`. The single-particle
workflow retains its gun-specific π⁰ anchor selection.

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
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
PhotonAnalysisTree/src/build.sh
TopologyEventDisplay/src/build.sh

photon_prefix="$PWD/PhotonAnalysisTree/install/dev/ana.565"
display_prefix="$PWD/TopologyEventDisplay/install"
export ROOT_INCLUDE_PATH="$display_prefix/include:$photon_prefix/include:${ROOT_INCLUDE_PATH:-}"
export LD_LIBRARY_PATH="$display_prefix/lib64:$photon_prefix/lib64:$photon_prefix/dependencies/Pi0Reconstruction/lib:${LD_LIBRARY_PATH:-}"
```

The two dump macros resolve their libraries using this environment. To build
against a different compatible PhotonAnalysisTree install, set
`PHOTON_ANALYSIS_TREE_PREFIX` when building the display and use the same prefix
in the exports above. The ROOT-only viewing macros can also open older dumps
with the `pi0_candidates` tree.

## Pythia workflow

The input manifest has the same suffix format used by
`PhotonAnalysisTree/workflows/photon_candidate_selection`.  Each suffix resolves the
four synchronized catalog files `DST_CALO_CLUSTER_*`, `DST_MBD_EPD_*`,
`DST_TRUTH_JET_*`, and `G4Hits_*`.

Example for manifest row 0 and one event:

```sh
root -l -b -q 'TopologyEventDisplay/workflows/pythia/RunPythiaTopologyEventDump.C("PhotonAnalysisTree/input/minimum_bias/segments.list",0,1,"output/topology_pythia.root",1,true)'
```

The Pythia macro uses the candidate-map defaults of `E > 0.1 GeV` for stored
clusters and `|z| < 60 cm` for the signal HepMC vertex. By default each meson's
partner threshold follows the cluster threshold. For production comparisons,
pass the settings recorded in the candidate map metadata. For example, the
Jet12 submit settings use `E > 0.12 GeV` for clusters, `E > 0.5 GeV` for both
truth partners, and `0.12 < E <= 0.5 GeV` for the missing energy band.

```sh
root -l -b -q 'TopologyEventDisplay/workflows/pythia/RunPythiaTopologyEventDump.C("PhotonAnalysisTree/input/jet12/segments.list",0,1,"output/topology_jet12.root",40,true,0.7,-1.0,1.1,0.12,0.5,0.0,0.0,0.5,0.15,0,90.0,0.5,0.5,0.10,0.20,0.45,0.65,0.12,0.5)'
```

The final eight arguments set π⁰/η partner thresholds, π⁰ and η mass windows,
and the missing-energy bounds. For a large sample, first run manifest ranges with
`write_detail=false`. The event, candidate, anchor and cluster origin tables
support selection of interesting ranges for a detailed rerun. The local
`event` index is relative to the input range. The map and display share
`event_uid`, derived from the input segment and entry within the file.

## Single-particle workflow

A reconstructed single-particle DST normally already contains the truth-cell
and truth-tower nodes:

```sh
root -l -b -q 'TopologyEventDisplay/workflows/single_particle/RunSingleParticleTopologyEventDump.C("input.root","output/topology_gun.root",10)'
```

After `n_events`, the optional arguments start with `cluster_node`,
`write_detail` and `rebuild_truth_nodes`. Set `rebuild_truth_nodes=true` for a
stripped DST that lacks the truth-cell and truth-tower nodes. The gun uses the
same ROOT schema and renderer. Its `cluster_origins` table is empty and its
anchor `applicability` is -1, the legacy-mode sentinel.

## Viewing

Create one multipage PDF, optionally filtering on topology and pathway:

```sh
root -l -b -q 'TopologyEventDisplay/display/MakeTopologyEventBook.C("output/topology_pythia.root","output/topology_book.pdf",-1,-1,-1)'
```

The next three arguments are `topology_filter`, `pathway_filter`, and
`max_events`. Use `-1` to leave a filter disabled. The topology codes are
0 other, 1 separated, 2 merged, 3 missing, and 4 single-contaminated. In a
Pythia dump these filters require an evaluable anchor. Pathway codes are
1 G4 primary, 2 generator, 3 single-particle gun, 4 legacy G4 η daughter, and
5 G4 secondary.

The four inclusive range arguments are `vertex_z_min`, `vertex_z_max`,
`truth_meson_pt_min`, and `truth_meson_pt_max`. Three optional arguments
follow. `meson_pdg_filter` is 111 for π⁰, 221 for η, or 0 for either.
`origin_filter` and `applicability_filter` use the codes in the [candidate
definitions](../PhotonAnalysisTree/docs/photon_candidate_definitions.md),
with -1 disabling the filter. Origin and applicability alone can find
clusters without anchor rows. When combined with a candidate or topology
filter, they apply to the cluster of the same anchor. Filters select events for
viewing. The PDF shows the full selected event.

This example selects up to 20 events containing a π⁰ candidate with
`5 <= pT <= 15 GeV/c` at `-10 <= z <= 10 cm`.

```sh
root -l -b -q 'TopologyEventDisplay/display/MakeTopologyEventBook.C("output/topology_pythia.root","output/selected.pdf",-1,-1,20,-10,10,5,15,111)'
```

For sequential inspection in an interactive ROOT session:

```cpp
.L TopologyEventDisplay/display/BrowseTopologyEvents.C
BrowseTopologyEvents("output/topology_pythia.root",
                     "output/current_event.pdf", -1, -1, -1,
                     -10.0, 10.0, 5.0, 15.0, 111);
```

Every event starts with an overview page containing x-y, z-r and eta-phi
views.  One color consistently identifies one meson family. The parent line is
solid, while its two daughter-photon branches use fine and coarse dashes.  The
open cluster ring uses the same family color.  Cluster markers below 0.1 GeV use a
separate compressed size scale, while markers at or above 0.1 GeV start from a substantially larger size; the fill color
encodes the anchor topology.  Each anchor then gets a clearly labeled detail
page showing only that family's G4 truth segments in x-y and z-r, together with
eta-phi, the classification inputs, a local tower energy map and a
self-contained legend.  Stars in eta-phi use the daughter-photon momentum
eta and phi, matching the vertex-corrected cluster coordinates.  The separate
CEMC projection eta and phi remain available for acceptance and missing-partner
diagnostics, but are not used as star coordinates.  Detail pages are grouped by
meson candidate. Nested tower outlines identify the anchor and the truth
representative for each daughter, including representatives below the cluster
threshold. Within a meson family, a direct photon with
stored daughters stops at its first daughter-production vertex instead of being
projected to the calorimeter.  Pre-CEMC descendants are drawn as thin branch
lines: successive daughter-production vertices on one parent track are joined
in time order and clipped at the CEMC inner surface.  Fallback projection to the
calorimeter is reserved for terminal meson and direct-photon tracks; terminal
shower descendants without a child vertex are omitted.  Truth segments are
straight display guides between G4 vertices or to the calorimeter radius; they
are not magnetic-field-propagated charged tracks.

## Cluster origin and topology

Pythia stores a common-composition row for each cluster above the map energy
threshold, including clusters outside anchor η acceptance. The selected origin
uses the common majority rule, with a strict fraction above 0.5 and the original
cluster truth deposit as denominator. η parent contributions have their own
origin classes. Direct π⁰ and η diphotons use the matching topology evaluator
for that meson. See the [candidate definitions](../PhotonAnalysisTree/docs/photon_candidate_definitions.md)
for the full origin and topology codes.

An anchor associates a cluster with a concrete diphoton candidate. The
cluster's `applicability` controls whether to count that topology in the event
summary. The values are 0 unknown, 1 evaluable, and 2 not applicable. For example,
a confirmed non-diphoton π⁰ has a cluster origin row with applicability 2.
A diphoton can have an anchor and still be unknown because its anchor match,
projection, or interaction information is unavailable.

Matching can use a complete match or a partial match with member-energy coverage
of at least 0.5 by default, with complete hit ancestry. A failed match in
another cluster leaves the current anchor's applicability intact. The stored
candidate-wide minimum coverage is diagnostic. A ray with a known lack of
intersection with the CEMC cylinder can produce a missing partner classified
as acceptance. An unavailable projection produces unknown applicability.

Truth partner selection uses the available positive-energy cluster with the
largest daughter deposit, with ties broken by the smallest cluster ID. The
partner energy and optional photon recovery cuts are applied to that
representative. The default recovery requirement is disabled. The PDF and ROOT
output retain the representative and the nearby matching diagnostic.
The single-contaminated category requires exactly one photon with its first
daughter-production radius strictly below the default 90 cm boundary.

## Intermediate ROOT schema

Schema 7 includes these classification tables.

- `metadata` records the algorithm versions, input range, thresholds, detail
  flag and processing counters.
- `events` records local and source event IDs, π⁰ and η candidate populations,
  topology counts and unevaluable anchor counts.
- `cluster_origins` records the common class, selected particle identity and PDG,
  fraction, original and unresolved truth deposit, quality, applicability and
  topology reason. It also retains candidate-wide and anchor-local match metrics.
- `meson_candidates` stores the meson PDG, canonical parent and photon identities,
  pathway, truth kinematics, projection and interaction diagnostics, and the two
  truth representatives. A representative can be below the cluster threshold.
- `anchor_decisions` connects each candidate and anchor cluster. It stores the
  topology code, applicability and reasons, truth partner, and direct-match
  diagnostics. Interpret a Pythia topology code together with applicability.

Detailed dumps also populate six tables. `candidate_cluster_truth` stores direct
daughter matches for evaluated candidates and above-threshold clusters.
`clusters`, `cluster_contributors`, and `cluster_towers` provide reconstructed
geometry, contributors and tower allocation. `truth_particles` and
`truth_segments` provide the G4 geometry and family labels. With Pythia, the
evaluator runs on anchored candidates in both detail modes. The gun can
evaluate all candidates when detail is enabled.

## Validation

Run the ROOT validator on an intermediate dump. To compare a Pythia dump with a
candidate map, use the same manifest and thresholds and run the second command.

```sh
root -l -b -q 'TopologyEventDisplay/validation/CheckTopologyEventDump.C("output/topology_jet12.root")'
python3 TopologyEventDisplay/validation/compare_dump_to_candidate_map.py output/topology_jet12.root output/candidate_map.root
```

The validator checks counts, references, evaluation status, missing-category
priority and direct energy-deposit closure. The comparison uses
`(event_uid, cluster_id)` to verify origin, the selected identity, matching and
topology diagnostics, and truth partners against the map.
