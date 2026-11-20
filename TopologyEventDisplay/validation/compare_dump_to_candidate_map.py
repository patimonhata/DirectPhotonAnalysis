#!/usr/bin/env python3
"""Compare a schema-7 Pythia display dump to the corresponding schema-10 candidate map."""

import argparse
import math
from collections import Counter
from pathlib import Path

import ROOT


ORIGIN_FIELDS = {
    "common_main": "main", "selected_pdg": "selected_pdg", "selected_embedding": "selected_embedding",
    "selected_source": "selected_source", "selected_id": "selected_id", "common_fraction": "fraction",
    "original_truth_total": "total", "unresolved_edep": "unresolved", "quality": "quality",
    "applicability": "applicability", "topology_reason": "topology_reason",
    "candidate_matching_status": "matching_status", "candidate_minimum_coverage": "minimum_coverage",
    "anchor_matching_status": "anchor_matching_status", "anchor_matching_coverage": "anchor_matching_coverage",
    "anchor_matching_usable": "anchor_matching_usable",
}
ANCHOR_FIELDS = {
    "topology": "topology", "reason": "reason", "missing_category": "missing_category",
    "missing_detail": "missing_detail", "partner_photon_index": "partner_photon_index", "main_fraction": "main_fraction",
    "truth_partner_cluster_id": "truth_partner_cluster_id", "truth_partner_cluster_energy": "truth_partner_cluster_e",
    "truth_partner_direct_edep": "truth_partner_direct_edep", "truth_partner_invariant_mass": "truth_partner_mass",
}


def equal(left, right):
    if isinstance(right, str) and len(right) == 1 and isinstance(left, (bool, int)):
        right = ord(right)  # PyROOT may expose vector<unsigned char> elements as one-character strings.
    if isinstance(left, float) or isinstance(right, float):
        return math.isclose(left, right, rel_tol=1e-6, abs_tol=1e-9)
    return left == right


def rows(tree, key_fields, value_fields):
    result = {}
    for row in tree:
        key = tuple(int(getattr(row, field)) for field in key_fields)
        if key in result:
            raise ValueError(f"duplicate {tree.GetName()} key {key}")
        result[key] = {field: getattr(row, field) for field in value_fields}
    return result


def check_fields(label, actual, expected):
    for field, value in expected.items():
        if not equal(actual[field], value):
            raise ValueError(f"{label}, {field}: dump={actual[field]!r}, map={value!r}")


def compare(dump, reference):
    dm, mm = dump.Get("metadata"), reference.Get("metadata")
    dm.GetEntry(0)
    mm.GetEntry(0)
    if dm.schema_version != 7 or mm.schema_version != 10 or not dm.has_input_provenance or dm.sample_mode != 1:
        raise ValueError("expected a schema-7 Pythia dump with input provenance and a schema-10 candidate map")
    if Path(str(dm.manifest_path)).resolve() != Path(str(mm.input_manifest)).resolve():
        raise ValueError("dump and map must use the same input manifest")
    fields = ("min_cluster_energy", "max_abs_vertex_z", "pi0_partner_min_energy", "eta_partner_min_energy",
              "pi0_mass_min", "pi0_mass_max", "eta_mass_min", "eta_mass_max", "missing_energy_min", "missing_energy_max",
              "min_photon_energy_recovery", "common_composition_algorithm_version", "eta_topology_algorithm_version")
    settings = {name: name for name in fields}
    settings.update(anchor_cluster_eta_max="candidate_abs_eta_max", evaluator_algorithm_version="pi0_topology_algorithm_version")
    check_fields("metadata", {name: getattr(dm, name) for name in settings}, {name: getattr(mm, branch) for name, branch in settings.items()})

    events = rows(dump.Get("events"), ("event",), ("event_uid",))
    event_uids = {row["event_uid"] for row in events.values()}
    if not event_uids or len(event_uids) != len(events):
        raise ValueError("expected unique dump event UIDs and at least one event")
    origins = rows(dump.Get("cluster_origins"), ("event_uid", "cluster_id"), ORIGIN_FIELDS)
    candidates = rows(dump.Get("meson_candidates"), ("event", "candidate_id"),
                      ("parent_pdg_id", "parent_embedding", "parent_source", "parent_id", "pathway"))
    anchors = {
        (events[(event,)]["event_uid"], cluster): row
        for (event, cluster), row in rows(dump.Get("anchor_decisions"), ("event", "cluster_id"), ("candidate_id", *ANCHOR_FIELDS)).items()
    }
    event_id_by_uid = {row["event_uid"]: key[0] for key, row in events.items()}
    counts = Counter()
    for row in reference.Get("event_tree"):
        uid = row.event_uid
        if uid not in event_id_by_uid:
            continue
        event_uids.remove(uid)
        counts["events"] += 1
        for index, cluster_id in enumerate(row.split_cluster_id):
            key = (uid, cluster_id)
            if key not in origins:
                raise ValueError(f"missing dump origin {key}")
            expected = {name: getattr(row, "split_cluster_common_" + branch)[index] for name, branch in ORIGIN_FIELDS.items()}
            check_fields(f"origin {key}", origins.pop(key), expected)
            counts["clusters"] += 1
            for meson, pdg in (("pi0", 111), ("eta", 221)):
                prefix = f"split_cluster_{meson}_anchor_"
                if getattr(row, prefix + "valid")[index] in (0, "\x00"):
                    continue
                if key not in anchors:
                    raise ValueError(f"missing {meson} dump anchor {key}")
                anchor = anchors.pop(key)
                expected = {name: getattr(row, prefix + branch)[index] for name, branch in ANCHOR_FIELDS.items()}
                check_fields(f"{meson} anchor {key}", anchor, expected)
                candidate_index = getattr(row, prefix + "candidate_index")[index]
                candidate = candidates[(event_id_by_uid[uid], anchor["candidate_id"])]
                expected = {
                    name: getattr(row, f"{meson}_topology_candidate_{name}")[candidate_index]
                    for name in ("parent_embedding", "parent_source", "parent_id", "pathway")
                }
                check_fields(f"{meson} candidate {key}", candidate, dict(expected, parent_pdg_id=pdg))
                counts[f"{meson} anchors"] += 1
    if event_uids or origins or anchors:
        raise ValueError(f"unmatched dump records: {len(event_uids)} events, {len(origins)} clusters, {len(anchors)} anchors")
    print(", ".join(f"{count} {label}" for label, count in counts.items()) + " agree")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dump")
    parser.add_argument("candidate_map")
    args = parser.parse_args()
    try:
        with ROOT.TFile.Open(args.dump, "READ") as dump, ROOT.TFile.Open(args.candidate_map, "READ") as reference:
            compare(dump, reference)
    except (AttributeError, KeyError, OSError, ValueError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
