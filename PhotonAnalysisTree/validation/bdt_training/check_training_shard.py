"""Compare an extraction shard with its original schema-11 maps, including excluded rows."""

import argparse
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


def check_sample(tree, meta, offset):
    features = {
        "cluster_Et": "et", "cluster_Eta": "eta", "cluster_weta_cogx": "shower_w_eta_cogx",
        "cluster_wphi_cogx": "shower_w_phi_cogx", "e11_over_e33": "shower_e11_over_e33",
        "e32_over_e35": "shower_e32_over_e35", **{f"cluster_et{i}": f"shower_et{i}" for i in range(1, 5)},
    }
    sample = meta["sample_name"]
    photonjet = sample.startswith("photonjet")
    common = ["id", "common_main", "common_selected_pdg", "common_applicability", "pi0_anchor_valid", "pi0_anchor_topology",
              "eta_anchor_valid", "eta_anchor_topology", "shower_valid"]
    if photonjet:
        common += ["truth_iso_valid", "truth_pass_isolated"]
    event = ["event_uid", "mbd_vertex_z", "event_weight_valid", "sample_stitching_valid", "sample_stitching_pass", "weight_numerator_pb"]
    branches = event + ["split_cluster_" + key for key in [*common, *features.values()]]
    for index in range(meta["map_begin"], meta["map_end"]):
        path = Path(meta["source_map_root"]) / sample / f"map_{index:06}.root"
        with uproot.open(path) as source:
            for chunk in source["event_tree"].iterate(branches, step_size="32 MB"):
                et = chunk["split_cluster_et"]
                keep = (abs(chunk.mbd_vertex_z) < 60) & (abs(chunk.split_cluster_eta) < 0.7) & (et > 3) & (et < 35)

                def selected(key):
                    values, _ = ak.broadcast_arrays(chunk[key], et)
                    return ak.to_numpy(ak.flatten(values[keep]))

                ids = selected("split_cluster_id")
                rows = tree.arrays(entry_start=offset, entry_stop=offset + len(ids), library="np")
                offset += len(ids)
                if len(ids) == 0:
                    continue

                def equal(key, expected):
                    np.testing.assert_array_equal(rows[key], expected, err_msg=f"{path}: {key}")

                equal("sample_name", np.full(len(ids), sample))
                equal("cluster_id", ids)
                equal("event_uid", selected("event_uid"))
                equal("vertexz", selected("mbd_vertex_z"))
                valid = selected("split_cluster_shower_valid").astype(bool)
                for name, original in features.items():
                    values = selected("split_cluster_" + original)
                    equal(name, values)
                    valid &= np.isfinite(values)
                equal("inputs_valid", valid)
                eligible = (selected("event_weight_valid") & selected("sample_stitching_valid") & selected("sample_stitching_pass")).astype(bool)
                equal("event_eligible", eligible)
                equal("weight_pb", np.where(eligible, selected("weight_numerator_pb") / meta["sample_sum_generator_weights"], 0))
                main = selected("split_cluster_common_main")
                pdg = selected("split_cluster_common_selected_pdg")
                overview = np.select([main == 1, main == 2, (main == 3) | (np.isin(main, [7, 8]) & (pdg == 221))], [0, 1, 2], 3)
                equal("common_main", main)
                equal("common_origin", overview)
                applicability = selected("split_cluster_common_applicability")
                equal("topology_applicability", applicability)
                topology = np.full(len(ids), -1)
                for code, parent in [(2, "pi0"), (3, "eta")]:
                    prefix = f"split_cluster_{parent}_anchor_"
                    use = (main == code) & (applicability == 1) & (selected(prefix + "valid") == 1)
                    topology[use] = selected(prefix + "topology")[use]
                equal("topology", topology)
                iso_valid = selected("split_cluster_truth_iso_valid").astype(bool) if photonjet else np.zeros(len(ids), bool)
                isolated = iso_valid & selected("split_cluster_truth_pass_isolated").astype(bool) if photonjet else iso_valid
                equal("truth_iso_valid", iso_valid)
                equal("truth_pass_isolated", isolated)

                base = np.full(len(ids), -1)
                if photonjet:
                    base[eligible & (main == 1) & isolated] = 1
                else:
                    base[eligible & (((overview == 3) & (main != 4)) | np.isin(topology, [2, 4]))] = 0
                    base[eligible & (topology == 3)] = 1
                equal("label_v1", base)
                if not photonjet:
                    base[eligible & (main == 4)] = 1
                equal("label_v2", base)
                if not photonjet:
                    base[eligible & (topology == 1)] = 1
                equal("label_v3", base)

                sample_hash = 14695981039346656037
                for byte in sample.encode("ascii"):
                    sample_hash = ((sample_hash ^ byte) * 1099511628211) & ((1 << 64) - 1)
                h = (np.uint64(sample_hash) ^ rows["event_uid"]) + np.uint64(0x9e3779b97f4a7c15)
                h = (h ^ (h >> 30)) * np.uint64(0xbf58476d1ce4e5b9)
                h = (h ^ (h >> 27)) * np.uint64(0x94d049bb133111eb)
                bucket = (h ^ (h >> 31)) % 100
                equal("split", np.select([bucket < 70, bucket < 85], [0, 1], 2))
    return offset


def check(path):
    with uproot.open(path) as output:
        tree = output["clusters"]
        offset = 0
        for meta in ak.to_list(output["metadata"].arrays()):
            offset = check_sample(tree, meta, offset)
        assert offset == tree.num_entries, f"expected {offset} rows, found {tree.num_entries}"
        print(f"Checked all {offset} rows in {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shards", type=Path, nargs="+")
    for shard in parser.parse_args().shards:
        check(shard)
