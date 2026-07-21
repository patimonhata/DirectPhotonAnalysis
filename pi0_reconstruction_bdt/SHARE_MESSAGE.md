Hi, I prepared an adapter that adds our PPG15/PPG12 `base_v3E` SPLIT
photon-ID BDT score to the vector-based `event_tree` produced by
`Pi0Reconstruction`.

The package contains the ROOT macros, a one-command wrapper, a 20-event test
input, and a validated reference output. The input ROOT file is not modified;
the adapter writes a separate output file with these new branches:

- `cluster_bdt_base_v3E_split`: BDT score for each cluster
- `cluster_bdt_base_v3E_split_valid`: input-quality flag for each score

Quick start:

```bash
tar -xzf pi0_reconstruction_bdt_adapter.tar.gz
cd pi0_reconstruction_bdt
./run_add_bdt.sh INPUT.root OUTPUT.root
```

Please read `README.md` before physics use. In particular, the supplied smoke
file contains only clusters below 5 GeV, outside the model's documented
performance bins, and the external shower-shape implementation still needs a
line-by-line compatibility check with PhotonAna.
