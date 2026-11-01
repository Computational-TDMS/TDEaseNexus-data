# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas"]
# ///
"""Derive per-group detection counts (n_WT, n_M) from the TopPICR expression matrix.

The toppicr-py skill's shortlist tool expects an occupancy table indexed by
featureId with columns n_<group>. The pipeline does not emit one, so this
script derives it directly from expression.csv (value > 0 == detected) and
joins the feature metadata needed for tier inspection.

Provenance: expression.csv and features.csv come from
toppic_run.py --output-dir analysis/differential-m-vs-wt (same directory).
"""

import re

import pandas as pd

expr = pd.read_csv("expression.csv", index_col=0)
feat = pd.read_csv("features.csv", index_col=0)

groups = {
    "WT": [c for c in expr.columns if re.search(r"WT0[1-5]", c)],
    "M": [c for c in expr.columns if re.search(r"M0[1-5]", c)],
}
assert len(groups["WT"]) == 5 and len(groups["M"]) == 5, groups

occ = pd.DataFrame(index=expr.index)
for name, cols in groups.items():
    occ[f"n_{name}"] = (expr[cols] > 0).sum(axis=1)
occ["n_total"] = (expr > 0).sum(axis=1)

meta_cols = [
    "gene", "proteoform", "uniprot_acc", "mass", "rt", "spectral_count",
    "n_unexpected_mods", "proteoform_level", "collision",
    "first_aa", "last_aa", "prot_length", "ann_type",
]
occ = occ.join(feat[[c for c in meta_cols if c in feat.columns]])
occ.index.name = "featureId"
occ.to_csv("features_with_occupancy.csv")
print(f"wrote features_with_occupancy.csv: {occ.shape[0]} features x {occ.shape[1]} cols")
print(occ[["n_WT", "n_M"]].value_counts().sort_index().to_string())
