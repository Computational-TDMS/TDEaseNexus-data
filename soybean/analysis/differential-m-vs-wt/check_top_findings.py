# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas"]
# ///
"""Targeted checks on the top differential findings (WT vs M)."""

import json
import re

import pandas as pd

pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 70)
expr = pd.read_csv("expression.csv", index_col=0)
feat = pd.read_csv("features.csv", index_col=0)
occ = pd.read_csv("features_with_occupancy.csv", index_col=0)
cmp = pd.DataFrame(json.load(open("comparison_WT_vs_M_full.json"))["candidates"]).set_index("featureId")
t = pd.read_csv("tiers.csv", index_col=0)
wt = [c for c in expr.columns if "WT0" in c]
mut = [c for c in expr.columns if re.search(r"_M0", c)]

print("== top tier-A by |FC|: raw intensities (Mcps), M then WT ==")
sig = t[t.significant == True]
top = sig.reindex(sig.log2FoldChange.abs().sort_values(ascending=False).index).head(8)
for fid, r in top.iterrows():
    row = expr.loc[fid]
    print(f"{fid:16s} {str(r.gene):16s} FC={r.log2FoldChange:+5.2f} nWT={r.n_WT} nM={r.n_M} | M: "
          + " ".join(f"{row[c]/1e6:5.1f}" for c in mut) + " | WT: "
          + " ".join(f"{row[c]/1e6:5.1f}" for c in wt))

print("\n== EKN / 100527608 features (the +317.39 family) ==")
for g in ["EKN", "100527608"]:
    ids = [i for i in feat.index if str(feat.loc[i, "gene"]) == g]
    print(f"-- {g}: {len(ids)} features")
    for i in ids:
        f = feat.loc[i]
        fc = cmp.loc[i, "log2FoldChange"] if i in cmp.index else float("nan")
        fdr = cmp.loc[i, "bhFdr"] if i in cmp.index else float("nan")
        print(f"  {i:16s} mass={f['mass']:8.1f} fc={fc:+6.2f} fdr={fdr:.2g} "
              f"nWT={occ.loc[i,'n_WT']} nM={occ.loc[i,'n_M']}  {str(f['proteoform'])[:110]}")
