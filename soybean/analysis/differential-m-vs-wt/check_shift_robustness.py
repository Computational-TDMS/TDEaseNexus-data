# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas"]
# ///
"""Robustness of the shifted-vs-unshifted direction split + tolerant shift families."""

import json
import re

import pandas as pd

t = pd.read_csv("tiers.csv", index_col=0)
feat = pd.read_csv("features.csv", index_col=0)
feat = feat[~feat.index.duplicated()]
SHIFT = re.compile(r"\[([+-]\d+(?:\.\d+)?)\]")

t["has_shift"] = [len(SHIFT.findall(str(p))) > 0 for p in t["proteoform"]]
t["direction"] = ["up in M" if fc > 0 else "up in WT" for fc in t.log2FoldChange]
t["tier"] = "ns"
sig = t.significant == True
t.loc[sig & (t.min_det >= 4), "tier"] = "A"
t.loc[sig & (t.min_det == 3), "tier"] = "B"
t.loc[sig & (t.min_det < 3), "tier"] = "C"

print("== shift vs direction, all significant ==")
print(pd.crosstab(t[t.tier != "ns"].has_shift, t[t.tier != "ns"].direction).to_string())
print("\n== tier A only ==")
print(pd.crosstab(t[t.tier == "A"].has_shift, t[t.tier == "A"].direction).to_string())

# which shift families drive it: tolerant clustering at 0.05 Da
rows = []
for i, r in t.iterrows():
    for v in SHIFT.findall(str(r["proteoform"])):
        rows.append({"featureId": i, "v": float(v), "gene": r["gene"],
                     "tier": r["tier"], "direction": r["direction"],
                     "fc": r["log2FoldChange"], "spectral": r["spectralCount"]})
sh = pd.DataFrame(rows)
sh["family"] = None
centers = []
for v in sorted(sh.v.unique()):
    hit = next((c for c in centers if abs(c - v) <= 0.05), None)
    if hit is None:
        centers.append(v)
        hit = v
    sh.loc[sh.v == v, "family"] = hit
fam = (sh.groupby("family")
         .agg(n_feat=("featureId", "nunique"), n_genes=("gene", "nunique"),
              n_sigA=("tier", lambda s: (s == "A").sum()),
              n_upM=("direction", lambda s: (s == "up in M").sum()),
              n_upWT=("direction", lambda s: (s == "up in WT").sum()),
              max_abs_fc=("fc", lambda s: s.abs().max()),
              spectral=("spectral", "sum"),
              genes=("gene", lambda s: ";".join(sorted(set(map(str, s))))[:70]))
         .sort_values("n_feat", ascending=False))
fam.to_csv("shift_families_tol005.csv")
print("\n== top shift families (0.05 Da tolerance, top 14) ==")
print(fam.head(14).to_string())

print("\n== 317.39 family detail ==")
print(sh[sh.family.between(317.38, 317.40)].sort_values("fc", ascending=False).to_string(index=False))

# excluding the ~16 Da family, does the split hold?
no16 = t[(t.tier != "ns") & ~t.index.isin(sh[sh.family.between(15.9, 16.1)].featureId)]
print("\n== significant, excluding ~16 Da features ==")
print(pd.crosstab(no16.has_shift, no16.direction).to_string())

# 16 Da family split
d16 = t[(t.tier != "ns") & t.index.isin(sh[sh.family.between(15.9, 16.1)].featureId)]
print("\n== ~16 Da features only ==")
print(pd.crosstab(d16.has_shift, d16.direction).to_string())
