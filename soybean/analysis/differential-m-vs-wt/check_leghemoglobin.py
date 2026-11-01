# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas"]
# ///
"""Characterise the LBA / LGB1 proteoform-level direction split."""

import json
import re

import pandas as pd

expr = pd.read_csv("expression.csv", index_col=0)
feat = pd.read_csv("features.csv", index_col=0)
feat = feat[~feat.index.duplicated()]
occ = pd.read_csv("features_with_occupancy.csv", index_col=0)
cmp = pd.DataFrame(json.load(open("comparison_WT_vs_M_full.json"))["candidates"]).set_index("featureId")
t = pd.read_csv("tiers.csv", index_col=0)
t["direction"] = ["up in M" if fc > 0 else "up in WT" for fc in t.log2FoldChange]

SHIFT = re.compile(r"\[([+-]\d+(?:\.\d+)?)\]")
wt = [c for c in expr.columns if "WT0" in c]
mut = [c for c in expr.columns if re.search(r"_M0", c)]

for gene in ["LBA", "LGB1"]:
    ids = [i for i in t.index if str(t.loc[i, "gene"]) == gene]
    d = t.loc[ids].copy()
    d["n_shift"] = [len(SHIFT.findall(str(p))) for p in d["proteoform"]]
    d["has_shift"] = d["n_shift"] > 0
    d["full_length"] = [bool(feat.loc[i, "first_aa"] == 1 and feat.loc[i, "last_aa"] == feat.loc[i, "prot_length"])
                        for i in d.index]
    print(f"== {gene}: {len(d)} tested features, {int((d.significant==True).sum())} significant")
    g = d[d.significant == True]
    print(g.groupby(["has_shift", "direction"]).size().to_string())
    print("  median FC | shifted:", round(g[g.has_shift].log2FoldChange.median(), 2),
          "| unshifted:", round(g[~g.has_shift].log2FoldChange.median(), 2))
    print("  full-length features:", int(g.full_length.sum()), "of", len(g))
    # total intensity per group for the gene, using all its features
    tot_m = expr.loc[ids, mut].sum().sum()
    tot_w = expr.loc[ids, wt].sum().sum()
    print(f"  gene-level summed intensity M {tot_m/1e9:.1f}e9 vs WT {tot_w/1e9:.1f}e9 "
          f"(log2 ratio {__import__('math').log2(tot_m/tot_w):+.2f})")
    # spectra-level
    print("  spectral_count sum:", int(d.spectralCount.sum()))
    print()

# how many genes show both directions among significant features
sig = t[t.significant == True]
both = (sig.groupby("gene").direction.nunique() == 2).sum()
print("genes with significant features in BOTH directions:", int(both),
      "of", sig.gene.nunique())

# shift-status vs direction across all significant features
sig = sig.copy()
sig["has_shift"] = [len(SHIFT.findall(str(p))) > 0 for p in sig["proteoform"]]
print("\nsignificant features with a numeric shift:", int(sig.has_shift.sum()), "of", len(sig))
print(pd.crosstab(sig.has_shift, sig.direction).to_string())
