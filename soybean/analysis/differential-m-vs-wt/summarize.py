# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas"]
# ///
"""Quantitative summary of the WT vs M comparison (soybean nodule top-down).

Inputs (all in this directory, produced by toppicr-py):
  expression.csv, features.csv           -- toppic_run.py build
  comparison_WT_vs_M_full.json           -- toppic_plot.py compare --top 5000
  features_with_occupancy.csv            -- build_occupancy.py
  tiers.csv, gene_concordance.csv        -- toppic_shortlist.py

Writes: summary.json, sample_detection.csv, presence_only.csv,
shift_families.csv, shift_feature_rows.csv, top_effect_all_tiers.csv
"""

import json
import re

import pandas as pd

CMP = "comparison_WT_vs_M_full.json"
cmp = json.load(open(CMP))
cand = pd.DataFrame(cmp["candidates"]).set_index("featureId")
expr = pd.read_csv("expression.csv", index_col=0)
feat = pd.read_csv("features.csv", index_col=0)
occ = pd.read_csv("features_with_occupancy.csv", index_col=0)
tiers = pd.read_csv("tiers.csv", index_col=0)
# derive the tier label the shortlist printed (shortlist keeps min_det + significant)
tiers["tier"] = "not significant"
sig_m = tiers["significant"] == True
tiers.loc[sig_m & (tiers.min_det >= 4), "tier"] = "A"
tiers.loc[sig_m & (tiers.min_det == 3), "tier"] = "B"
tiers.loc[sig_m & (tiers.min_det < 3), "tier"] = "C"
tiers["direction"] = "up in M"
tiers.loc[tiers.log2FoldChange < 0, "direction"] = "up in WT"

wt = [c for c in expr.columns if re.search(r"WT0[1-5]", c)]
mut = [c for c in expr.columns if re.search(r"M0[1-5]", c)]

out = {}
out["design"] = {
    "groupA_WT": wt, "groupB_M": mut,
    "n_features_matrix": int(expr.shape[0]),
    "tested": cmp["testedCount"], "eligible": cmp["eligibleFeatureCount"],
    "bh_significant": cmp["bhSignificantCount"],
    "up_in_M": cmp["enrichedInGroupB"], "up_in_WT": cmp["enrichedInGroupA"],
    "method": cmp["method"], "definitions": cmp["definitions"],
}

# --- per-sample detection -------------------------------------------------
det = pd.DataFrame({
    "sample": expr.columns,
    "group": ["M" if re.search(r"M0[1-5]", c) else "WT" for c in expr.columns],
    "features_detected": [(expr[c] > 0).sum() for c in expr.columns],
    "total_intensity": [expr[c].sum() for c in expr.columns],
    "median_positive_intensity": [expr.loc[expr[c] > 0, c].median() for c in expr.columns],
})
det.to_csv("sample_detection.csv", index=False)
out["sample_detection"] = det.to_dict("records")

# --- presence-only (untested) features ------------------------------------
pres = occ[(occ["n_WT"] == 0) | (occ["n_M"] == 0)].copy()
pres["side"] = ["only_M" if r.n_M > 0 else "only_WT" for r in pres.itertuples()]
pres["all_in_group"] = [
    "both groups have >=1" if (r.n_WT > 0 and r.n_M > 0) else
    ("M only" if r.n_M > 0 else "WT only") for r in pres.itertuples()
]
pres_out = pres[pres["all_in_group"] != "both groups have >=1"].copy()
pres_out.sort_values("spectral_count", ascending=False).to_csv("presence_only.csv")
summary_pres = {
    "only_M": int((pres_out["all_in_group"] == "M only").sum()),
    "only_WT": int((pres_out["all_in_group"] == "WT only").sum()),
    "n_genes_only_M": int(pres_out.loc[pres_out["all_in_group"] == "M only", "gene"].nunique()),
    "n_genes_only_WT": int(pres_out.loc[pres_out["all_in_group"] == "WT only", "gene"].nunique()),
    "max_spectral_count_only_M": int(pres_out.loc[pres_out["all_in_group"] == "M only", "spectral_count"].max()),
    "max_spectral_count_only_WT": int(pres_out.loc[pres_out["all_in_group"] == "WT only", "spectral_count"].max()),
}

# features detected in >=1 sample of each group but not testable (<2 per group)
sub = occ[(occ.n_WT >= 1) & (occ.n_M >= 1) & ((occ.n_WT < 2) | (occ.n_M < 2))]
summary_pres["partial_both_groups_untested"] = int(len(sub))

out["presence_only"] = summary_pres

# --- unnamed numeric shifts ----------------------------------------------
SHIFT = re.compile(r"\[([+-]\d+(?:\.\d+)?)\]")


def shifts(s):
    return [float(x) for x in SHIFT.findall(str(s))]


rows = []
for fid, r in feat.iterrows():
    vals = shifts(r.get("proteoform", ""))
    if not vals:
        continue
    for v in vals:
        rows.append({
            "featureId": fid, "shift": v, "shift_key": f"{v:+.2f}",
            "gene": r.get("gene"), "uniprot_acc": r.get("uniprot_acc"),
            "n_unexpected_mods": r.get("n_unexpected_mods"),
        })
sh = pd.DataFrame(rows)

# attach statistics where the feature was tested
stat_cols = ["log2FoldChange", "bhFdr", "pValue", "spectralCount"]
sh = sh.merge(cand[stat_cols], left_on="featureId", right_index=True, how="left")
sh = sh.merge(occ[["n_WT", "n_M"]], left_on="featureId", right_index=True, how="left")
sh = sh.merge(tiers[["tier"]], left_on="featureId", right_index=True, how="left")
sh["significant"] = (sh["bhFdr"] < 0.05) & (sh["log2FoldChange"].abs() > 1)
sh["direction"] = [
    "untested" if pd.isna(fc) else ("up in M" if fc > 0 else "up in WT")
    for fc in sh["log2FoldChange"]
]
sh.to_csv("shift_feature_rows.csv", index=False)

fam = (sh.groupby("shift_key")
         .agg(n_features=("featureId", "nunique"),
              n_genes=("gene", "nunique"),
              n_sig=("significant", "sum"),
              n_up_M=("direction", lambda s: (s == "up in M").sum()),
              n_up_WT=("direction", lambda s: (s == "up in WT").sum()),
              max_abs_fc=("log2FoldChange", lambda s: s.abs().max()),
              total_spectral=("spectralCount", "sum"),
              genes=("gene", lambda s: ";".join(sorted(set(map(str, s.dropna()))))[:160]))
         .sort_values(["n_features", "n_genes"], ascending=False))
fam.to_csv("shift_families.csv")

# complement analysis: does +v pair with -v (or with -v + constant)?
vals = sorted(sh["shift"].unique())
comps = []
for v in vals:
    if v <= 0:
        continue
    exact = any(abs(x + v) < 0.01 for x in vals)
    near = sorted([x for x in vals if x < 0], key=lambda x: abs(x + v))[:2]
    comps.append({"shift": v, "n_features": int((sh["shift"] == v).sum()),
                  "exact_complement": exact,
                  "nearest_negative": near})
out["shift_complements"] = comps

out["shift_families_top"] = fam.head(15).reset_index().to_dict("records")
out["shift_overview"] = {
    "features_with_numeric_shift": int(sh["featureId"].nunique()),
    "distinct_shift_values": int(sh["shift"].nunique()),
    "families_recurring_multi_protein": int(((fam.n_genes >= 2) & (fam.n_features >= 2)).sum()),
    "features_with_shift_and_sig": int(sh.loc[sh.significant, "featureId"].nunique()),
}

# --- effect-size view across ALL tiers -----------------------------------
top = cand.copy()
top["abs_log2FC"] = top["log2FoldChange"].abs()
top = top.join(occ[["n_WT", "n_M"]])
top = top.join(tiers[["tier"]])
top = top.drop(columns=[c for c in ("proteoform",) if c in top.columns]).join(
    feat[["proteoform", "n_unexpected_mods", "proteoform_level", "mass"]])
top.sort_values("abs_log2FC", ascending=False).head(60).to_csv("top_effect_all_tiers.csv")
out["top_effect_head"] = top.sort_values("abs_log2FC", ascending=False).head(15)[
    ["gene", "log2FoldChange", "bhFdr", "n_WT", "n_M", "tier", "spectralCount",
     "n_unexpected_mods"]].reset_index().to_dict("records")

# --- tier composition ----------------------------------------------------
tc = (tiers.groupby(["tier", "direction"]).size().rename("n").reset_index())
tc.to_csv("tier_direction_counts.csv", index=False)
out["tier_direction"] = tc.to_dict("records")

# --- gene-level concentration -------------------------------------------
sig = tiers[tiers["significant"] == True]
gene = (sig.groupby("gene")
        .agg(n_sig=("tier", "size"),
             n_up_M=("direction", lambda s: (s == "up in M").sum()),
             n_up_WT=("direction", lambda s: (s == "up in WT").sum()),
             med_fc=("log2FoldChange", "median"),
             max_abs_fc=("log2FoldChange", lambda s: s.abs().max()))
        .sort_values("n_sig", ascending=False))
gene.to_csv("gene_significant_counts.csv")
out["top_genes"] = gene.head(12).reset_index().to_dict("records")
out["tiers_columns"] = list(tiers.columns)
out["n_genes_with_sig"] = int(gene.shape[0])

json.dump(out, open("summary.json", "w"), indent=1, default=str)
print(json.dumps({k: out[k] for k in
                  ["sample_detection", "presence_only", "shift_overview",
                   "top_effect_head", "top_genes"]}, indent=1, default=str)[:6000])
