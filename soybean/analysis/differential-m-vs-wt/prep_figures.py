# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas"]
# ///
"""Prepare chart-forge data payloads for the WT-vs-M figures.

Sources: comparison_WT_vs_M_full.json (toppic_plot.py compare), tiers.csv
(toppic_shortlist.py), features.csv / expression.csv (toppic_run.py).
Writes figs/*.json plus figs/manifest.json recording sources and transforms.
"""

import json
import math
import re
from pathlib import Path

import pandas as pd

OUT = Path("figs")
OUT.mkdir(exist_ok=True)
SHIFT = re.compile(r"\[([+-]\d+(?:\.\d+)?)\]")

t = pd.read_csv("tiers.csv", index_col=0)
feat = pd.read_csv("features.csv", index_col=0)
feat = feat[~feat.index.duplicated()]
t["tier"] = "not significant"
sig = t.significant == True
t.loc[sig & (t.min_det >= 4), "tier"] = "tier A"
t.loc[sig & (t.min_det == 3), "tier"] = "tier B"
t.loc[sig & (t.min_det < 3), "tier"] = "tier C"
t["shifts"] = [SHIFT.findall(str(p)) for p in t["proteoform"]]
t["has_shift"] = [len(s) > 0 for s in t["shifts"]]
t["is317"] = [any(abs(float(v) - 317.39) <= 0.05 for v in ss) for ss in t["shifts"]]
t["direction"] = ["up in M" if fc > 0 else "up in WT" for fc in t.log2FoldChange]

# ---------- figure 1: volcano with evidence tiers -------------------------
series_order = ["tier A", "tier B", "tier C", "not significant", "+317.39 family"]
points = {k: [] for k in series_order}
for fid, r in t.iterrows():
    key = "+317.39 family" if r.is317 else r.tier
    rows = feat.loc[fid]
    span = ""
    if r.shifts:
        span = "; ".join(f"{float(s):+.3f}" for s in r.shifts)
    points[key].append({
        "x": round(float(r.log2FoldChange), 3),
        "y": round(-math.log10(float(r.bhFdr)), 3),
        "hover": (f"<b>{r.gene}</b> · {fid}<br>log2FC {r.log2FoldChange:+.2f} "
                  f"({'higher in M' if r.log2FoldChange > 0 else 'higher in WT'})"
                  f"<br>BH FDR {r.bhFdr:.2e}<br>detected WT {r.n_WT}/5 · M {r.n_M}/5"
                  f"<br>spectral count {int(r.spectralCount)}"
                  + (f"<br>numeric shift {span} Da" if span else "")),
    })
volcano = {
    "data": {
        "series": [{"label": f"{k} ({len(points[k])})", "points": points[k]}
                   for k in series_order],
        "vlines": [{"v": -1, "label": "log2FC -1"}, {"v": 1, "label": "log2FC +1"}],
        "hlines": [{"v": round(-math.log10(0.05), 3), "label": "BH FDR 0.05"}],
        "xlabel": "log2 fold change  (M - WT)",
        "ylabel": "-log10 BH-adjusted FDR",
        "legend": True, "height": 560, "pointR": 2.9,
    },
    "meta": {
        "title": "WT vs M: 590 testable proteoform features, coloured by detection tier",
        "source": "comparison_WT_vs_M_full.json + tiers.csv (toppicr-py)",
        "note": ("Only features with >=2 detections in BOTH groups are testable; "
                 "presence-only features are outside this plot."),
    },
}
json.dump(volcano, open(OUT / "volcano_tiers.json", "w"), indent=1)

# ---------- figure 2: proteoform-state split ------------------------------
def split(df):
    g = df[df.significant == True]
    return {
        "shifted · up in M": int(((g.has_shift) & (g.direction == "up in M")).sum()),
        "shifted · up in WT": int(((g.has_shift) & (g.direction == "up in WT")).sum()),
        "no numeric shift · up in M": int(((~g.has_shift) & (g.direction == "up in M")).sum()),
        "no numeric shift · up in WT": int(((~g.has_shift) & (g.direction == "up in WT")).sum()),
    }


groups = {
    "LBA (leghemoglobin A)": split(t[t.gene == "LBA"]),
    "LGB1 (leghemoglobin C3)": split(t[t.gene == "LGB1"]),
    "all genes": split(t),
}
series = ["shifted · up in M", "shifted · up in WT", "no numeric shift · up in M",
          "no numeric shift · up in WT"]
rows = []
for label, d in groups.items():
    tot = sum(d.values())
    for s in series:
        rows.append({"label": f"{label.split(' (')[0]} · {tot} sig", "value": d[s], "series": s,
                     "hover": f"<b>{label}</b> · {s}: {d[s]} significant features"})
states = {
    "data": {"rows": rows, "series": series,
             "xlabel": "gene group (significant features, BH FDR < 0.05 and |log2FC| > 1)",
             "ylabel": "significant proteoform features",
             "valueLabels": True, "height": 420},
    "meta": {
        "title": "Within one gene, shifted and unshifted proteoforms move in opposite directions",
        "source": "tiers.csv + features.csv (toppicr-py)",
        "note": ("'shifted' = the proteoform string carries a numeric bracket mass shift "
                 "(chemistry unresolved); 'no numeric shift' may still carry a named "
                 "terminal modification."),
    },
}
json.dump(states, open(OUT / "proteoform_states.json", "w"), indent=1)

# ---------- figure 3: recurring shift families ---------------------------
fam_rows = []
for fid, r in t.iterrows():
    for v in r.shifts:
        fam_rows.append({"featureId": fid, "v": float(v), "gene": r.gene,
                         "direction": r.direction, "sig": bool(r.significant),
                         "fc": r.log2FoldChange})
sh = pd.DataFrame(fam_rows)
centers = []
sh["family"] = None
for v in sorted(sh.v.unique()):
    hit = next((c for c in centers if abs(c - v) <= 0.05), None)
    if hit is None:
        centers.append(v)
        hit = v
    sh.loc[sh.v == v, "family"] = hit
fam = (sh.groupby("family")
       .agg(n_feat=("featureId", "nunique"), n_gene=("gene", "nunique"),
            n_sig=("sig", "sum"),
            upM=("direction", lambda s: int((s == "up in M").sum())),
            upWT=("direction", lambda s: int((s == "up in WT").sum())),
            maxfc=("fc", lambda s: float(s.abs().max())))
       .query("n_feat >= 3").sort_values(["n_feat", "n_gene"], ascending=False))
fam.to_csv("shift_families_tol005.csv")
sel = fam.head(9)
fam_rows_out = []
for f, r in sel.iterrows():
    lab = f"{f:+.2f}"
    fam_rows_out.append({"label": f"{lab} · {int(r.n_gene)} prot", "value": int(r.upM),
                         "series": "up in M", "hover": f"<b>{lab}</b> · up in M: {int(r.upM)}"})
    fam_rows_out.append({"label": f"{lab} · {int(r.n_gene)} prot", "value": int(r.upWT),
                         "series": "up in WT", "hover": f"<b>{lab}</b> · up in WT: {int(r.upWT)}"})
families = {
    "data": {"rows": fam_rows_out, "series": ["up in M", "up in WT"],
             "xlabel": "unresolved numeric shift family (mass, proteins carrying it)",
             "ylabel": "features carrying that shift",
             "valueLabels": True, "height": 420},
    "meta": {"title": "Recurring unresolved shift families across unrelated proteins",
             "source": "features.csv proteoform strings, clustered at 0.05 Da",
             "note": "Counts are all tested features; significance is not implied."},
}
json.dump(families, open(OUT / "shift_families.json", "w"), indent=1)

manifest = {
    "sources": {
        "expression.csv": "toppic_run.py build (1794 features x 10 samples)",
        "features.csv": "toppic_run.py build",
        "comparison_WT_vs_M_full.json": "toppic_plot.py compare --group-a WT --group-b M --top 5000",
        "tiers.csv": "toppic_shortlist.py --groups WT,M",
    },
    "tests": {"fdr": 0.05, "abs_log2fc": 1.0, "min_detected_per_group": 2,
              "tier_a": "min detections >= 4 in both groups",
              "tier_b": "min detections == 3", "tier_c": "min detections <= 2"},
    "counts": {"tested": int(len(t)), "significant": int((t.significant == True).sum())},
}
json.dump(manifest, open(OUT / "manifest.json", "w"), indent=1)
print("written", [p.name for p in OUT.iterdir()])
print(json.dumps({k: sum(r["value"] for r in rows if r["label"].startswith(k.splitlines()[0])) for k in []}))
print(fam.head(10).to_string())
