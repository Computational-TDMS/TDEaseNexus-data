# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Basis table: observed vs bare-backbone mass per EKN PrSM row.

bare = Proteoform mass - reported residual (the unmodified backbone of the
reported span/terminal form).  raw_delta = Precursor mass - bare.
This is the residual recomputed on the *raw* observed precursor, before
TopPIC's isotope adjustment, so 1-Da adjustments become visible.
"""

import csv
import statistics

rows = [r for r in csv.DictReader(open("ekn_proteoforms.csv", encoding="utf-8"))
        if "Q2KMJ4" in r["Protein accession"]]


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


out = []
for r in rows:
    pm, prec = f(r["Proteoform mass"]), f(r["Precursor mass"])
    un = r["Unexpected modifications"] or ""
    shift = float(un.split(":")[0]) if un.startswith(("+", "-")) else 0.0
    bare = pm - shift
    out.append({
        "file": r["data_file"][-3:], "span": f"{r['First residue']}-{r['Last residue']}",
        "nform": r["Protein N-terminal form"], "residual": shift,
        "bare": round(bare, 4), "raw_precursor": prec,
        "raw_delta": round(prec - bare, 4),
        "adjusted_delta": round(f(r["Adjusted precursor mass"]) - bare, 4),
        "charge": r["Charge"], "fragments": r["# Matched fragments"],
        "scan": r["Scan(s)"], "spectrum": r["Spectrum ID"], "q_proteoform": r["Proteoform-level Q-value"],
    })

with open("ekn_basis.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)

none_101 = [o for o in out if o["span"] == "1-101" and o["nform"] == "NONE"]
allnone = [o for o in out if o["nform"] == "NONE" and o["residual"] > 300]
print(f"{'file':>4} {'span':>7} {'nform':>17} {'resid':>10} {'bare':>11} {'raw':>11} {'rawDelta':>9} {'adjDelta':>9} frag")
for o in sorted(out, key=lambda o: (o["nform"], o["span"], o["residual"])):
    print(f"{o['file']:>4} {o['span']:>7} {o['nform']:>17} {o['residual']:>10.4f} {o['bare']:>11.4f} "
          f"{o['raw_precursor']:>11.4f} {o['raw_delta']:>9.4f} {o['adjusted_delta']:>9.4f} {o['fragments']:>4}")

vals = [o["raw_delta"] for o in allnone]
adj = [o["adjusted_delta"] for o in allnone]
print(f"\nNONE-frame +300-family rows: n={len(vals)} raw mean={statistics.mean(vals):.4f} "
      f"sd={statistics.pstdev(vals):.4f} min={min(vals):.4f} max={max(vals):.4f}")
print(f"adjusted mean={statistics.mean(adj):.4f} sd={statistics.pstdev(adj):.4f} "
      f"min={min(adj):.4f} max={max(adj):.4f}")
