# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Separate genuine discriminators from degenerate ions in the two EKN cards.

A row is a genuine discriminator only when the two models predict materially
different m/z for that ion (|delta| > 0.05 Da).  Rows whose predictions differ
by <0.01 Da are degenerate: a matched/rejected split there is a knife-edge of
the correlation gate, not evidence.  Numbers are read from card cells.
"""

import json


def analyse(path):
    c = json.load(open(path))
    print(f"\n===== {path} =====")
    print("models:", [m["id"] for m in c["models"]])
    print("compositionDeltaDa:", c["verdict"]["compositionDeltaDa"])
    genuine, degenerate = [], []
    for row in c["verdict"]["ions"]:
        cells = row.get("cellsByScan", {}).get("10876") or row.get("cells") or {}
        a, b = cells.get("dual_acyl"), cells.get("mass_at_nterm")
        if not a or not b:
            continue
        d = abs(a["theoreticalMz"] - b["theoreticalMz"])
        entry = (row["label"], row["role"], d, sorted(row["matchedBy"]),
                 a.get("state"), b.get("state"),
                 a.get("errorPpm"), b.get("errorPpm"))
        (genuine if d > 0.05 else degenerate).append(entry)
    print(f"rows: {len(genuine)} genuine (|dmz|>0.05), {len(degenerate)} degenerate")
    print("\n-- genuine discriminators --")
    print(f"{'ion':>8} {'role':<16} {'|dmz|':>9} {'matched':<24} {'dual':<18} {'alt':<18}")
    for lab, role, d, mb, sa, sb, ea, eb in sorted(genuine, key=lambda x: -x[2]):
        print(f"{lab:>8} {role:<16} {d:>9.3f} {','.join(mb):<24} "
              f"{sa+' '+format(ea or 0,'.1f'):<18} {sb+' '+format(eb or 0,'.1f'):<18}")
    print("\n-- degenerate examples --")
    for lab, role, d, mb, sa, sb, ea, eb in sorted(degenerate, key=lambda x: -x[2])[:6]:
        print(f"{lab:>8} {role:<16} {d:>9.5f} {','.join(mb):<24} {sa:<18} {sb:<18}")


for p in ("card_scan11831_shared.json", "card_full_shared.json"):
    analyse(p)
