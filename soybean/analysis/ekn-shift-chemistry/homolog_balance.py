# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Same mass balance for the EKN homolog C6T4V4 (100527608).

Arithmetic only: does the dual-acylation model close on the homolog's rows too?
No spectrum is measured here.
"""

import csv
import json
import subprocess

MODEL_RESIDUAL = {"NONE": None, "NME": None, "NME_ACETYLATION": None}


def mass(proforma):
    out = subprocess.run(["tdeasenexus-toolbox", "parse-proforma", "--proforma", proforma],
                         capture_output=True, text=True).stdout
    return json.loads(out)["result"]["monoisotopicMass"]


myr = mass("[Myristoyl]-A") - mass("A")
palm = mass("C[Palmitoyl]") - mass("C")
met = mass("AM") - mass("A")
ace = mass("[Acetyl]-A") - mass("A")
MODEL_RESIDUAL = {"NONE": myr + palm - met, "NME": myr + palm, "NME_ACETYLATION": myr + palm - ace}

rows = [r for r in csv.DictReader(open("ekn_proteoforms.csv", encoding="utf-8"))
        if "C6T4V4" in r["Protein accession"]]
print(f"{'file':>4} {'span':>7} {'nform':>17} {'span (Da)':>21} {'residual':>10} {'model':>9} {'diff':>8}")
kept = []
for r in rows:
    try:
        res = float((r["Unexpected modifications"] or "0").split(":")[0])
    except ValueError:
        res = 0.0
    if res < 300:
        continue
    expect = MODEL_RESIDUAL[r["Protein N-terminal form"]]
    kept.append((r["data_file"][-3:], r["First residue"] + "-" + r["Last residue"],
                 r["Protein N-terminal form"], r["Unexpected modifications"], res, expect, res - expect))
for k in sorted(kept, key=lambda x: abs(x[6])):
    print(f"{k[0]:>4} {k[1]:>7} {k[2]:>17} {k[3]:>21} {k[4]:>10.4f} {k[5]:>9.4f} {k[6]:>+8.4f}")
close = [k for k in kept if abs(k[6]) < 0.02]
print(f"\nrows within 0.02 Da of the model: {len(close)} of {len(kept)}")
