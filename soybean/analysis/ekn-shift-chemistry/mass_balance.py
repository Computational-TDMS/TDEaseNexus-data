# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Close the EKN mass balance per source row against the dual-acylation model.

Backbones and the model mass come from `parse-proforma` cells; residuals are
read from the TopPIC rows (ekn_basis.csv).  Nothing here is a measurement of
chemistry: it is the arithmetic that has to close before spectra are read.
"""

import csv
import json
import subprocess

MODEL = ("[Myristoyl]-GN[Palmitoyl]CGSNPKTNEGPEAVPEPVIEEVKVEQKESSEANVETKLK"
         "DQTPNDSEIKSLGTLLNEKVEEAPKTEEVTAEPKPEEPKTDEVKVQEEKPKAEEAKAET")


def mass(proforma):
    out = subprocess.run(["tdeasenexus-toolbox", "parse-proforma", "--proforma", proforma],
                         capture_output=True, text=True).stdout
    return json.loads(out)["result"]["monoisotopicMass"]


model_mass = mass(MODEL)
myr = mass("[Myristoyl]-A") - mass("A")
palm = mass("AC[Palmitoyl]A") - mass("ACA") if False else mass("C[Palmitoyl]") - mass("C")
print(f"model (NME + myristoyl-G1 + palmitoyl-C3) = {model_mass:.5f}")
print(f"myristoyl {myr:.5f}  palmitoyl {palm:.5f}  sum {myr + palm:.5f}")

met_res = mass("AM") - mass("A")
ace = mass("[Acetyl]-A") - mass("A")
FRAME_OFFSET = {"NONE": met_res, "NME": 0.0, "NME_ACETYLATION": ace}
print(f"frame offsets: Met {met_res:.5f}, acetyl {ace:.5f}")

rows = list(csv.DictReader(open("ekn_basis.csv", encoding="utf-8")))
keep = [r for r in rows if float(r["residual"]) > 300 or r["nform"] == "NME_ACETYLATION"]
print(f"\n{'span':>7} {'nform':>17} {'observed':>11} {'bare':>11} {'residual':>10} "
      f"{'model delta':>12} {'diff':>8}")
out = []
for r in keep:
    obs, bare, res = float(r["raw_precursor"]), float(r["bare"]), float(r["residual"])
    # expected residual in this row's frame: components minus the processing
    # already modelled by the frame (span-independent, so truncations drop out)
    expect = myr + palm - FRAME_OFFSET[r["nform"]]
    out.append({"span": r["span"], "nform": r["nform"], "observed_raw_precursor": obs,
                "bare_backbone": bare, "reported_residual": res,
                "model_residual_in_frame": round(expect, 4),
                "diff_Da": round(res - expect, 4), "file": r["file"]})
    print(f"{r['span']:>7} {r['nform']:>17} {obs:>11.4f} {bare:>11.4f} {res:>10.4f} "
          f"{expect:>12.4f} {res - expect:>+8.4f}")

with open("ekn_mass_balance.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)
print(f"\nwrote ekn_mass_balance.csv ({len(out)} rows)")
