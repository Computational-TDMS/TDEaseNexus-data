# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Extract the EKN (Q2KMJ4) object rows from the soybean-all TopPIC tables.

Writes ekn_proteoforms.csv (one row per proteoform row in the 10 proteoform
tables) and ekn_prsms.csv (every PrSM row for Q2KMJ4 with spectrum/scan ids).
Reads only; source tables stay untouched.
"""

import csv
import glob
import os

ACCESSIONS = ("Q2KMJ4", "C6T4V4")  # EKN and its uncharacterised N-terminal homolog
RUN = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def read_table(path):
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("Data file name"))
    return csv.DictReader(lines[start:], delimiter="\t")


def keep(row):
    acc = row.get("Protein accession", "")
    return any(a in acc for a in ACCESSIONS)


proto_cols = ["data_file", "Prsm ID", "Proteoform ID", "Spectrum ID", "Scan(s)",
              "Retention time", "Charge", "Precursor mass", "Adjusted precursor mass",
              "Feature ID", "Protein accession", "First residue", "Last residue",
              "Previous amino acid", "Next amino acid", "Proteoform",
              "Proteoform mass", "Protein N-terminal form", "Fixed PTMs",
              "# Unexpected modifications", "Unexpected modifications",
              "# Matched masses", "# Matched fragments", "E-value",
              "Spectrum-level Q-value", "Proteoform-level Q-value", "Database protein sequence"]

prsm_cols = ["data_file", "Prsm ID", "Spectrum ID", "Scan(s)", "Retention time",
             "Charge", "Precursor mass", "Adjusted precursor mass", "Proteoform ID",
             "Protein accession", "First residue", "Last residue",
             "Previous amino acid", "Next amino acid", "Proteoform", "Proteoform mass",
             "Protein N-terminal form", "# Unexpected modifications",
             "Unexpected modifications", "# Masses", "# Matched masses",
             "# Matched fragments", "MIScore", "E-value", "Spectrum-level Q-value",
             "Proteoform-level Q-value"]

out_p, out_s = [], []
for path in sorted(glob.glob(os.path.join(RUN, "*_ms2_toppic_proteoform.tsv"))):
    for row in read_table(path):
        if keep(row):
            row["data_file"] = os.path.basename(path).split("_ms2_toppic")[0]
            out_p.append({c: row.get(c, "") for c in proto_cols})
for path in sorted(glob.glob(os.path.join(RUN, "*_ms2_toppic_prsm.tsv"))):
    for row in read_table(path):
        if keep(row):
            row["data_file"] = os.path.basename(path).split("_ms2_toppic")[0]
            out_s.append({c: row.get(c, "") for c in prsm_cols})


def write(path, rows, cols):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(path, len(rows))


write("ekn_proteoforms.csv", out_p, proto_cols)
write("ekn_prsms.csv", out_s, prsm_cols)

# sequence from the source FASTA
seq = None
with open(os.path.join(RUN, "Glycine.fasta"), encoding="utf-8", errors="replace") as fh:
    hit = False
    for line in fh:
        if line.startswith(">"):
            hit = "Q2KMJ4" in line
        elif hit:
            seq = (seq or "") + line.strip()
print("EKN length", len(seq), "seq", seq)
open("ekn_sequence.txt", "w").write(seq + "\n")
