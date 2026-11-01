# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas"]
# ///
import pandas as pd

f = pd.read_csv("shift_families_tol005.csv")
print("families >=3 features:", len(f))
print("families >=3 features & >=2 proteins:", int(((f.n_feat >= 3) & (f.n_gene >= 2)).sum()))
print("features in >=2-protein families:", int(f.loc[f.n_gene >= 2, "n_feat"].sum()))
t = pd.read_csv("tiers.csv", index_col=0)
sig = t[t.significant == True]
print("sig median |FC|:", round(sig.log2FoldChange.abs().median(), 2))
top = sig.reindex(sig.log2FoldChange.abs().sort_values(ascending=False).index)
print("top20 |FC| tiers:", top.head(20).min_det.apply(
    lambda d: "A" if d >= 4 else ("B" if d == 3 else "C")).value_counts().to_dict())
print("sig with |FC|>2:", int((sig.log2FoldChange.abs() > 2).sum()))
