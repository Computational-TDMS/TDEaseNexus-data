# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Make a small inline-viewable copy of the automatic qc_report.html.

The genuine report embeds all 138,622 TopFD MS1 features (20.3 MB of its
30 MB data payload). For an inline copy, the Feature map tab is drawn from an
evenly strided 4,000-feature subsample; every other panel keeps its complete
data. The original report is untouched.
"""

import json
import re

SRC = "qc_report.html"
DST = "qc_report_inline.html"
KEEP = 4000

html = open(SRC, encoding="utf-8").read()
start = html.index("window.DATAByPanel = ") + len("window.DATAByPanel = ")
depth = 0
for i in range(start, len(html)):
    if html[i] == "{":
        depth += 1
    elif html[i] == "}":
        depth -= 1
        if depth == 0:
            end = i + 1
            break
panels = json.loads(html[start:end])

fm = panels["featuremap"]["data"]
feats = fm["features"]
stride = max(1, len(feats) // KEEP)
fm["features"] = feats[::stride]
fm["shownInInlineCopy"] = len(fm["features"])
fm["note"] = f"inline copy: {len(fm['features'])} of {len(feats)} bars shown"
print(f"featuremap: {len(feats)} -> {len(fm['features'])} (stride {stride})")

new_json = json.dumps(panels, separators=(",", ":"), ensure_ascii=False)
html = html[:start] + new_json + html[end:]

caption = ("All TopFD MS1 features are shown as horizontal elution-window bars at average m/z; "
           "bar width is the observed RT interval, colour is apex intensity, and charge is available on hover.")
new_caption = (caption + f" <strong>Inline copy:</strong> {len(fm['features'])} of {len(feats)} "
               "bars drawn (even stride); the full report carries all of them.")
assert caption in html
html = html.replace(caption, new_caption)

banner = (
    '<div style="margin:0 0 14px;padding:10px 14px;border-left:4px solid #b45309;'
    'background:#fffbeb;font:13px/1.45 system-ui,sans-serif;color:#3f2d0b">'
    "<b>Inline copy</b> of <code>qc_report.html</code> from this analysis directory "
    "(toppicr-py build + volcano groups <code>WT0[1-5]</code> vs <code>M0[1-5]</code>). "
    "Only the Feature map tab is subsampled for size; every other panel is complete."
    "</div>"
)
html = html.replace("<body>", "<body>\n" + banner, 1)
open(DST, "w", encoding="utf-8").write(html)
print(f"wrote {DST}: {len(html)/1e6:.2f} MB")
