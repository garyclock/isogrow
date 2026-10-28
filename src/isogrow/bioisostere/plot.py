"""Property-space plot: ΔtPSA vs ΔlogP (SwissBioisostere's interactive plot, static here)."""
from __future__ import annotations

import html
import os

import pandas as pd


def _color(score: float | None) -> str:
    if score is None:
        return "#95a5a6"
    if score >= 0.3:
        return "#c0392b"      # clearly better
    if score > 0:
        return "#e67e22"      # mildly better
    if score == 0:
        return "#7f8c8d"
    return "#2c7bb6"          # worse


def plot_property_space(
    df: pd.DataFrame,
    out_html: str,
    x: str = "Δ tPSA",
    y: str = "Δ logP",
    title: str = "Property space (ΔtPSA vs ΔlogP)",
    top_table: int = 60,
) -> str:
    """Render a self-contained HTML with an SVG scatter (point size ~ frequency)."""
    d = df.dropna(subset=[x, y]).copy()
    if d.empty:
        raise ValueError("no rows with both coordinates")

    W, H, m = 720, 520, 64
    xs, ys = list(d[x]), list(d[y])
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    xpad = (xmax - xmin) * 0.08 or 1
    ypad = (ymax - ymin) * 0.08 or 1
    xmin, xmax = xmin - xpad, xmax + xpad
    ymin, ymax = ymin - ypad, ymax + ypad

    def sx(v):
        return m + (v - xmin) / (xmax - xmin) * (W - 2 * m)

    def sy(v):
        return H - m - (v - ymin) / (ymax - ymin) * (H - 2 * m)

    pts = []
    freq_max = max(1, d["Frequency"].max()) if "Frequency" in d else 1
    for _, r in d.iterrows():
        cx, cy = sx(r[x]), sy(r[y])
        f = r.get("Frequency", 1) or 1
        rad = 3 + 6 * (float(f) / float(freq_max)) ** 0.5
        c = _color(r.get("success_score"))
        tip = f"{r['Candidate Fragments']} | freq={f} | Δpchembl={r.get('Δ pchembl')} | success={r.get('success_score')}"
        pts.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{rad:.1f}" fill="{c}" fill-opacity="0.75" '
                   f'stroke="#333" stroke-width="0.4"><title>{html.escape(tip)}</title></circle>')

    # axes
    def xticks():
        out = []
        for i in range(6):
            v = xmin + (xmax - xmin) * i / 5
            px = sx(v)
            out.append(f'<line x1="{px:.1f}" y1="{H-m}" x2="{px:.1f}" y2="{H-m+5}" stroke="#666"/>'
                       f'<text x="{px:.1f}" y="{H-m+20}" font-size="11" text-anchor="middle" fill="#333">{v:.1f}</text>')
        return "".join(out)

    def yticks():
        out = []
        for i in range(6):
            v = ymin + (ymax - ymin) * i / 5
            py = sy(v)
            out.append(f'<line x1="{m-5}" y1="{py:.1f}" x2="{m}" y2="{py:.1f}" stroke="#666"/>'
                       f'<text x="{m-8}" y="{py+4:.1f}" font-size="11" text-anchor="end" fill="#333">{v:.1f}</text>')
        return "".join(out)

    svg = f"""<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px;background:#fff">
<rect x="{m}" y="{m}" width="{W-2*m}" height="{H-2*m}" fill="#fbfcfd" stroke="#d0d7de"/>
<line x1="{m}" y1="{H-m}" x2="{W-m}" y2="{H-m}" stroke="#333"/>
<line x1="{m}" y1="{m}" x2="{m}" y2="{H-m}" stroke="#333"/>
{xticks()}{yticks()}
<text x="{(W)/2}" y="{H-14}" font-size="13" text-anchor="middle" fill="#333">{html.escape(x)}</text>
<text x="18" y="{H/2}" font-size="13" text-anchor="middle" fill="#333" transform="rotate(-90 18 {H/2})">{html.escape(y)}</text>
{''.join(pts)}
</svg>"""

    cols = [c for c in ("Candidate Fragments", "Frequency", "# Better", "# Similar", "# Worse",
                        "Δ logP", "Δ tPSA", "Δ MW", "Δ pchembl", "success_score") if c in d.columns]
    trs = []
    for _, r in d.head(top_table).iterrows():
        trs.append("<tr>" + "".join(f"<td>{html.escape(str(r[c]))}</td>" for c in cols) + "</tr>")
    ths = "".join(f"<th>{html.escape(c)}</th>" for c in cols)

    doc = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>body{{font-family:system-ui,"Microsoft YaHei",sans-serif;max-width:960px;margin:28px auto;padding:0 18px;color:#1f2328}}
h1{{font-size:22px}} table{{border-collapse:collapse;width:100%;font-size:12.5px;margin-top:16px}}
th,td{{border:1px solid #d0d7de;padding:5px 8px;text-align:left}} th{{background:#f6f8fa}}
tr:nth-child(even) td{{background:#fbfcfd}} code{{font-family:ui-monospace,Consolas,monospace}}
.legend span{{display:inline-block;margin-right:14px;font-size:12px}}</style></head><body>
<h1>{html.escape(title)}</h1>
<p style="color:#57606a;font-size:13px">点大小≈频次；颜色：<span style="color:#c0392b">红=更好</span>
<span style="color:#e67e22">橙=略好</span> <span style="color:#7f8c8d">灰=相似</span> <span style="color:#2c7bb6">蓝=更差</span>（悬停看详情）</p>
{svg}
<table><tr>{ths}</tr>{''.join(trs)}</table>
</body></html>"""
    os.makedirs(os.path.dirname(out_html) or ".", exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(doc)
    return out_html
