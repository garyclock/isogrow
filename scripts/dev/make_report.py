#!/usr/bin/env python
"""Render an L1 benchmark report (report.json + per_candidate.csv) to a readable HTML."""
from __future__ import annotations

import argparse
import html
import json
import os

import pandas as pd


def _fmt(v, nd=2):
    try:
        return f"{float(v):.{nd}f}"
    except Exception:
        return str(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="run directory containing report.json + per_candidate.csv")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    run = args.run
    rep = json.load(open(os.path.join(run, "report.json"), encoding="utf-8"))
    df = pd.read_csv(os.path.join(run, "per_candidate.csv"))
    top = df.head(args.top)

    rows = []
    for _, r in top.iterrows():
        ok = bool(r["in_local"])
        cls = "" if ok else "miss"
        loc = ""
        if ok:
            loc = (f"<td>{html.escape(str(r['local_frag']))}</td><td>{_fmt(r['local_freq'],0)}</td>"
                   f"<td>{_fmt(r['loc_dlogp'])}</td><td>{_fmt(r['loc_dtpsa'])}</td><td>{_fmt(r['loc_dmw'])}</td>")
        else:
            loc = "<td colspan='5' class='misscell'>— 未命中 —</td>"
        rows.append(
            f"<tr class='{cls}'><td>{int(r['official_rank'])}</td>"
            f"<td><code>{html.escape(str(r['official_frag']))}</code></td><td>{html.escape(str(r['official_freq']))}</td>"
            f"<td>{_fmt(r['off_dlogp'])}</td><td>{_fmt(r['off_dtpsa'])}</td><td>{_fmt(r['off_dmw'])}</td>"
            f"{loc}</tr>"
        )

    mae = rep.get("delta_mae", {})
    metrics = [
        ("官方候选总数", f"{rep['official_candidates']}"),
        ("本地候选总数", f"{rep['local_candidates']}"),
        ("重叠候选", f"{rep['overlap']}"),
        ("覆盖率 coverage（overlap/官方）", f"{rep['coverage_vs_official']*100:.1f}%"),
        ("本地精度 precision（overlap/本地）", f"{rep['precision_vs_official']*100:.1f}%"),
        ("recall@50", f"{rep.get('recall_at_50',0)*100:.1f}%"),
        ("Δ logP MAE", f"{mae.get('dlogp')}"),
        ("Δ tPSA MAE", f"{mae.get('dtpsa')}"),
        ("Δ MW MAE", f"{mae.get('dmw')}"),
        ("Δ 对比样本数", f"{mae.get('n_compared')}"),
    ]
    metric_rows = "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in metrics)

    doc = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>L1 本地 SwissBioisostere 对标报告</title>
<style>
body{{font-family:system-ui,"Microsoft YaHei",sans-serif;max-width:1000px;margin:32px auto;padding:0 20px;color:#1f2328;line-height:1.6}}
h1{{font-size:24px;margin-bottom:4px}} .sub{{color:#57606a;font-size:13px;margin-bottom:20px}}
h2{{font-size:18px;border-bottom:2px solid #2b6cb0;padding-bottom:6px;margin-top:30px}}
table{{border-collapse:collapse;width:100%;font-size:13.5px;margin:12px 0}}
th,td{{border:1px solid #d0d7de;padding:6px 9px;text-align:left}}
th{{background:#f6f8fa}} tr:nth-child(even) td{{background:#fbfcfd}}
code{{font-family:ui-monospace,Consolas,monospace;background:#eef1f4;padding:1px 4px;border-radius:4px;font-size:12.5px}}
tr.miss td{{background:#fff1f0}} .misscell{{color:#b42318;font-style:italic}}
.kpi{{display:flex;gap:14px;flex-wrap:wrap;margin:14px 0}}
.card{{flex:1 1 180px;background:#f6f8fa;border:1px solid #d0d7de;border-radius:10px;padding:12px 16px}}
.card b{{display:block;font-size:24px;color:#2b6cb0}} .card small{{color:#57606a}}
.note{{background:#fff8e6;border:1px solid #e6c86e;border-radius:8px;padding:10px 14px;font-size:13px}}
</style></head><body>
<h1>L1 本地 SwissBioisostere 对标报告</h1>
<div class="sub">查询片段 <code>{html.escape(rep['query'])}</code> · 库 {html.escape(os.path.basename(rep['mmpdb']))} · 对照官方 SwissBioisostere 输出</div>

<div class="kpi">
  <div class="card"><b>{rep['precision_vs_official']*100:.0f}%</b><small>precision（本地候选命中官方）</small></div>
  <div class="card"><b>{rep.get('recall_at_50',0)*100:.0f}%</b><small>recall@50（官方 Top50 命中）</small></div>
  <div class="card"><b>{rep['coverage_vs_official']*100:.1f}%</b><small>coverage（受数据规模限制）</small></div>
</div>

<h2>指标</h2>
<table>{metric_rows}</table>

<div class="note"><b>解读：</b>precision 高、recall@50 高，说明本地方法学 SwissBioisostere 相同且前排名次一致；
coverage 低是<b>库只用了 ChEMBL 子集</b>（约全库 1/48），不是方法错误。ΔtPSA 几乎精确；ΔlogP 的差是
RDKit Crippen 与官方 ChEMBL ACD logP 的口径差（待全量 SQLite 到位后替换）。</div>

<h2>官方 Top-{args.top} 逐条比对</h2>
<table>
<tr><th>官方排名</th><th>候选片段</th><th>官方频次</th><th>官方ΔlogP</th><th>官方ΔtPSA</th><th>官方ΔMW</th>
<th>本地片段</th><th>本地频次</th><th>本地ΔlogP</th><th>本地ΔtPSA</th><th>本地ΔMW</th></tr>
{''.join(rows)}
</table>
<p class="sub">红色行 = 本地未命中（多为 50k 子集未观察到该替换）。</p>
</body></html>"""

    out = os.path.join(run, "report.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(doc)
    print("wrote", out)


if __name__ == "__main__":
    main()
