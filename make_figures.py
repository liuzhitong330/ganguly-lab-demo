"""
Two SVG figures for the Ganguly lab demo:
1. circuit_heatmap.svg  — all 12 genes × 7 regions, 3 functional groups highlighted
2. circuit_bars.svg     — focused bar chart showing 3 key circuits:
                          thalamocortical spindle, corticostriatal learning,
                          and GABAergic consolidation
"""
import csv, math
from pathlib import Path

OUT = Path(__file__).parent

STRUCTURES = ["Motor_ctx","PFC","Striatum","Thalamus","Hippocampus","Subst_nigra","Midbrain"]
STRUCT_LABELS = ["Motor\ncortex","PFC","Striatum","Thalamus","Hippocampus","Subst.\nnigra","Midbrain"]

GROUPS = [
    ("sleep spindle / oscillation", ["Cacna1g","Cacna1h","Hcn1","Kcnq2"],  "#2980b9"),
    ("GABAergic consolidation",     ["Gad1","Gad2","Pvalb","Sst"],          "#c0392b"),
    ("corticostriatal learning",    ["Drd1","Drd2","Foxp1","Slc17a7"],      "#27ae60"),
]
GENE_COLOR = {g: c for _, gs, c in GROUPS for g in gs}
ALL_GENES  = [g for _, gs, _ in GROUPS for g in gs]

rows_by_gene = {}
with open(OUT / "expression.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows_by_gene[r["gene"]] = [float(r[s]) for s in STRUCTURES]


def row_norm(vals):
    mn, mx = min(vals), max(vals)
    return [(v - mn) / (mx - mn) if mx > mn else 0.0 for v in vals]


def blue_orange(v):
    """0 → cool white-blue, 1 → deep orange"""
    if v <= 0.5:
        t = v / 0.5
        r = int(215 + (255 - 215) * t)
        g = int(230 + (255 - 230) * t)
        b = 255
    else:
        t = (v - 0.5) / 0.5
        r = 255
        g = int(255 + (100 - 255) * t)
        b = int(255 + (20 - 255) * t)
    return f"rgb({r},{g},{b})"


# ── Figure 1: heatmap ─────────────────────────────────────────────────────────
CELL_W = 72
CELL_H = 26
LEFT   = 95
TOP    = 108
n_genes   = len(ALL_GENES)
n_structs = len(STRUCTURES)
W = LEFT + n_structs * CELL_W + 70
H = TOP + n_genes * CELL_H + 100

THAL_I   = STRUCTURES.index("Thalamus")
MCTX_I   = STRUCTURES.index("Motor_ctx")
STR_I    = STRUCTURES.index("Striatum")

highlights = [
    (THAL_I,  "#f0f8ff", "#1a5c8a",  "spindles"),
    (MCTX_I,  "#fff8f0", "#c05000",  "slow osc."),
    (STR_I,   "#f5fff0", "#1a6b3a",  "learning"),
]

col_bg = ""
for si, bg, _, _ in highlights:
    x = LEFT + si * CELL_W
    col_bg += (f'<rect x="{x}" y="{TOP-56}" width="{CELL_W}" '
               f'height="{n_genes*CELL_H+60}" fill="{bg}" opacity="0.65"/>')

cells = ""
dividers = ""
gi_abs = 0
for _, gs, _ in GROUPS:
    for gene in gs:
        if gene not in rows_by_gene:
            gi_abs += 1; continue
        nv = row_norm(rows_by_gene[gene])
        for si, v in enumerate(nv):
            x = LEFT + si * CELL_W
            y = TOP + gi_abs * CELL_H
            cells += (f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" '
                      f'fill="{blue_orange(v)}" stroke="white" stroke-width="1.5"/>')
            raw = rows_by_gene[gene][si]
            if v > 0.78 and raw > 0.5:
                fc = "white" if v > 0.91 else "#333"
                cells += (f'<text x="{x+CELL_W/2:.0f}" y="{y+CELL_H/2+4:.0f}" '
                          f'text-anchor="middle" font-size="8" fill="{fc}">{raw:.1f}</text>')
        gi_abs += 1
    dividers += (f'<line x1="{LEFT}" y1="{TOP+gi_abs*CELL_H}" '
                 f'x2="{LEFT+n_structs*CELL_W}" y2="{TOP+gi_abs*CELL_H}" '
                 f'stroke="#bbb" stroke-width="1.5"/>')

# Column headers
col_hdrs = ""
for si, lbl in enumerate(STRUCT_LABELS):
    x = LEFT + si * CELL_W + CELL_W // 2
    hi = next(((fc, sub) for idx, _, fc, sub in highlights if idx == si), None)
    fc = hi[0] if hi else "#444"
    fw = "bold" if hi else "normal"
    parts = lbl.split("\n")
    for li, part in enumerate(parts):
        y = TOP - 18 - (len(parts) - 1 - li) * 13
        col_hdrs += (f'<text x="{x}" y="{y}" text-anchor="middle" font-size="10" '
                     f'fill="{fc}" font-weight="{fw}">{part}</text>')
    if hi:
        col_hdrs += (f'<text x="{x}" y="{TOP-4}" text-anchor="middle" '
                     f'font-size="8.5" fill="{hi[0]}" font-style="italic">↑ {hi[1]}</text>')

# Row labels
row_lbls = ""
gi_abs = 0
for _, gs, _ in GROUPS:
    for gene in gs:
        y = TOP + gi_abs * CELL_H + CELL_H // 2 + 4
        row_lbls += (f'<text x="{LEFT-5}" y="{y}" text-anchor="end" font-size="10.5" '
                     f'fill="{GENE_COLOR[gene]}" font-weight="600" '
                     f'font-family="ui-monospace,Menlo,monospace">{gene}</text>')
        gi_abs += 1

# Colorbar
cb_x = LEFT + n_structs * CELL_W + 10
cb_h = n_genes * CELL_H
colorbar = (f'<defs><linearGradient id="cb1" x1="0" y1="1" x2="0" y2="0">'
            f'<stop offset="0%" stop-color="{blue_orange(0)}"/>'
            f'<stop offset="50%" stop-color="{blue_orange(0.5)}"/>'
            f'<stop offset="100%" stop-color="{blue_orange(1)}"/>'
            f'</linearGradient></defs>'
            f'<rect x="{cb_x}" y="{TOP}" width="13" height="{cb_h}" '
            f'fill="url(#cb1)" stroke="#ccc" stroke-width="0.5"/>'
            f'<text x="{cb_x+6}" y="{TOP-4}" text-anchor="middle" font-size="8" fill="#555">high</text>'
            f'<text x="{cb_x+6}" y="{TOP+cb_h+10}" text-anchor="middle" font-size="8" fill="#555">low</text>')

svg1 = f"""<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{W//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Motor Learning and Sleep Consolidation Circuit Genes Across Brain Regions
  </text>
  <text x="{W//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Allen Mouse Brain Atlas · ISH expression energy (row-normalized per gene)
  </text>
  <text x="{W//2}" y="53" text-anchor="middle" font-size="10" fill="#555">
    Spindle/oscillation genes (blue) · GABAergic (red) · Corticostriatal learning (green)
  </text>
  {col_bg}{colorbar}{col_hdrs}{cells}{dividers}{row_lbls}
</svg>"""

with open(OUT / "circuit_heatmap.svg", "w") as f:
    f.write(svg1)
print("Wrote circuit_heatmap.svg")


# ── Figure 2: three-circuit bar chart ────────────────────────────────────────
# Three mini-panels side by side, one per circuit
# Panel A: Thalamocortical spindle — Cacna1g (thalamus) vs Hcn1 (motor ctx)
# Panel B: Corticostriatal learning — Slc17a7 (PFC) vs Drd1/Foxp1 (striatum)
# Panel C: GABA consolidation — Gad1 (striatum) vs Pvalb (motor ctx) vs Sst (PFC)

FW2, FH2 = 700, 400
PAD_T2, PAD_B2 = 90, 70
PANEL_W = (FW2 - 60) // 3
PANEL_H = FH2 - PAD_T2 - PAD_B2

CIRCUITS = [
    {
        "title": "Thalamocortical spindles",
        "subtitle": "Cacna1g (thalamus) drives spindle\nHcn1 (motor ctx) sets rebound",
        "bars": [
            ("Cacna1g", "Thalamus",   "#1a5c8a", "Cav3.1"),
            ("Hcn1",    "Motor_ctx",  "#2980b9", "HCN1"),
            ("Cacna1h", "Thalamus",   "#5dade2", "Cav3.2"),
            ("Kcnq2",   "Motor_ctx",  "#85c1e9", "KCNQ2"),
        ],
        "note": "Thalamus = spindle pacemaker\nMotor ctx = slow-osc recipient",
    },
    {
        "title": "Corticostriatal reversal",
        "subtitle": "PFC→striatum during learning;\nreversal marks consolidation",
        "bars": [
            ("Slc17a7", "PFC",       "#1a6b3a", "VGLUT1"),
            ("Drd1",    "Striatum",  "#27ae60", "DRD1"),
            ("Foxp1",   "Striatum",  "#52be80", "FOXP1"),
            ("Drd2",    "Striatum",  "#7dcea0", "DRD2"),
        ],
        "note": "VGLUT1 in PFC = cortical drive\nDRD1/FOXP1 = striatal recipients",
    },
    {
        "title": "GABAergic consolidation",
        "subtitle": "GABA gates sleep spindle coupling;\nstroke disrupts this balance",
        "bars": [
            ("Gad1",  "Striatum",  "#922b21", "GAD67"),
            ("Gad2",  "Striatum",  "#c0392b", "GAD65"),
            ("Pvalb", "Motor_ctx", "#e74c3c", "PV"),
            ("Sst",   "PFC",       "#f1948a", "SST"),
        ],
        "note": "Striatum = massive GABA output\nPV/SST = cortical IN subtypes",
    },
]

panels = ""
for pi, circ in enumerate(CIRCUITS):
    px = 20 + pi * PANEL_W
    # Panel box
    panels += (f'<rect x="{px}" y="{PAD_T2-5}" width="{PANEL_W-10}" '
               f'height="{PANEL_H+10}" fill="#fafafa" rx="6" stroke="#e0e0e0" stroke-width="1"/>')
    # Title
    panels += (f'<text x="{px+PANEL_W//2-5}" y="{PAD_T2-18}" text-anchor="middle" '
               f'font-size="11" font-weight="600" fill="#222">{circ["title"]}</text>')
    # Subtitle
    for li, line in enumerate(circ["subtitle"].split("\n")):
        panels += (f'<text x="{px+PANEL_W//2-5}" y="{PAD_T2-6+li*12}" text-anchor="middle" '
                   f'font-size="8.5" fill="#777" font-style="italic">{line}</text>')

    # Bars
    n_bars = len(circ["bars"])
    bar_area_w = PANEL_W - 28
    bar_w = bar_area_w / n_bars * 0.65
    bar_gap = bar_area_w / n_bars

    max_val = max(rows_by_gene[gene][STRUCTURES.index(reg)]
                  for gene, reg, col, lbl in circ["bars"])

    for bi, (gene, reg, col, short_lbl) in enumerate(circ["bars"]):
        si = STRUCTURES.index(reg)
        raw = rows_by_gene[gene][si]
        bx = px + 14 + bi * bar_gap + bar_gap * 0.1
        bh = raw / (max_val * 1.15) * PANEL_H
        by = PAD_T2 + PANEL_H - bh

        panels += (f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" '
                   f'fill="{col}" opacity="0.88" rx="3"/>')
        # Value
        panels += (f'<text x="{bx+bar_w/2:.1f}" y="{by-4:.1f}" text-anchor="middle" '
                   f'font-size="9" fill="{col}">{raw:.1f}</text>')
        # Gene label
        panels += (f'<text transform="rotate(-40,{bx+bar_w/2:.1f},{PAD_T2+PANEL_H+14})" '
                   f'x="{bx+bar_w/2:.1f}" y="{PAD_T2+PANEL_H+14}" text-anchor="end" '
                   f'font-size="9.5" fill="{col}" font-family="ui-monospace,Menlo,monospace">'
                   f'{gene}</text>')
        # Region label below gene
        reg_short = reg.replace("Motor_ctx","M1").replace("Subst_nigra","SN")
        panels += (f'<text transform="rotate(-40,{bx+bar_w/2:.1f},{PAD_T2+PANEL_H+28})" '
                   f'x="{bx+bar_w/2:.1f}" y="{PAD_T2+PANEL_H+28}" text-anchor="end" '
                   f'font-size="8" fill="#888">[{reg_short}]</text>')

    # Y-axis
    panels += (f'<line x1="{px+12}" y1="{PAD_T2}" x2="{px+12}" y2="{PAD_T2+PANEL_H}" '
               f'stroke="#ccc" stroke-width="1"/>')
    panels += (f'<line x1="{px+12}" y1="{PAD_T2+PANEL_H}" x2="{px+PANEL_W-12}" '
               f'y2="{PAD_T2+PANEL_H}" stroke="#ccc" stroke-width="1"/>')

    # Note box
    ny = PAD_T2 + PANEL_H + 55
    for li, line in enumerate(circ["note"].split("\n")):
        panels += (f'<text x="{px+PANEL_W//2-5}" y="{ny+li*12}" text-anchor="middle" '
                   f'font-size="8.5" fill="#555">{line}</text>')

svg2 = f"""<svg viewBox="0 0 {FW2} {FH2}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW2//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Three Circuits — Spindle Generation, Corticostriatal Learning, GABAergic Consolidation
  </text>
  <text x="{FW2//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Raw ISH expression energy for key genes in their highest-expressing region (Allen Mouse Brain Atlas)
  </text>
  {panels}
</svg>"""

with open(OUT / "circuit_bars.svg", "w") as f:
    f.write(svg2)
print("Wrote circuit_bars.svg")
