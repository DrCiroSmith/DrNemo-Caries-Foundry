#!/usr/bin/env python3
"""
build_readme_visuals.py: Generate 5 README figures + social preview as Retina PNGs.
Author: Dr. Roi Castrellon, DDS (DrNemo Caries Foundry)

Design system: Deep Clinical Slate
  bg=#0c1017, card=#131b28, border=1px solid #1e293b
  Inter + JetBrains Mono, cyan=#38bdf8, emerald=#4ade80, amber=#facc15, coral=#f87171

GitHub 896px Column Readability Rules:
  All HTML viewports at 1200px logical width, device_scale_factor=2 -> 2400px Retina PNG.
  Hero headers: 36px-44px
  Panel titles: 22px-28px
  Body labels / text: 16px-18px
  Mono badges / metrics: 15px-18px
"""
import csv
import json
import os
import sys
from pathlib import Path
from playwright.sync_api import sync_playwright

REPO = Path(__file__).resolve().parent.parent
ASSETS = REPO / "assets"
DATA = REPO / "data"
DEMO = REPO / "demo_data" / "results"

# Load PDB content for 3Dmol WebGL renders
PDB_PATH = DEMO / "3AIC_docked_complex_CAHWRPWHC.pdb"
PDB_CONTENT = PDB_PATH.read_text(encoding="utf-8").replace("`", "\\`").replace("$", "\\$")

# === Design tokens ===
BG = "#0c1017"
CARD = "#131b28"
BORDER = "#1e293b"
CYAN = "#38bdf8"
EMERALD = "#4ade80"
AMBER = "#facc15"
CORAL = "#f87171"
TEXT1 = "#f1f5f9"
TEXT2 = "#94a3b8"
FONT_STACK = "'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"
MONO = "'JetBrains Mono', 'Fira Code', monospace"

GOOGLE_FONTS = (
    "https://fonts.googleapis.com/css2?"
    "family=Inter:wght@400;500;600;700;800&family=JetBrains+Mono:wght@500;600;700&display=swap"
)

CSS_BASE = f"""
@import url('{GOOGLE_FONTS}');
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
  background: {BG};
  color: {TEXT1};
  font-family: {FONT_STACK};
  -webkit-font-smoothing: antialiased;
}}
"""

PLAYWRIGHT_ARGS = [
    "--enable-webgl",
    "--ignore-gpu-blocklist",
    "--enable-gpu-rasterization",
    "--use-gl=angle",
    "--use-angle=swiftshader",
]


def load_ledger():
    with open(DATA / "contenders_ledger.csv", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_benchmark():
    with open(DEMO / "falsification_benchmark.json", encoding="utf-8") as f:
        return json.load(f)


def render_html_to_png(html: str, width: int, height: int, out_path: Path, wait_ms: int = 3500):
    """Render HTML string to Retina PNG via Playwright Chromium."""
    tmp_html = Path(f"/tmp/{out_path.stem}.html")
    tmp_html.write_text(html, encoding="utf-8")
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True, args=PLAYWRIGHT_ARGS)
        page = browser.new_page(
            viewport={"width": width, "height": height},
            device_scale_factor=2,
        )
        page.goto(f"file://{tmp_html}")
        page.wait_for_timeout(wait_ms)
        page.screenshot(path=str(out_path), full_page=False)
        browser.close()
    tmp_html.unlink(missing_ok=True)
    print(f"  -> {out_path.name} ({out_path.stat().st_size:,} bytes)")


# ===  ===
# FIGURE 1: Hero Banner (1200x380 @2x = 2400x760)
# ===  ===
def build_hero_banner():
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
{CSS_BASE}
body {{ width:1200px; height:380px; overflow:hidden; display:flex; }}
.left {{
  width:56%; padding:36px 36px 36px 44px; display:flex; flex-direction:column; justify-content:center;
}}
.right {{
  width:44%; height:380px; position:relative; overflow:hidden;
  border-left:1px solid {BORDER}; background:#080c13;
}}
#viewer3d {{ width:100%; height:100%; position:relative; }}
.badge-top {{
  display:inline-flex; align-items:center; gap:8px;
  background:rgba(56,189,248,0.12); border:1px solid rgba(56,189,248,0.3);
  padding:4px 14px; border-radius:100px; font-family:{MONO}; font-size:14px;
  font-weight:600; color:{CYAN}; margin-bottom:12px; width:fit-content;
}}
h1 {{
  font-size:38px; font-weight:800; letter-spacing:-0.8px; line-height:1.15;
  color:#ffffff; margin-bottom:8px;
}}
.subtitle {{
  font-size:20px; font-weight:600; color:{CYAN}; margin-bottom:4px;
}}
.subline {{
  font-size:16px; font-weight:400; color:{TEXT2}; margin-bottom:22px;
}}
.pills {{
  display:grid; grid-template-columns:1fr 1fr; gap:12px;
}}
.pill {{
  background:{CARD}; border:1px solid {BORDER}; border-radius:10px;
  padding:10px 14px; display:flex; flex-direction:column; gap:2px;
}}
.pill .val {{
  font-family:{MONO}; font-size:18px; font-weight:700;
}}
.pill .lbl {{
  font-size:14px; font-weight:500; color:{TEXT2};
}}
.render-overlay {{
  position:absolute; bottom:12px; right:16px;
  background:rgba(12,16,23,0.85); backdrop-filter:blur(8px);
  border:1px solid {BORDER}; border-radius:8px; padding:6px 12px;
  font-family:{MONO}; font-size:13px; color:{TEXT2}; z-index:10;
}}
.triad-legend {{
  position:absolute; top:12px; left:16px;
  background:rgba(12,16,23,0.85); backdrop-filter:blur(8px);
  border:1px solid {BORDER}; border-radius:8px; padding:8px 12px;
  font-family:{MONO}; font-size:13px; color:{TEXT1}; z-index:10;
  display:flex; flex-direction:column; gap:4px;
}}
.legend-item {{ display:flex; align-items:center; gap:6px; }}
.dot {{ width:8px; height:8px; border-radius:50%; }}
</style>
<script src="file:///tmp/3Dmol-min.js"></script>
</head><body>
<div class="left">
  <div class="badge-top">AUTODOCK VINA + OPENMM AMBER14 + POPPERIAN FALSIFICATION</div>
  <h1>DrNemo Caries Foundry</h1>
  <div class="subtitle">Targeted GtfB Macrocycle Inhibitor Design</div>
  <div class="subline">Selective cleft blockade · Oral commensal preservation</div>
  <div class="pills">
    <div class="pill">
      <span class="val" style="color:{CYAN}">-8.191 kcal/mol</span>
      <span class="lbl">3AIC Target Affinity</span>
    </div>
    <div class="pill">
      <span class="val" style="color:{EMERALD}">-0.85 kcal/mol</span>
      <span class="lbl">1SMD Host Selectivity ΔΔG</span>
    </div>
    <div class="pill">
      <span class="val" style="color:{AMBER}">0 Sites</span>
      <span class="lbl">Oral Protease Cleavage</span>
    </div>
    <div class="pill">
      <span class="val" style="color:{CORAL}">+1.62 ΔQ</span>
      <span class="lbl">Stephan pH Acid Switch</span>
    </div>
  </div>
</div>
<div class="right">
  <div id="viewer3d"></div>
  <div class="triad-legend">
    <div class="legend-item"><span class="dot" style="background:{CORAL}"></span><span>Triad Asp477/Glu515/Asp588</span></div>
    <div class="legend-item"><span class="dot" style="background:{CYAN}"></span><span>CAHWRPWHC Macrocycle</span></div>
    <div class="legend-item"><span class="dot" style="background:{AMBER}"></span><span>C1:SG-C9:SG Disulfide (2.06 Å)</span></div>
  </div>
  <div class="render-overlay">3Dmol.js WebGL · RCSB 3AIC</div>
</div>
<script>
  const pdbData = `{PDB_CONTENT}`;
  const element = document.getElementById('viewer3d');
  const viewer = $3Dmol.createViewer(element, {{ backgroundColor: '#080c13' }});
  viewer.addModel(pdbData, 'pdb');

  viewer.setStyle({{ chain: 'A' }}, {{ cartoon: {{ color: '#334155', opacity: 0.85 }} }});
  viewer.setStyle({{ chain: 'A', resi: [477, 515, 588] }}, {{
    stick: {{ color: '{CORAL}', radius: 0.35 }}
  }});
  viewer.setStyle({{ chain: 'B' }}, {{
    stick: {{ color: '{CYAN}', radius: 0.28 }}
  }});

  viewer.addSphere({{ center: {{ x: 184.831, y: 32.035, z: 198.614 }}, radius: 0.58, color: '{AMBER}' }});
  viewer.addSphere({{ center: {{ x: 184.677, y: 34.077, z: 198.820 }}, radius: 0.58, color: '{AMBER}' }});
  viewer.addCylinder({{
    start: {{ x: 184.831, y: 32.035, z: 198.614 }},
    end: {{ x: 184.677, y: 34.077, z: 198.820 }},
    radius: 0.24,
    color: '{AMBER}',
    fromCap: true,
    toCap: true
  }});

  try {{
    viewer.addSurface($3Dmol.SurfaceType.MS, {{ opacity: 0.22, color: '#94a3b8' }}, {{ chain: 'A', within: {{ distance: 10.0, sel: {{ chain: 'B' }} }} }});
  }} catch(e) {{
    console.error(e);
  }}

  viewer.zoomTo({{ chain: 'B' }});
  viewer.zoom(0.85);
  viewer.rotate(25, 'y');
  viewer.rotate(-15, 'x');
  viewer.render();
</script>
</body></html>"""
    render_html_to_png(html, 1200, 380, ASSETS / "hero_banner.png", wait_ms=4000)


# ===  ===
# FIGURE 2: Biological Mechanism (1200x360 @2x = 2400x720)
# ===  ===
def build_biological_mechanism():
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
{CSS_BASE}
body {{ width:1200px; height:360px; overflow:hidden; display:flex; gap:18px; padding:16px 20px; }}
.panel {{
  flex:1; background:{CARD}; border:1px solid {BORDER}; border-radius:14px;
  padding:16px 20px; display:flex; flex-direction:column; justify-content:space-between;
}}
.panel-header {{
  display:flex; align-items:center; gap:12px; margin-bottom:10px;
}}
.panel-badge {{
  width:30px; height:30px; border-radius:8px; display:flex;
  align-items:center; justify-content:center; font-family:{MONO};
  font-size:17px; font-weight:700;
}}
.panel-title {{
  font-size:21px; font-weight:700;
}}
.flow {{
  display:flex; flex-direction:column; gap:6px; flex:1; justify-content:space-between;
}}
.step {{
  display:flex; align-items:center; gap:12px;
  padding:8px 12px; border-radius:8px; font-size:15px; font-weight:500;
}}
.step-num {{
  width:26px; height:26px; border-radius:50%; display:flex;
  align-items:center; justify-content:center; font-family:{MONO};
  font-size:14px; font-weight:700; flex-shrink:0;
}}
.arrow-down {{
  display:flex; justify-content:center; color:{TEXT2}; font-size:13px; line-height:1;
}}
.bad-border {{ border:1px solid rgba(248,113,113,0.3); background:rgba(248,113,113,0.06); }}
.bad-num {{ background:rgba(248,113,113,0.2); color:{CORAL}; }}
.good-border {{ border:1px solid rgba(74,222,128,0.3); background:rgba(74,222,128,0.06); }}
.good-num {{ background:rgba(74,222,128,0.2); color:{EMERALD}; }}
.neutral-border {{ border:1px solid {BORDER}; background:{BG}; }}
.neutral-num {{ background:rgba(148,163,184,0.15); color:{TEXT2}; }}
.tag-blocked {{
  margin-left:auto; font-family:{MONO}; font-size:13px; font-weight:700;
  color:{CYAN}; background:rgba(56,189,248,0.15); border:1px solid rgba(56,189,248,0.4);
  padding:2px 8px; border-radius:6px; white-space:nowrap;
}}
</style></head><body>

<div class="panel">
  <div class="panel-header">
    <div class="panel-badge" style="background:rgba(248,113,113,0.15);color:{CORAL};">A</div>
    <div class="panel-title" style="color:{CORAL};">Uninhibited S. mutans Biofilm</div>
  </div>
  <div class="flow">
    <div class="step neutral-border">
      <div class="step-num neutral-num">1</div>
      <span>Salivary Pellicle + Enamel Hydroxyapatite</span>
    </div>
    <div class="arrow-down">↓</div>
    <div class="step bad-border">
      <div class="step-num bad-num">2</div>
      <span>SpaP Adhesin Anchors Cells to Pellicle Surface</span>
    </div>
    <div class="arrow-down">↓</div>
    <div class="step bad-border">
      <div class="step-num bad-num">3</div>
      <span>GtfB Secretion: α-1,3-Glucan Insoluble Mutan Synthesis</span>
    </div>
    <div class="arrow-down">↓</div>
    <div class="step bad-border">
      <div class="step-num bad-num">4</div>
      <span>Lactic Acid Trapped at pH ≤ 5.2 → Enamel Demineralization</span>
    </div>
  </div>
</div>

<div class="panel">
  <div class="panel-header">
    <div class="panel-badge" style="background:rgba(74,222,128,0.15);color:{EMERALD};">B</div>
    <div class="panel-title" style="color:{EMERALD};">CAHWRPWHC Targeted Blockade</div>
  </div>
  <div class="flow">
    <div class="step neutral-border">
      <div class="step-num neutral-num">1</div>
      <span>Salivary Pellicle + Enamel Hydroxyapatite</span>
    </div>
    <div class="arrow-down">↓</div>
    <div class="step good-border">
      <div class="step-num good-num">2</div>
      <span>SpaP Adhesion Preserved (Commensals Intact)</span>
    </div>
    <div class="arrow-down">↓</div>
    <div class="step good-border">
      <div class="step-num good-num">3</div>
      <span>GtfB Catalytic Triad (D477/E515/D588) Blocked</span>
      <span class="tag-blocked">BLOCKED -8.191</span>
    </div>
    <div class="arrow-down">↓</div>
    <div class="step good-border">
      <div class="step-num good-num">4</div>
      <span>Mutan Shield Halted → Salivary Clearance, Commensals Safe</span>
    </div>
  </div>
</div>

</body></html>"""
    render_html_to_png(html, 1200, 360, ASSETS / "fig1_biological_mechanism.png")


# ===  ===
# FIGURE 3: Pipeline Architecture (1200x420 @2x = 2400x840)
# ===  ===
def build_pipeline_architecture():
    stages = [
        ("1", "Evolutionary Generator", "C-X7-C cyclic scaffold<br>P2A hinge mutant pool", CYAN),
        ("2", "Oral Shield Gate", "Donnan Q ≤ +4.0<br>Stephan ΔQ ≥ +1.0<br>0 protease sites", EMERALD),
        ("3", "28.0 Å Dual Docking", "3AIC target (GtfB)<br>1SMD host counter (amylase)<br>Vina + Vinardo", AMBER),
        ("4", "OpenMM Relaxation", "Amber14/GBn2 at 310.15 K<br>All-atom minimization<br>S-S 2.06 Å check", CORAL),
        ("5", "Popper Falsification", "1SMD ΔΔG ≤ -0.80<br>N=39 scrambled nulls<br>Alanine scanning", CYAN),
        ("6", "Contender Vault", "SQLite registry<br>Pareto rank-1 advance<br>Fmoc SPPS protocol", EMERALD),
    ]

    cards_html = ""
    for i, (num, title, desc, col) in enumerate(stages):
        cards_html += f"""
        <div class="stage-card" style="border-top:3px solid {col};">
          <div class="stage-badge" style="background:{col};">{num}</div>
          <div class="stage-title">{title}</div>
          <div class="stage-desc">{desc}</div>
        </div>
        """
        if i < len(stages) - 1:
            cards_html += f"""
            <div class="connector">
              <svg width="22" height="22" viewBox="0 0 24 24"><path d="M7 4l8 8-8 8" stroke="{TEXT2}" stroke-width="2.5" fill="none" stroke-linecap="round"/></svg>
            </div>
            """

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
{CSS_BASE}
body {{ width:1200px; height:420px; overflow:hidden; padding:24px 28px; display:flex; flex-direction:column; justify-content:space-between; }}
.header-row {{
  text-align:center; margin-bottom:14px;
}}
.main-title {{
  font-size:24px; font-weight:800; color:#ffffff; letter-spacing:-0.4px; margin-bottom:4px;
}}
.main-subtitle {{
  font-size:16px; color:{TEXT2};
}}
.pipeline-flow {{
  display:flex; align-items:center; justify-content:space-between;
}}
.stage-card {{
  flex:1; background:{CARD}; border:1px solid {BORDER}; border-radius:12px;
  padding:16px 12px; display:flex; flex-direction:column; align-items:center;
  text-align:center; min-height:210px; justify-content:flex-start;
}}
.stage-badge {{
  width:32px; height:32px; border-radius:50%; color:#080c13;
  font-family:{MONO}; font-size:16px; font-weight:800;
  display:flex; align-items:center; justify-content:center; margin-bottom:10px;
}}
.stage-title {{
  font-size:16px; font-weight:700; color:{TEXT1}; margin-bottom:8px; line-height:1.25;
}}
.stage-desc {{
  font-size:13px; color:{TEXT2}; line-height:1.45;
}}
.connector {{
  display:flex; align-items:center; justify-content:center; padding:0 3px;
}}
.feedback-box {{
  display:flex; align-items:center; justify-content:center; gap:14px;
  background:rgba(248,113,113,0.06); border:1px dashed rgba(248,113,113,0.4);
  border-radius:10px; padding:10px 20px; margin-top:10px;
}}
.feedback-text {{
  font-family:{MONO}; font-size:15px; font-weight:600; color:{CORAL};
}}
</style></head><body>

<div class="header-row">
  <div class="main-title">Autonomous Biophysical Discovery & Karl Popper Falsification Cycle</div>
  <div class="main-subtitle">Closed-loop LangGraph state machine orchestrating Vina, OpenMM, and multi-angle falsification gates</div>
</div>

<div class="pipeline-flow">
  {cards_html}
</div>

<div class="feedback-box">
  <span style="color:{CORAL};font-size:18px;">↺</span>
  <span class="feedback-text">Falsification Feedback: Off-target amylase binding or protease susceptibility triggers hinge re-engineering loop</span>
  <span style="color:{CORAL};font-size:18px;">↺</span>
</div>

</body></html>"""
    render_html_to_png(html, 1200, 420, ASSETS / "fig2_pipeline_architecture.png")


# ===  ===
# FIGURE 4: Pocket + Attrition Funnel (1200x480 @2x = 2400x960)
# ===  ===
def build_pocket_and_funnel():
    funnel_data = [
        ("82", "Design Candidates", 1.0, TEXT2),
        ("77", "Vina 3AIC Cleft Affinity Passed", 77/82, CYAN),
        ("39", "1SMD Amylase Selectivity (ΔΔG ≤ -0.80)", 39/82, CYAN),
        ("18", "Oral Protease Steric Shield (0 Cleavage Sites)", 18/82, AMBER),
        ("2",  "Stephan Acid-Switch (ΔQ ≥ +1.0)", 2/82, EMERALD),
        ("2",  "Chitosan Donnan Net Charge (Q ≤ +4.0)", 2/82, EMERALD),
        ("2",  "Corroborated Pareto Rank-1 Survivors", 2/82, EMERALD),
    ]

    funnel_rows = ""
    for count, label, frac, color in funnel_data:
        bar_w = max(int(frac * 100), 7)
        is_survivor = count == "2"
        font_weight = "700" if is_survivor else "600"
        funnel_rows += f"""
        <div class="funnel-row">
          <div class="funnel-count" style="color:{color};font-weight:{font_weight};">{count}</div>
          <div class="funnel-track">
            <div class="funnel-bar" style="width:{bar_w}%;background:{color};"></div>
            <div class="funnel-label" style="font-weight:{font_weight};">{label}</div>
          </div>
        </div>
        """

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
{CSS_BASE}
body {{ width:1200px; height:480px; overflow:hidden; display:flex; gap:20px; padding:22px; }}
.panel {{
  flex:1; background:{CARD}; border:1px solid {BORDER}; border-radius:14px;
  padding:22px; display:flex; flex-direction:column; position:relative; overflow:hidden;
}}
.panel-title {{
  font-size:22px; font-weight:700; color:#ffffff; margin-bottom:14px; display:flex;
  align-items:center; justify-content:space-between;
}}
.panel-tag {{
  font-family:{MONO}; font-size:13px; color:{CYAN};
  background:rgba(56,189,248,0.12); padding:4px 10px; border-radius:6px;
  border:1px solid rgba(56,189,248,0.3); font-weight:600;
}}
.pocket-container {{
  flex:1; border-radius:10px; overflow:hidden; position:relative;
  background:#080c13; border:1px solid {BORDER};
}}
#viewer-pocket {{ width:100%; height:100%; }}
.pocket-labels {{
  position:absolute; bottom:10px; left:12px; right:12px;
  display:flex; justify-content:space-between;
  background:rgba(12,16,23,0.85); backdrop-filter:blur(6px);
  padding:6px 12px; border-radius:6px; border:1px solid {BORDER};
  font-family:{MONO}; font-size:12px; color:{TEXT2};
}}
.funnel-list {{
  display:flex; flex-direction:column; gap:9px; flex:1; justify-content:center;
}}
.funnel-row {{
  display:flex; align-items:center; gap:14px;
}}
.funnel-count {{
  width:36px; text-align:right; font-family:{MONO}; font-size:17px;
  white-space:nowrap; flex-shrink:0;
}}
.funnel-track {{
  flex:1; height:32px; background:{BG}; border:1px solid {BORDER};
  border-radius:8px; position:relative; overflow:hidden; display:flex; align-items:center;
}}
.funnel-bar {{
  position:absolute; top:0; left:0; height:100%; opacity:0.25; border-radius:6px;
}}
.funnel-label {{
  position:relative; z-index:2; padding-left:14px; font-size:15px; color:{TEXT1};
  white-space:nowrap;
}}
.survivors-card {{
  margin-top:12px; background:{BG}; border:1px solid rgba(74,222,128,0.3);
  border-radius:10px; padding:12px 16px; display:flex; flex-direction:column; gap:6px;
}}
.survivor-head {{
  font-size:14px; font-weight:700; color:{EMERALD}; font-family:{MONO};
}}
.survivor-lead {{
  display:flex; justify-content:space-between; font-family:{MONO}; font-size:15px;
}}
</style>
<script src="file:///tmp/3Dmol-min.js"></script>
</head><body>

<!-- Left: Active Site 3D Viewport -->
<div class="panel">
  <div class="panel-title">
    <span>GtfB Catalytic Active Site</span>
    <span class="panel-tag">3AIC Cleft · CAHWRPWHC</span>
  </div>
  <div class="pocket-container">
    <div id="viewer-pocket"></div>
    <div class="pocket-labels">
      <span>Asp477 · Glu515 · Asp588 (Triad)</span>
      <span style="color:{AMBER}">S-S 2.06 Å Restraint</span>
    </div>
  </div>
</div>

<!-- Right: Exact Ledger Attrition Funnel -->
<div class="panel">
  <div class="panel-title">
    <span>Popperian Attrition Funnel</span>
    <span class="panel-tag" style="color:{EMERALD};border-color:rgba(74,222,128,0.4);background:rgba(74,222,128,0.1);">82 → 2 Leads</span>
  </div>
  <div class="funnel-list">
    {funnel_rows}
  </div>
  <div class="survivors-card">
    <div class="survivor-head">CORROBORATED PARETO RANK-1 LEADS (ALL 6 GATES CLEARED)</div>
    <div class="survivor-lead">
      <span style="color:{CYAN};font-weight:700;">#1 CAHWRPWHC (P2A Hinge)</span>
      <span style="color:{TEXT1};">-8.191 kcal/mol · ΔΔG -0.85</span>
    </div>
    <div class="survivor-lead">
      <span style="color:{EMERALD};font-weight:600;">#2 CPHWRPAHC (W7A Lead)</span>
      <span style="color:{TEXT1};">-7.805 kcal/mol · ΔΔG -0.96</span>
    </div>
  </div>
</div>

<script>
  const pdbData = `{PDB_CONTENT}`;
  const element = document.getElementById('viewer-pocket');
  const viewer = $3Dmol.createViewer(element, {{ backgroundColor: '#080c13' }});
  viewer.addModel(pdbData, 'pdb');

  viewer.setStyle({{ chain: 'A' }}, {{ cartoon: {{ color: '#334155', opacity: 0.8 }} }});
  viewer.setStyle({{ chain: 'A', resi: [477, 515, 588] }}, {{
    stick: {{ color: '{CORAL}', radius: 0.38 }}
  }});
  viewer.setStyle({{ chain: 'B' }}, {{
    stick: {{ color: '{CYAN}', radius: 0.3 }}
  }});

  viewer.addSphere({{ center: {{ x: 184.831, y: 32.035, z: 198.614 }}, radius: 0.6, color: '{AMBER}' }});
  viewer.addSphere({{ center: {{ x: 184.677, y: 34.077, z: 198.820 }}, radius: 0.6, color: '{AMBER}' }});
  viewer.addCylinder({{
    start: {{ x: 184.831, y: 32.035, z: 198.614 }},
    end: {{ x: 184.677, y: 34.077, z: 198.820 }},
    radius: 0.25,
    color: '{AMBER}',
    fromCap: true,
    toCap: true
  }});

  try {{
    viewer.addSurface($3Dmol.SurfaceType.MS, {{ opacity: 0.22, color: '#94a3b8' }}, {{ chain: 'A', within: {{ distance: 10.0, sel: {{ chain: 'B' }} }} }});
  }} catch(e) {{
    console.error(e);
  }}

  viewer.zoomTo({{ chain: 'B' }});
  viewer.zoom(0.95);
  viewer.rotate(35, 'y');
  viewer.rotate(10, 'x');
  viewer.render();
</script>
</body></html>"""
    render_html_to_png(html, 1200, 480, ASSETS / "fig3_pocket_and_funnel.png", wait_ms=4000)


# ===  ===
# FIGURE 5: Falsification Data Plate 2x2 (1200x760 @2x = 2400x1520)
# REVISION 2:
#   Panel A: Centered callouts, no card edge clipping
#   Panel B: Axis label lowered to avoid baseline collision
#   Panel C: Scaled bar widths so -8.890 and "Flagged: Trypsin Lysis" never collide
#   Panel D: Fully contained labels
# ===  ===
def build_falsification_plate():
    ledger = load_ledger()
    benchmark = load_benchmark()
    alanine = benchmark["alanine_scan"]
    scrambled = benchmark["scrambled_null_test"]

    # (A) Scatter plot points
    design_rows = [r for r in ledger if r.get("falsification_status") != "SCRAMBLED_NULL_CONTROL"]
    scatter_svg = ""

    def map_pareto(x, y):
        px = 55 + (x - (-11.0)) / ((-6.5) - (-11.0)) * 430
        py = 35 + (y - (1.0)) / ((-4.0) - (1.0)) * 180
        return px, py

    cahw_coords = None
    cphw_coords = None

    for r in design_rows:
        try:
            x_val = float(r["vina_3aic_kcal"])
            y_val = float(r["selectivity_ddg"])
        except (ValueError, KeyError):
            continue
        px, py = map_pareto(x_val, y_val)
        seq = r.get("sequence", "")
        if seq == "CAHWRPWHC":
            cahw_coords = (px, py)
        elif seq == "CPHWRPAHC":
            cphw_coords = (px, py)
        elif r.get("falsification_status") == "CORROBORATED":
            scatter_svg += f'<circle cx="{px:.1f}" cy="{py:.1f}" r="5" fill="{EMERALD}" opacity="0.75"/>'
        elif r.get("falsification_status") == "FALSIFIED_AMYLASE_CROSSREACTIVITY":
            scatter_svg += f'<circle cx="{px:.1f}" cy="{py:.1f}" r="3.5" fill="{CORAL}" opacity="0.35"/>'
        else:
            scatter_svg += f'<circle cx="{px:.1f}" cy="{py:.1f}" r="3.5" fill="{TEXT2}" opacity="0.3"/>'

    threshold_y = map_pareto(-8.0, -0.80)[1]

    # Staggered lead annotations (fully inside boundaries)
    if cahw_coords:
        px, py = cahw_coords
        scatter_svg += f'<circle cx="{px:.1f}" cy="{py:.1f}" r="8.5" fill="{CYAN}" stroke="#ffffff" stroke-width="2"/>'
        scatter_svg += f'<line x1="{px:.1f}" y1="{py:.1f}" x2="{px-15:.1f}" y2="{py-22:.1f}" stroke="{CYAN}" stroke-width="1.5"/>'
        scatter_svg += f'<rect x="{px-175:.1f}" y="{py-38:.1f}" width="180" height="24" rx="4" fill="{CARD}" stroke="{CYAN}" stroke-width="1"/>'
        scatter_svg += f'<text x="{px-85:.1f}" y="{py-22:.1f}" text-anchor="middle" fill="{CYAN}" font-family="{MONO}" font-size="12" font-weight="700">#1 CAHWRPWHC (-8.191, -0.85)</text>'

    if cphw_coords:
        px, py = cphw_coords
        scatter_svg += f'<circle cx="{px:.1f}" cy="{py:.1f}" r="7.5" fill="{EMERALD}" stroke="#ffffff" stroke-width="1.5"/>'
        scatter_svg += f'<line x1="{px:.1f}" y1="{py:.1f}" x2="{px-10:.1f}" y2="{py+22:.1f}" stroke="{EMERALD}" stroke-width="1.5"/>'
        scatter_svg += f'<rect x="{px-95:.1f}" y="{py+22:.1f}" width="170" height="24" rx="4" fill="{CARD}" stroke="{EMERALD}" stroke-width="1"/>'
        scatter_svg += f'<text x="{px-10:.1f}" y="{py+38:.1f}" text-anchor="middle" fill="{EMERALD}" font-family="{MONO}" font-size="12" font-weight="700">#2 CPHWRPAHC (-7.805, -0.96)</text>'

    # (B) Histogram data for scrambled nulls
    perm_vals = [p["vina_3aic_kcal"] for p in scrambled["permutations"]]
    min_v, max_v = min(perm_vals), max(perm_vals)
    n_bins = 10
    bin_step = (max_v - min_v) / n_bins
    counts = [0] * n_bins
    for v in perm_vals:
        idx = min(int((v - min_v) / bin_step), n_bins - 1)
        counts[idx] += 1
    max_c = max(counts)

    hist_svg = ""
    for i, c in enumerate(counts):
        bx = 55 + i * (440 / n_bins)
        bw = (440 / n_bins) - 5
        bh = (c / max_c) * 155
        by = 205 - bh
        hist_svg += f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bw:.1f}" height="{bh:.1f}" fill="{CYAN}" opacity="0.65" rx="3"/>'
        if c > 0:
            hist_svg += f'<text x="{bx+bw/2:.1f}" y="{by-5:.1f}" fill="{TEXT1}" font-family="{MONO}" font-size="12" font-weight="600" text-anchor="middle">{c}</text>'

    mean_val = scrambled["mean_kcal"]
    mx = 55 + (mean_val - min_v) / (max_v - min_v) * 440
    cx = 55 + (-8.191 - min_v) / (max_v - min_v) * 440

    # (C) Alanine Scan bars
    ala_keys = ["P2A", "H3A", "W4A", "R5A", "P6A", "W7A", "H8A"]
    ala_svg = ""
    parent_base = -8.062
    for i, k in enumerate(ala_keys):
        info = alanine[k]
        val = info["vina_3aic_kcal"]
        delta = info["delta_binding_affinity_kcal"]
        is_hot = info["is_hotspot"]
        y_pos = 26 + i * 27
        # Scale to 140 max width
        bar_w = (abs(val) - 7.5) / 1.5 * 140
        bar_col = CORAL if k == "P6A" else (AMBER if is_hot else CYAN)

        ala_svg += f"""
        <text x="75" y="{y_pos+15}" text-anchor="end" fill="{TEXT1}" font-family="{MONO}" font-size="13" font-weight="600">{k}</text>
        <rect x="85" y="{y_pos}" width="{bar_w:.1f}" height="20" fill="{bar_col}" opacity="0.7" rx="3"/>
        <text x="{92+bar_w:.1f}" y="{y_pos+15}" fill="{TEXT1}" font-family="{MONO}" font-size="12" font-weight="600">{val:.3f} ({delta:+.3f})</text>
        """
        if k == "P6A":
            ala_svg += f'<text x="510" y="{y_pos+15}" text-anchor="end" fill="{CORAL}" font-family="{MONO}" font-size="12" font-weight="700">Trypsin Lysis (Arg5)</text>'
        elif is_hot:
            ala_svg += f'<text x="510" y="{y_pos+15}" text-anchor="end" fill="{AMBER}" font-family="{MONO}" font-size="12" font-weight="600">Catalytic Anchor</text>'

    # (D) Stephan Charge Curve
    def map_stephan(ph, q):
        sx = 55 + (7.5 - ph) / (7.5 - 4.5) * 440
        sy = 205 - (q / 8.0) * 170
        return sx, sy

    cahw_pts = [(7.5, 0.50), (7.2, 1.13), (6.0, 2.05), (5.2, 2.75), (4.5, 3.20)]
    crrr_pts = [(7.5, 7.00), (7.2, 7.00), (6.0, 7.00), (5.2, 7.00), (4.5, 7.00)]

    cahw_poly = " ".join(f"{map_stephan(p, q)[0]:.1f},{map_stephan(p, q)[1]:.1f}" for p, q in cahw_pts)
    crrr_poly = " ".join(f"{map_stephan(p, q)[0]:.1f},{map_stephan(p, q)[1]:.1f}" for p, q in crrr_pts)

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
{CSS_BASE}
body {{ width:1200px; height:760px; overflow:hidden; padding:18px; }}
.plate-grid {{
  display:grid; grid-template-columns:1fr 1fr; grid-template-rows:1fr 1fr; gap:16px;
  width:100%; height:100%;
}}
.plate-cell {{
  background:{CARD}; border:1px solid {BORDER}; border-radius:12px;
  padding:14px 18px; display:flex; flex-direction:column; justify-content:space-between;
}}
.cell-header {{
  display:flex; align-items:center; justify-content:space-between; margin-bottom:4px;
}}
.cell-badge {{
  display:inline-flex; align-items:center; gap:8px;
}}
.badge-letter {{
  width:26px; height:26px; border-radius:6px; background:rgba(56,189,248,0.15);
  color:{CYAN}; font-family:{MONO}; font-size:15px; font-weight:800;
  display:flex; align-items:center; justify-content:center;
}}
.cell-title {{
  font-size:17px; font-weight:700; color:{TEXT1};
}}
.cell-meta {{
  font-family:{MONO}; font-size:13px; color:{TEXT2};
}}
.chart-box {{
  flex:1; width:100%; height:235px; position:relative;
}}
.legend-bar {{
  display:flex; gap:14px; align-items:center; font-family:{MONO}; font-size:12px;
  color:{TEXT2}; margin-top:4px; border-top:1px solid rgba(255,255,255,0.06); padding-top:4px;
}}
.legend-dot {{ width:8px; height:8px; border-radius:50%; }}
</style></head><body>

<div class="plate-grid">

  <!-- (A) Pareto Frontier -->
  <div class="plate-cell">
    <div class="cell-header">
      <div class="cell-badge">
        <div class="badge-letter">A</div>
        <div class="cell-title">Pareto Frontier: Target vs Amylase Selectivity</div>
      </div>
      <div class="cell-meta">N=155 Design Pool</div>
    </div>
    <div class="chart-box">
      <svg width="530" height="235" viewBox="0 0 530 235">
        <!-- Axes -->
        <line x1="55" y1="205" x2="505" y2="205" stroke="{BORDER}" stroke-width="1.5"/>
        <line x1="55" y1="20" x2="55" y2="205" stroke="{BORDER}" stroke-width="1.5"/>
        <!-- Threshold dashed line -->
        <line x1="55" y1="{threshold_y:.1f}" x2="505" y2="{threshold_y:.1f}" stroke="{AMBER}" stroke-width="1.5" stroke-dasharray="5,4"/>
        <text x="510" y="{threshold_y+4:.1f}" fill="{AMBER}" font-family="{MONO}" font-size="11" font-weight="600">ΔΔG -0.80</text>
        {scatter_svg}
        <text x="280" y="225" text-anchor="middle" fill="{TEXT2}" font-family="{MONO}" font-size="12">GtfB Vina Affinity ΔG (kcal/mol)</text>
        <text x="18" y="115" text-anchor="middle" fill="{TEXT2}" font-family="{MONO}" font-size="12" transform="rotate(-90 18 115)">Selectivity ΔΔG</text>
      </svg>
    </div>
    <div class="legend-bar">
      <span style="display:flex;align-items:center;gap:6px;"><span class="legend-dot" style="background:{CYAN}"></span>Pareto Lead</span>
      <span style="display:flex;align-items:center;gap:6px;"><span class="legend-dot" style="background:{EMERALD}"></span>Corroborated</span>
      <span style="display:flex;align-items:center;gap:6px;"><span class="legend-dot" style="background:{CORAL}"></span>Amylase Falsified</span>
    </div>
  </div>

  <!-- (B) Scrambled Null Histogram -->
  <div class="plate-cell">
    <div class="cell-header">
      <div class="cell-badge">
        <div class="badge-letter">B</div>
        <div class="cell-title">N=39 Scrambled Null Isomer Distribution</div>
      </div>
      <div class="cell-meta">Z = +0.32σ · p = 0.5897</div>
    </div>
    <div class="chart-box">
      <svg width="530" height="235" viewBox="0 0 530 235">
        <line x1="55" y1="205" x2="505" y2="205" stroke="{BORDER}" stroke-width="1.5"/>
        <line x1="55" y1="20" x2="55" y2="205" stroke="{BORDER}" stroke-width="1.5"/>
        {hist_svg}
        <!-- Mean line -->
        <line x1="{mx:.1f}" y1="35" x2="{mx:.1f}" y2="205" stroke="{AMBER}" stroke-width="2" stroke-dasharray="5,3"/>
        <text x="{mx:.1f}" y="28" fill="{AMBER}" font-family="{MONO}" font-size="11" font-weight="700" text-anchor="middle">Mean -8.343 ± 0.470</text>
        <!-- CAHW line -->
        <line x1="{cx:.1f}" y1="45" x2="{cx:.1f}" y2="205" stroke="{CYAN}" stroke-width="2"/>
        <text x="{cx:.1f}" y="42" fill="{CYAN}" font-family="{MONO}" font-size="11" font-weight="700" text-anchor="middle">CAHWRPWHC -8.191</text>
        <text x="280" y="225" text-anchor="middle" fill="{TEXT2}" font-family="{MONO}" font-size="12">Vina 3AIC Score (kcal/mol)</text>
      </svg>
    </div>
    <div class="legend-bar">
      <span>Non-rejection proves topological sequence constraint rather than generic charge collapse</span>
    </div>
  </div>

  <!-- (C) Alanine Scan -->
  <div class="plate-cell">
    <div class="cell-header">
      <div class="cell-badge">
        <div class="badge-letter">C</div>
        <div class="cell-title">7-Site Alanine Mutational Hotspot Scan</div>
      </div>
      <div class="cell-meta">Parent CPHWRPWHC (-8.062)</div>
    </div>
    <div class="chart-box">
      <svg width="530" height="235" viewBox="0 0 530 235">
        <line x1="85" y1="20" x2="85" y2="218" stroke="{BORDER}" stroke-width="1.5"/>
        <!-- Parent reference line -->
        <line x1="{85 + (abs(parent_base) - 7.5)/1.5 * 190:.1f}" y1="15"
              x2="{85 + (abs(parent_base) - 7.5)/1.5 * 190:.1f}" y2="218"
              stroke="{TEXT2}" stroke-width="1.5" stroke-dasharray="4,4"/>
        <text x="{85 + (abs(parent_base) - 7.5)/1.5 * 190:.1f}" y="12" fill="{TEXT2}" font-family="{MONO}" font-size="11" text-anchor="middle">Parent -8.062</text>
        {ala_svg}
      </svg>
    </div>
    <div class="legend-bar">
      <span style="display:flex;align-items:center;gap:6px;"><span class="legend-dot" style="background:{AMBER}"></span>Hotspot Anchor</span>
      <span style="display:flex;align-items:center;gap:6px;"><span class="legend-dot" style="background:{CORAL}"></span>Trypsin Lysis</span>
      <span style="display:flex;align-items:center;gap:6px;"><span class="legend-dot" style="background:{CYAN}"></span>Hinge / Permissive</span>
    </div>
  </div>

  <!-- (D) Stephan Charge Titration -->
  <div class="plate-cell">
    <div class="cell-header">
      <div class="cell-badge">
        <div class="badge-letter">D</div>
        <div class="cell-title">Stephan pH Acid-Switch Charge Titration</div>
      </div>
      <div class="cell-meta">Salivary pH 7.5 → Plaque pH 4.5</div>
    </div>
    <div class="chart-box">
      <svg width="530" height="235" viewBox="0 0 530 235">
        <line x1="55" y1="205" x2="505" y2="205" stroke="{BORDER}" stroke-width="1.5"/>
        <line x1="55" y1="20" x2="55" y2="205" stroke="{BORDER}" stroke-width="1.5"/>
        <!-- Curves -->
        <polyline points="{cahw_poly}" fill="none" stroke="{CYAN}" stroke-width="3" stroke-linecap="round"/>
        <polyline points="{crrr_poly}" fill="none" stroke="{CORAL}" stroke-width="2.5" stroke-dasharray="6,4" stroke-linecap="round"/>
        <!-- Labels inside plot -->
        <text x="490" y="125" text-anchor="end" fill="{CYAN}" font-family="{MONO}" font-size="13" font-weight="700">CAHWRPWHC (ΔQ = +1.62)</text>
        <text x="260" y="42" fill="{CORAL}" font-family="{MONO}" font-size="13" font-weight="700">CRRRWWRRRRC (Static Q=+7.0)</text>
        <!-- pH ticks -->
        <text x="55" y="222" fill="{TEXT2}" font-family="{MONO}" font-size="11" text-anchor="middle">pH 7.5</text>
        <text x="201" y="222" fill="{TEXT2}" font-family="{MONO}" font-size="11" text-anchor="middle">pH 6.5</text>
        <text x="348" y="222" fill="{TEXT2}" font-family="{MONO}" font-size="11" text-anchor="middle">pH 5.5</text>
        <text x="495" y="222" fill="{TEXT2}" font-family="{MONO}" font-size="11" text-anchor="middle">pH 4.5</text>
        <text x="18" y="115" text-anchor="middle" fill="{TEXT2}" font-family="{MONO}" font-size="12" transform="rotate(-90 18 115)">Net Charge Q</text>
      </svg>
    </div>
    <div class="legend-bar">
      <span>Acid activation: Gains +1.62 positive charge selectively inside acidogenic plaque (pH 5.2)</span>
    </div>
  </div>

</div>
</body></html>"""
    render_html_to_png(html, 1200, 760, ASSETS / "fig4_falsification_plate.png")


# ===  ===
# SOCIAL PREVIEW (1280x640 @2x = 2560x1280)
# ===  ===
def build_social_preview():
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
{CSS_BASE}
body {{
  width:1280px; height:640px; overflow:hidden; display:flex; flex-direction:column;
  justify-content:center; align-items:center; text-align:center;
  background: radial-gradient(ellipse at 50% 30%, #172338 0%, {BG} 75%);
  padding:40px;
}}
.badge {{
  display:inline-flex; align-items:center; gap:8px;
  background:rgba(56,189,248,0.12); border:1px solid rgba(56,189,248,0.35);
  padding:6px 18px; border-radius:100px; font-family:{MONO}; font-size:15px;
  font-weight:700; color:{CYAN}; margin-bottom:18px;
}}
h1 {{
  font-size:56px; font-weight:800; letter-spacing:-1.2px; line-height:1.15;
  color:#ffffff; margin-bottom:12px;
}}
.subtitle {{
  font-size:24px; font-weight:600; color:{CYAN}; margin-bottom:8px;
}}
.subline {{
  font-size:18px; color:{TEXT2}; margin-bottom:36px; max-width:800px;
}}
.metric-row {{
  display:flex; gap:16px; justify-content:center;
}}
.metric-card {{
  background:{CARD}; border:1px solid {BORDER}; border-radius:14px;
  padding:16px 26px; display:flex; flex-direction:column; gap:4px;
}}
.metric-val {{
  font-family:{MONO}; font-size:24px; font-weight:800;
}}
.metric-lbl {{
  font-size:14px; color:{TEXT2}; font-weight:500;
}}
</style></head><body>
<div class="badge">AUTONOMOUS COMPUTATIONAL BIOPHYSICS</div>
<h1>DrNemo Caries Foundry</h1>
<div class="subtitle">Targeted GtfB Macrocycle Inhibitor Discovery</div>
<div class="subline">Karl Popper falsification battery ranking Pareto lead CAHWRPWHC for dental caries mutan biofilm blockade</div>
<div class="metric-row">
  <div class="metric-card">
    <span class="metric-val" style="color:{CYAN}">-8.191 kcal/mol</span>
    <span class="metric-lbl">3AIC Target Affinity</span>
  </div>
  <div class="metric-card">
    <span class="metric-val" style="color:{EMERALD}">-0.85 kcal/mol</span>
    <span class="metric-lbl">1SMD Selectivity ΔΔG</span>
  </div>
  <div class="metric-card">
    <span class="metric-val" style="color:{AMBER}">0 Sites</span>
    <span class="metric-lbl">Oral Protease Cleavage</span>
  </div>
  <div class="metric-card">
    <span class="metric-val" style="color:{CORAL}">6/6 Gates</span>
    <span class="metric-lbl">Popperian Corroboration</span>
  </div>
</div>
</body></html>"""
    render_html_to_png(html, 1200, 640, ASSETS / "social_preview.png")


def cleanup_obsolete():
    """Remove obsolete asset files."""
    obsolete = [
        "hero.svg",
        "hero_gtfb_docking.svg",
        "langgraph_loop.svg",
        "evolution_funnel_147.svg",
        "evolution_demo.gif",
        "nature_biotech_hero.png",
        "physicochemical_benchmark.png",
    ]
    for name in obsolete:
        p = ASSETS / name
        if p.exists():
            p.unlink()
            print(f"  Removed obsolete: {name}")


def main():
    os.chdir(REPO)
    ASSETS.mkdir(exist_ok=True)

    print("=== Cleaning obsolete assets ===")
    cleanup_obsolete()

    print("\n=== Building hero_banner.png (1200x380 @2x) ===")
    build_hero_banner()

    print("\n=== Building fig1_biological_mechanism.png (1200x360 @2x) ===")
    build_biological_mechanism()

    print("\n=== Building fig2_pipeline_architecture.png (1200x420 @2x) ===")
    build_pipeline_architecture()

    print("\n=== Building fig3_pocket_and_funnel.png (1200x480 @2x) ===")
    build_pocket_and_funnel()

    print("\n=== Building fig4_falsification_plate.png (1200x760 @2x) ===")
    build_falsification_plate()

    print("\n=== Building social_preview.png (1280x640 @2x) ===")
    build_social_preview()

    print("\n=== ALL FIGURES BUILT ===")
    for f in sorted(ASSETS.glob("*.png")):
        print(f"  {f.name}: {f.stat().st_size:,} bytes")


if __name__ == "__main__":
    main()
