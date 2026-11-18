#!/usr/bin/env python3
"""
DrNemo Caries Foundry — LinkedIn Animated GIF Generator
Author: Dr. Roi Castrellon, DDS (DrNemo Caries Foundry)

Generates 4 high-fidelity animated GIFs (<8 MB each) into /mnt/c/dev/job_packets/linkedin_drnemo/:
1. 01_3d_gtfb_pocket_orbit.gif (3Dmol.js WebGL rotation of 3AIC active cleft + CAHWRPWHC)
2. 02_popper_falsification_loop.gif (6-angle Karl Popper falsification cycle)
3. 03_contender_funnel_52.gif (52-contender evolutionary funnel to Pareto champion)
4. 04_evolution_progression.gif (20-generation evolution scrubber and convergence scatter)
"""

import os
import sys
import time
import shutil
import subprocess
from pathlib import Path
from playwright.sync_api import sync_playwright

REPO_ROOT = Path(__file__).parent.parent.resolve()
DOCS_HTML = REPO_ROOT / "docs" / "index.html"
OUT_DIR = Path("/mnt/c/dev/job_packets/linkedin_drnemo")
TMP_FRAMES = Path("/tmp/foundry_gif_frames")

OUT_DIR.mkdir(parents=True, exist_ok=True)


def encode_gif_with_palette(frames_pattern: str, output_gif: Path, fps: int = 15, scale_w: int = 720):
    """Encodes frames into an ultra-clean, high-quality GIF using ffmpeg palettegen/paletteuse."""
    palette_png = output_gif.parent / f"{output_gif.stem}_palette.png"
    
    # 1. Generate palette
    cmd_pal = [
        "ffmpeg", "-y", "-framerate", str(fps),
        "-i", frames_pattern,
        "-vf", f"fps={fps},scale={scale_w}:-1:flags=lanczos,palettegen=stats_mode=diff",
        str(palette_png)
    ]
    subprocess.run(cmd_pal, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    # 2. Render GIF using palette
    cmd_gif = [
        "ffmpeg", "-y", "-framerate", str(fps),
        "-i", frames_pattern,
        "-i", str(palette_png),
        "-lavfi", f"fps={fps},scale={scale_w}:-1:flags=lanczos [x]; [x][1:v] paletteuse=dither=bayer:bayer_scale=3",
        str(output_gif)
    ]
    subprocess.run(cmd_gif, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    if palette_png.exists():
        palette_png.unlink()

    sz_mb = output_gif.stat().st_size / (1024 * 1024)
    print(f"  [OK] Generated {output_gif.name} ({sz_mb:.2f} MB, <8 MB target)")
    return sz_mb


def generate_gif_01_3d_orbit():
    print("\n--- Generating 01_3d_gtfb_pocket_orbit.gif ---")
    frames_dir = TMP_FRAMES / "orbit"
    if frames_dir.exists():
        shutil.rmtree(frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)

    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            args=[
                "--enable-webgl",
                "--ignore-gpu-blocklist",
                "--enable-gpu-rasterization",
                "--use-gl=angle",
                "--use-angle=swiftshader"
            ]
        )
        page = browser.new_page(viewport={"width": 1280, "height": 900})
        page.goto(f"file://{DOCS_HTML}")
        page.wait_for_timeout(3500)

        # Stop spin, center on ligand, and orient view
        page.evaluate("""() => {
            if (window.viewer) {
                viewer.spin(false);
                isSpinning = false;
                viewer.zoomTo({ chain: 'B' });
                viewer.zoom(0.85);
                viewer.render();
            }
        }""")
        page.wait_for_timeout(1000)

        viewport = page.locator(".molecular-stage-container")
        box = viewport.bounding_box()

        num_frames = 36
        rot_per_frame = 360.0 / num_frames

        for f in range(num_frames):
            page.evaluate(f"""() => {{
                if (window.viewer) {{
                    viewer.rotate({rot_per_frame:.3f}, 'y');
                    viewer.render();
                }}
            }}""")
            time.sleep(0.02)
            frame_path = frames_dir / f"frame_{f:04d}.png"
            if box:
                page.screenshot(path=str(frame_path), clip=box, animations="allow")
            else:
                viewport.screenshot(path=str(frame_path), animations="allow")

        browser.close()

    out_gif = OUT_DIR / "01_3d_gtfb_pocket_orbit.gif"
    encode_gif_with_palette(str(frames_dir / "frame_%04d.png"), out_gif, fps=16, scale_w=720)


def generate_gif_02_popper_loop():
    print("\n--- Generating 02_popper_falsification_loop.gif ---")
    frames_dir = TMP_FRAMES / "popper"
    if frames_dir.exists():
        shutil.rmtree(frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)

    svg_file = REPO_ROOT / "assets" / "langgraph_loop.svg"
    svg_content = svg_file.read_text(encoding="utf-8")

    # Render animated SVG frames by creating an HTML harness with an animated scanning laser / pulse
    html_harness = f"""<!DOCTYPE html>
    <html>
    <head>
    <meta charset="utf-8">
    <style>
      body {{ margin:0; padding:0; background:#0B1120; display:flex; justify-content:center; align-items:center; height:100vh; overflow:hidden; font-family:sans-serif; }}
      .svg-wrap {{ position:relative; width:1100px; height:440px; }}
      .scanline {{ position:absolute; top:0; bottom:0; width:6px; background:linear-gradient(to right, rgba(56,189,248,0), rgba(56,189,248,0.8), rgba(16,185,129,0.9), rgba(56,189,248,0)); filter:blur(1px); box-shadow:0 0 16px rgba(56,189,248,0.9); }}
      .badge-tag {{ position:absolute; bottom:20px; right:30px; background:#10B981; color:#064E3B; font-weight:800; font-size:14px; padding:6px 14px; border-radius:9999px; box-shadow:0 4px 12px rgba(16,185,129,0.4); }}
    </style>
    </head>
    <body>
      <div class="svg-wrap" id="stage">
        {svg_content}
        <div class="scanline" id="scanner"></div>
        <div class="badge-tag" id="badge">6/6 GATES CORROBORATED</div>
      </div>
    </body>
    </html>"""

    harness_path = frames_dir / "harness.html"
    harness_path.write_text(html_harness, encoding="utf-8")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1180, "height": 480})
        page.goto(f"file://{harness_path}")
        page.wait_for_timeout(500)

        num_frames = 40
        for f in range(num_frames):
            pct = (f / num_frames) * 100.0
            page.evaluate(f"""() => {{
                const scan = document.getElementById('scanner');
                scan.style.left = '{pct:.1f}%';
                const badge = document.getElementById('badge');
                if ({pct} < 45) {{
                    badge.innerText = 'SCREENING CANDIDATE: CAHWRPWHC';
                    badge.style.background = '#38BDF8';
                    badge.style.color = '#0C4A6E';
                }} else if ({pct} < 85) {{
                    badge.innerText = 'AMYLASE ΔΔG: -0.85 kcal/mol (EXCLUDED)';
                    badge.style.background = '#F59E0B';
                    badge.style.color = '#78350F';
                }} else {{
                    badge.innerText = '🏆 CORROBORATED (Passed 6 of 6 Gates)';
                    badge.style.background = '#10B981';
                    badge.style.color = '#064E3B';
                }}
            }}""")
            frame_path = frames_dir / f"frame_{f:04d}.png"
            page.locator("#stage").screenshot(path=str(frame_path))

        browser.close()

    out_gif = OUT_DIR / "02_popper_falsification_loop.gif"
    encode_gif_with_palette(str(frames_dir / "frame_%04d.png"), out_gif, fps=15, scale_w=720)


def generate_gif_03_contender_funnel():
    print("\n--- Generating 03_contender_funnel_52.gif ---")
    frames_dir = TMP_FRAMES / "funnel"
    if frames_dir.exists():
        shutil.rmtree(frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)

    svg_file = REPO_ROOT / "assets" / "evolution_funnel_147.svg"
    svg_content = svg_file.read_text(encoding="utf-8")

    html_harness = f"""<!DOCTYPE html>
    <html>
    <head>
    <meta charset="utf-8">
    <style>
      body {{ margin:0; padding:0; background:#0B1120; display:flex; justify-content:center; align-items:center; height:100vh; overflow:hidden; }}
      .svg-wrap {{ position:relative; width:1100px; height:440px; }}
      .flow-beam {{ position:absolute; top:0; width:8px; height:100%; background:linear-gradient(to bottom, transparent, #38BDF8, #10B981, transparent); filter:blur(1px); }}
    </style>
    </head>
    <body>
      <div class="svg-wrap" id="stage">
        {svg_content}
        <div class="flow-beam" id="beam"></div>
      </div>
    </body>
    </html>"""

    harness_path = frames_dir / "harness.html"
    harness_path.write_text(html_harness, encoding="utf-8")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1180, "height": 480})
        page.goto(f"file://{harness_path}")
        page.wait_for_timeout(500)

        num_frames = 40
        for f in range(num_frames):
            pct = (f / num_frames) * 100.0
            page.evaluate(f"document.getElementById('beam').style.left = '{pct:.1f}%';")
            frame_path = frames_dir / f"frame_{f:04d}.png"
            page.locator("#stage").screenshot(path=str(frame_path))

        browser.close()

    out_gif = OUT_DIR / "03_contender_funnel_52.gif"
    encode_gif_with_palette(str(frames_dir / "frame_%04d.png"), out_gif, fps=15, scale_w=720)


def generate_gif_04_evolution_progression():
    print("\n--- Generating 04_evolution_progression.gif ---")
    frames_dir = TMP_FRAMES / "evolution"
    if frames_dir.exists():
        shutil.rmtree(frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1280, "height": 900})
        page.goto(f"file://{DOCS_HTML}")
        page.wait_for_timeout(3000)

        # Switch to tab-evolution
        page.evaluate("""() => {
            if (typeof switchTab === 'function') {
                switchTab('tab-evolution');
            } else {
                document.querySelectorAll('.tab-btn').forEach(b => {
                    if (b.innerText.includes('Evolution') || b.innerText.includes('Campaign')) b.click();
                });
            }
        }""")
        page.wait_for_timeout(800)

        panel = page.locator("#tab-evolution")

        num_frames = 40
        for f in range(num_frames):
            gen = int(1 + (f / (num_frames - 1)) * 19)
            page.evaluate(f"""() => {{
                const sl = document.getElementById('gen-slider');
                if (sl) {{
                    sl.value = {gen};
                    if (typeof onScrubGeneration === 'function') onScrubGeneration({gen});
                }}
            }}""")
            time.sleep(0.03)
            frame_path = frames_dir / f"frame_{f:04d}.png"
            panel.screenshot(path=str(frame_path))

        browser.close()

    out_gif = OUT_DIR / "04_evolution_progression.gif"
    encode_gif_with_palette(str(frames_dir / "frame_%04d.png"), out_gif, fps=15, scale_w=720)


def main():
    print("=" * 70)
    print("  Exporting 4 Animated LinkedIn GIFs into /mnt/c/dev/job_packets/linkedin_drnemo/")
    print("=" * 70)
    t0 = time.time()

    generate_gif_01_3d_orbit()
    generate_gif_02_popper_loop()
    generate_gif_03_contender_funnel()
    generate_gif_04_evolution_progression()

    el = time.time() - t0
    print("=" * 70)
    print(f"All 4 GIFs successfully exported in {el:.1f}s!")
    for g in ["01_3d_gtfb_pocket_orbit.gif", "02_popper_falsification_loop.gif", "03_contender_funnel_52.gif", "04_evolution_progression.gif"]:
        p = OUT_DIR / g
        sz = p.stat().st_size / (1024 * 1024)
        print(f"  - {p} ({sz:.2f} MB)")
    print("=" * 70)


if __name__ == "__main__":
    main()
