#!/usr/bin/env python3
"""Deterministic 60fps/30fps Frame-Stepped Screen Studio Media Generator.

Renders 570 exact 30.00 fps frames (19.00s) across 6 parallel Chromium workers,
with camera zoom rig, vector cursor with click rings, lower-third caption pill,
continuous 20-gen evolution scrub, and Chitosan Inspector mutant validation.
Downscales via ffmpeg Lanczos to razor-sharp 1600x1200 H.264 mp4.
"""

import math
import os
import shutil
import subprocess
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from playwright.sync_api import sync_playwright

TOTAL_FRAMES = 570
FPS = 30
VIEWPORT_W = 1600
VIEWPORT_H = 1200
NUM_WORKERS = 6

ROOT = Path("/home/drnemo/projects/portfolio_overhaul/DrNemo-Caries-Foundry")
DOCS_HTML = ROOT / "docs" / "index.html"
FRAMES_DIR = Path("/tmp/drnemo_frames")
OUTPUT_DIR = Path("/mnt/c/dev/job_packets/linkedin_drnemo")

def get_cumulative_rot(f: int) -> float:
    if f <= 150:
        return f * 0.55
    elif f <= 300:
        return 150 * 0.55 + (f - 150) * 0.25
    elif f <= 465:
        return 150 * 0.55 + 150 * 0.25
    else:
        return 150 * 0.55 + 150 * 0.25 + (f - 465) * 0.35

def render_chunk(chunk):
    chunk_idx, start_f, end_f = chunk
    print(f"[Worker {chunk_idx}] Starting frames {start_f:04d} -> {end_f:04d} ({end_f - start_f} frames)...")
    
    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            args=[
                "--enable-webgl",
                "--ignore-gpu-blocklist",
                "--enable-gpu-rasterization",
                "--hide-scrollbars"
            ]
        )
        page = browser.new_page(viewport={"width": VIEWPORT_W, "height": VIEWPORT_H}, device_scale_factor=1.0)
        page.goto(f"file://{DOCS_HTML}")
        page.wait_for_timeout(3200)

        # Inject CSS Rig, Cursor, and Caption Pill (without DOM re-parenting!)
        page.evaluate("""() => {
            viewer.spin(false);
            const main = document.querySelector('main');
            main.style.transformOrigin = '0 0';
            main.style.willChange = 'transform';
            document.body.style.overflow = 'hidden';
            document.body.style.margin = '0';

            // Add macOS Vector Cursor
            const cur = document.createElement('div');
            cur.id = 'demo-cursor';
            cur.style.cssText = 'position:fixed;pointer-events:none;z-index:999999;left:0;top:0;transform:translate3d(-100px,-100px,0);';
            cur.innerHTML = `
                <svg width="28" height="28" viewBox="0 0 24 24" fill="none" style="filter:drop-shadow(0 3px 6px rgba(0,0,0,0.5));">
                    <path d="M5.5 3.5l11.5 8.5-5.5 1.5-2.5 5.5-3.5-15.5z" fill="#0F172A" stroke="#FFFFFF" stroke-width="1.8" stroke-linejoin="round"/>
                </svg>
                <div id="demo-click-ring" style="position:absolute;top:0;left:0;width:32px;height:32px;border:2.5px solid rgba(13,148,136,0.9);border-radius:50%;transform:translate(-50%,-50%) scale(0);opacity:0;transition:none;"></div>
            `;
            document.body.appendChild(cur);

            // Add Caption Pill
            const pill = document.createElement('div');
            pill.id = 'demo-caption-pill';
            pill.style.cssText = 'position:fixed;bottom:32px;left:50%;transform:translateX(-50%);background:rgba(15,23,42,0.92);backdrop-filter:blur(12px);border:1px solid rgba(255,255,255,0.18);border-radius:9999px;padding:12px 28px;color:#F8FAFC;font-family:Inter,-apple-system,sans-serif;font-size:16px;font-weight:600;letter-spacing:-0.2px;box-shadow:0 12px 36px rgba(0,0,0,0.42);pointer-events:none;z-index:999999;display:flex;align-items:center;gap:8px;';
            pill.innerHTML = '<span id="demo-caption-text"></span>';
            document.body.appendChild(pill);
        }""")

        # Fast forward 3Dmol rotation to initial start_f state
        initial_rot = get_cumulative_rot(start_f)
        if initial_rot > 0:
            page.evaluate(f"viewer.rotate({initial_rot:.4f}, 'y'); viewer.render();")

        t_start = time.time()
        for f in range(start_f, end_f):
            delta_rot = 0.55 if f == 0 else (get_cumulative_rot(f) - get_cumulative_rot(f - 1))

            page.evaluate(f"""() => {{
                const f = {f};
                const main = document.querySelector('main');
                const header = document.querySelector('header');
                const footer = document.querySelector('footer');
                const cur = document.getElementById('demo-cursor');
                const ring = document.getElementById('demo-click-ring');
                const cap = document.getElementById('demo-caption-text');

                const easeInOutCubic = t => t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
                const lerp = (a, b, t) => a + (b - a) * t;

                let scale = 1.0;
                let tx = 800;
                let ty = 600;
                let cx = -100;
                let cy = -100;
                let ringScale = 0;
                let ringOpacity = 0;
                let caption = "";
                let headerOpacity = 1.0;
                let activeTab = 'evolution';
                let activeGen = 20;
                let activePreset = 'CRRRWWRRRRC';

                if (f <= 150) {{
                    // ACT 1: Full-Bleed 3D Crystallography Hook
                    scale = 1.38;
                    tx = 565.6;
                    ty = 585.2;
                    headerOpacity = 0.0;
                    caption = "S. mutans GtfB Pocket (RCSB 3AIC) · Evolved Cyclic Lead CPHWRPWHC (-14.02 kcal/mol)";
                    cx = -100;
                    cy = -100;
                }} else if (f <= 300) {{
                    // ACT 2: Ease Out to Full Studio & Continuous 20-Gen Scrub
                    caption = "Karl Popper 6-Angle Falsification · AutoDock Vina v1.2.7 Physics & Statistical Nulls";
                    if (f <= 185) {{
                        const t = (f - 150) / 35.0;
                        const e = easeInOutCubic(t);
                        scale = lerp(1.38, 1.0, e);
                        tx = lerp(565.6, 800.0, e);
                        ty = lerp(585.2, 600.0, e);
                        headerOpacity = e;
                        if (f >= 170) {{
                            const tc = (f - 170) / 15.0;
                            const ec = easeInOutCubic(tc);
                            cx = lerp(800.0, 1205.0, ec);
                            cy = lerp(600.0, 346.0, ec);
                        }}
                        activeGen = 1;
                    }} else if (f <= 285) {{
                        scale = 1.0;
                        tx = 800.0;
                        ty = 600.0;
                        headerOpacity = 1.0;
                        const t = (f - 186) / 99.0;
                        activeGen = Math.max(1, Math.min(20, Math.floor(1 + t * 19.99)));
                        cx = lerp(1205.0, 1430.0, t);
                        cy = 346.0;
                        ringScale = 0.85;
                        ringOpacity = 0.45;
                    }} else {{
                        scale = 1.0;
                        tx = 800.0;
                        ty = 600.0;
                        headerOpacity = 1.0;
                        activeGen = 20;
                        const t = (f - 285) / 15.0;
                        const e = easeInOutCubic(t);
                        cx = lerp(1430.0, 1292.7, e);
                        cy = lerp(346.0, 209.8, e);
                    }}
                }} else if (f <= 465) {{
                    // ACT 3: Zoom Into Chitosan Hydrogel Formulation Shield
                    caption = "Oral Protease Shield & Chitosan Gate · Rejects Lysis, Precipitation & Traps Pre-Folding";
                    activeTab = 'inspector';
                    if (f <= 335) {{
                        const t = (f - 300) / 35.0;
                        const e = easeInOutCubic(t);
                        scale = lerp(1.0, 1.42, e);
                        tx = lerp(800.0, 1290.2, e);
                        ty = lerp(600.0, 541.4, e);
                        // Cursor on Tab
                        cx = 800 + (1292.7 - tx) * scale;
                        cy = 600 + (209.8 - ty) * scale;
                        if (f >= 305 && f <= 315) {{
                            const rp = (f - 305) / 10.0;
                            ringScale = lerp(0.2, 1.4, rp);
                            ringOpacity = lerp(0.9, 0.0, rp);
                        }}
                    }} else if (f <= 395) {{
                        scale = 1.42;
                        tx = 1290.2;
                        ty = 541.4;
                        // Button '+ Asp': DOM offset (1206.0, 400.0) -> screen coordinates
                        const targetX = 800 + (1206.0 - tx) * scale;
                        const targetY = 600 + (400.0 - ty) * scale;
                        if (f <= 350) {{
                            const startX = 800 + (1292.7 - tx) * scale;
                            const startY = 600 + (209.8 - ty) * scale;
                            const t = (f - 335) / 15.0;
                            const e = easeInOutCubic(t);
                            cx = lerp(startX, targetX, e);
                            cy = lerp(startY, targetY, e);
                        }} else {{
                            cx = targetX;
                            cy = targetY;
                        }}
                        if (f >= 348 && f <= 358) {{
                            const rp = (f - 348) / 10.0;
                            ringScale = lerp(0.2, 1.4, rp);
                            ringOpacity = lerp(0.9, 0.0, rp);
                        }}
                        if (f >= 350) {{
                            activePreset = 'CRRRWWDRRRC'; // FAIL
                        }}
                    }} else {{
                        scale = 1.42;
                        tx = 1290.2;
                        ty = 541.4;
                        const prevX = 800 + (1206.0 - tx) * scale;
                        const prevY = 600 + (400.0 - ty) * scale;
                        const targetX = 800 + (1301.7 - tx) * scale;
                        const targetY = 600 + (370.0 - ty) * scale;
                        if (f <= 412) {{
                            const t = (f - 395) / 17.0;
                            const e = easeInOutCubic(t);
                            cx = lerp(prevX, targetX, e);
                            cy = lerp(prevY, targetY, e);
                        }} else {{
                            cx = targetX;
                            cy = targetY;
                        }}
                        if (f >= 410 && f <= 420) {{
                            const rp = (f - 410) / 10.0;
                            ringScale = lerp(0.2, 1.4, rp);
                            ringOpacity = lerp(0.9, 0.0, rp);
                        }}
                        if (f >= 412) {{
                            activePreset = 'CRRRWWRRRRC'; // PASS
                        }} else {{
                            activePreset = 'CRRRWWDRRRC'; // Still FAIL until clicked
                        }}
                    }}
                }} else {{
                    // ACT 4: Final Heroic View & Clean Loop
                    caption = "Open-Source MIT · Real RCSB Coordinates & 100% Empirical Physics (Link in Comments)";
                    activeTab = 'falsification';
                    activePreset = 'CRRRWWRRRRC';
                    if (f <= 505) {{
                        const t = (f - 465) / 40.0;
                        const e = easeInOutCubic(t);
                        scale = lerp(1.42, 1.08, e);
                        tx = lerp(1290.2, 800.0, e);
                        ty = lerp(541.4, 560.0, e);
                    }} else {{
                        scale = 1.08;
                        tx = 800.0;
                        ty = 560.0;
                    }}
                    cx = -100;
                    cy = -100;
                }}

                // Apply CSS Transform
                main.style.transform = `translate(${{800 - tx * scale}}px, ${{600 - ty * scale}}px) scale(${{scale}})`;
                header.style.opacity = headerOpacity;
                footer.style.opacity = headerOpacity;

                // Cursor & Ring
                cur.style.transform = `translate3d(${{cx}}px, ${{cy}}px, 0)`;
                ring.style.transform = `translate(-50%, -50%) scale(${{ringScale}})`;
                ring.style.opacity = ringOpacity;

                // Caption
                cap.innerText = caption;

                // Tab switching
                const tabInsp = document.getElementById('tab-inspector');
                const tabEvol = document.getElementById('tab-evolution');
                if (activeTab === 'inspector' && tabInsp && tabInsp.style.display !== 'flex') {{
                    switchTab('tab-inspector');
                }} else if (activeTab === 'evolution' && tabEvol && tabEvol.style.display !== 'flex') {{
                    switchTab('tab-evolution');
                }}

                // Slider scrub
                if (activeTab === 'evolution') {{
                    const slider = document.getElementById('gen-slider');
                    if (slider && parseInt(slider.value, 10) !== activeGen) {{
                        slider.value = activeGen;
                        onScrubGeneration(activeGen);
                    }}
                }}

                // Inspector preset
                if (activeTab === 'inspector') {{
                    const curInput = document.getElementById('peptide-input');
                    if (curInput && curInput.value !== activePreset) {{
                        setPreset(activePreset);
                    }}
                }}

                // 3D rotation
                if ({delta_rot:.4f} > 0) {{
                    viewer.rotate({delta_rot:.4f}, 'y');
                    viewer.render();
                }}
            }}""")

            page.screenshot(path=str(FRAMES_DIR / f"frame_{f:04d}.png"))

        browser.close()
        dt = time.time() - t_start
        print(f"[Worker {chunk_idx}] Completed {end_f - start_f} frames in {dt:.1f}s ({dt/(end_f-start_f)*1000:.1f}ms/frame)!")
        return end_f - start_f

def main():
    print("=" * 80)
    print("DETERMINISTIC 60FPS/30FPS 'SCREEN STUDIO' CINEMATIC MEDIA GENERATOR")
    print(f"Total Frames: {TOTAL_FRAMES} (19.00s @ {FPS} fps) | Workers: {NUM_WORKERS}")
    print("=" * 80)

    # 1. Clean frames directory
    FRAMES_DIR.mkdir(parents=True, exist_ok=True)
    for f in FRAMES_DIR.glob("*.png"):
        f.unlink()

    # 2. Divide frames across workers
    chunk_size = math.ceil(TOTAL_FRAMES / NUM_WORKERS)
    chunks = []
    for i in range(NUM_WORKERS):
        start = i * chunk_size
        end = min(TOTAL_FRAMES, (i + 1) * chunk_size)
        if start < TOTAL_FRAMES:
            chunks.append((i, start, end))

    print(f"Divided into {len(chunks)} chunks: {[(c[1], c[2]) for c in chunks]}")
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        results = list(executor.map(render_chunk, chunks))

    t_render = time.time() - t0
    total_rendered = sum(results)
    print(f"\nAll {total_rendered} frames rendered in {t_render:.1f}s ({t_render/total_rendered*1000:.1f}ms/frame across {NUM_WORKERS} workers)!")

    # Verify all frames exist
    missing = []
    for f in range(TOTAL_FRAMES):
        p = FRAMES_DIR / f"frame_{f:04d}.png"
        if not p.exists():
            missing.append(f)
    if missing:
        raise RuntimeError(f"Missing {len(missing)} frames: {missing[:10]}...")

    # 3. Capture 2400px Retina Screenshots with updated badge
    print("\nCapturing 2400px Retina Screenshots (01 & 02)...")
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True, args=["--enable-webgl", "--ignore-gpu-blocklist"])
        page = browser.new_page(viewport={"width": 1200, "height": 900}, device_scale_factor=2)
        page.goto(f"file://{DOCS_HTML}")
        page.wait_for_timeout(3500)

        # 02_live_webgl_studio_ui_2400px.png
        ui_shot = OUTPUT_DIR / "02_live_webgl_studio_ui_2400px.png"
        page.screenshot(path=str(ui_shot))
        print(f"  Captured: {ui_shot.name} ({ui_shot.stat().st_size / 1024.0:.1f} KB)")

        # 01_real_3aic_docked_pocket_2400px.png
        page.evaluate("""() => {
            const stage = document.querySelector('.molecular-stage-container');
            stage.style.position = 'fixed';
            stage.style.top = '0';
            stage.style.left = '0';
            stage.style.width = '100vw';
            stage.style.height = '100vh';
            stage.style.zIndex = '99999';
            viewer.resize();
            viewer.spin(false);
            viewer.zoomTo({ chain: 'B' });
            viewer.zoom(1.18);
            viewer.render();
        }""")
        page.wait_for_timeout(1000)
        pocket_shot = OUTPUT_DIR / "01_real_3aic_docked_pocket_2400px.png"
        page.screenshot(path=str(pocket_shot))
        print(f"  Captured: {pocket_shot.name} ({pocket_shot.stat().st_size / 1024.0:.1f} KB)")
        browser.close()

    # 4. High-bitrate Lanczos FFmpeg encoding
    mp4_out = OUTPUT_DIR / "drnemo_linkedin_demo_1080p.mp4"
    print(f"\nEncoding frames via ffmpeg to {mp4_out.name}...")
    cmd = [
        "ffmpeg", "-y",
        "-framerate", "30",
        "-i", str(FRAMES_DIR / "frame_%04d.png"),
        "-vf", "scale=1600:1200:flags=lanczos",
        "-c:v", "libx264",
        "-profile:v", "high",
        "-crf", "16",
        "-pix_fmt", "yuv420p",
        "-movflags", "+faststart",
        str(mp4_out)
    ]
    subprocess.run(cmd, check=True)
    print(f"Video encoded successfully! Size: {mp4_out.stat().st_size / (1024.0 * 1024.0):.2f} MB")

    # ffprobe specs
    ffprobe_cmd = [
        "ffprobe", "-v", "error",
        "-select_streams", "v:0",
        "-show_entries", "stream=width,height,duration,r_frame_rate,codec_name",
        "-of", "default=noprint_wrappers=1",
        str(mp4_out)
    ]
    res = subprocess.run(ffprobe_cmd, capture_output=True, text=True, check=True)
    print("\nffprobe specs:\n" + res.stdout)

    # 5. Extract 5 keyframes for vision audit
    print("\nExtracting 5 keyframes for vision audit...")
    keyframe_stamps = [
        ("t0_1", "0.1", "Frame 0.1s - Act 1 Crystallography Hook"),
        ("t3_5", "3.5", "Frame 3.5s - Act 1 3D Orbit Parallax"),
        ("t8_5", "8.5", "Frame 8.5s - Act 2 Continuous Evolution Scrub"),
        ("t12_5", "12.5", "Frame 12.5s - Act 3 Zoomed RED FAIL State"),
        ("t14_8", "14.8", "Frame 14.8s - Act 3 Zoomed GREEN PASS State")
    ]
    for tag, sec, desc in keyframe_stamps:
        kf_path = Path(f"/tmp/keyframe_{tag}.png")
        cmd_kf = [
            "ffmpeg", "-y",
            "-ss", sec,
            "-i", str(mp4_out),
            "-vframes", "1",
            str(kf_path)
        ]
        subprocess.run(cmd_kf, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"  Extracted {kf_path.name} ({desc})")

    print("\nCapture and keyframe extraction completed successfully!")

if __name__ == "__main__":
    main()
