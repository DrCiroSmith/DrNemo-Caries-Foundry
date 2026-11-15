#!/usr/bin/env python3
import os
import shutil
import subprocess
import time
from pathlib import Path
from playwright.sync_api import sync_playwright

OUTPUT_DIR = Path("/mnt/c/dev/job_packets/linkedin_drnemo")
HTML_PATH = Path("/home/drnemo/projects/portfolio_overhaul/DrNemo-Caries-Foundry/docs/index.html")
TEMP_VIDEO_DIR = Path("/tmp/playwright_linkedin_video")

# 1. Clean old synthetic/box-card slop
old_cards = [
    "01_3d_gtfb_docked_peptide_hero.png",
    "02_interactive_webgl_studio_showcase.png",
    "03_nature_style_evolution_convergence.png",
    "04_langgraph_biotech_architecture.png"
]
print("Cleaning old synthetic box cards from linkedin_drnemo...")
for card in old_cards:
    p = OUTPUT_DIR / card
    if p.exists():
        p.unlink()
        print(f"  Deleted old card: {card}")

TEMP_VIDEO_DIR.mkdir(parents=True, exist_ok=True)
for f in TEMP_VIDEO_DIR.glob("*"):
    f.unlink()

print("Launching Playwright for Retina screenshots and 1080p video recording...")
with sync_playwright() as p:
    browser = p.chromium.launch(
        headless=True,
        args=["--enable-webgl", "--ignore-gpu-blocklist", "--hide-scrollbars"]
    )

    # -------------------------------------------------------------
    # STEP 2A: Capture 2400px Retina Screenshots (device_scale_factor=2)
    # -------------------------------------------------------------
    context_retina = browser.new_context(
        viewport={"width": 1200, "height": 900},
        device_scale_factor=2
    )
    page = context_retina.new_page()
    page.goto(f"file://{HTML_PATH}")
    page.wait_for_timeout(3500) # Allow 3Dmol to load 3AIC crystal pocket & CRRRWWRRRRC

    # Capture 02_live_webgl_studio_ui_2400px.png (Full Studio UI)
    ui_shot_path = OUTPUT_DIR / "02_live_webgl_studio_ui_2400px.png"
    page.screenshot(path=str(ui_shot_path))
    print(f"Captured: {ui_shot_path.name} ({ui_shot_path.stat().st_size / 1024.0:.1f} KB)")

    # Capture 01_real_3aic_docked_pocket_2400px.png (Close-up Macro on Pocket)
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
    pocket_shot_path = OUTPUT_DIR / "01_real_3aic_docked_pocket_2400px.png"
    page.screenshot(path=str(pocket_shot_path))
    print(f"Captured: {pocket_shot_path.name} ({pocket_shot_path.stat().st_size / 1024.0:.1f} KB)")
    context_retina.close()

    # -------------------------------------------------------------
    # STEP 2B: Record 1080p Walkthrough Video (1440x1080, 16–18s)
    # -------------------------------------------------------------
    print("\nRecording live 1080p software walkthrough choreography...")
    context_video = browser.new_context(
        record_video_dir=str(TEMP_VIDEO_DIR),
        record_video_size={"width": 1440, "height": 1080},
        viewport={"width": 1440, "height": 1080}
    )
    vpage = context_video.new_page()
    vpage.goto(f"file://{HTML_PATH}")
    vpage.wait_for_timeout(2000) # allow initial load

    # Phase 1: 0–6s 3D smooth rotation of 3AIC canyon + CRRRWWRRRRC
    print("  [0-6s] 3D smooth rotation of real 3AIC crystal canyon + docked CRRRWWRRRRC...")
    vpage.evaluate("viewer.spin('y', 0.45);")
    vpage.wait_for_timeout(5500)

    # Phase 2: 6–11s Toggle electrostatics & scrub Gen 1->20 slider
    print("  [6-11s] Toggling electrostatics and scrubbing Gen 1->20 evolution slider...")
    vpage.locator("#btn-col-chain").click()
    vpage.wait_for_timeout(1000)
    
    vpage.locator("#btn-rep-surface").click()
    vpage.wait_for_timeout(1500)
    
    vpage.locator("#btn-rep-complex").click()
    vpage.wait_for_timeout(800)

    # Scrub slider
    for g in [1, 5, 10, 15, 20]:
        vpage.evaluate(f"onScrubGeneration({g}); document.getElementById('gen-slider').value = {g};")
        vpage.wait_for_timeout(450)

    # Phase 3: 11–17s Chitosan Inspector mutant test
    print("  [11-17s] Testing anionic mutant in Chitosan Inspector -> Fix to Gold Lead...")
    vpage.locator(".tab-btn:has-text('Chitosan Inspector')").click()
    vpage.wait_for_timeout(1000)

    # Test anionic mutant (click preset)
    vpage.locator("button.preset-btn:has-text('+ Asp')").click()
    vpage.wait_for_timeout(2200)

    # Fix back to Gold Lead
    vpage.locator("button.preset-btn:has-text('Gold Lead')").click()
    vpage.wait_for_timeout(2200)

    context_video.close()
    browser.close()

# -------------------------------------------------------------
# STEP 3: Encode with ffmpeg to drnemo_linkedin_demo_1080p.mp4
# -------------------------------------------------------------
webms = list(TEMP_VIDEO_DIR.glob("*.webm"))
if not webms:
    raise RuntimeError("No recorded video found in temp dir")

raw_webm = webms[0]
mp4_out = OUTPUT_DIR / "drnemo_linkedin_demo_1080p.mp4"
print(f"\nEncoding {raw_webm.name} to {mp4_out.name} via ffmpeg...")

cmd = [
    "ffmpeg", "-y",
    "-i", str(raw_webm),
    "-c:v", "libx264",
    "-preset", "slow",
    "-crf", "18",
    "-pix_fmt", "yuv420p",
    "-r", "30",
    "-movflags", "+faststart",
    str(mp4_out)
]
subprocess.run(cmd, check=True)

print(f"\nFinal Video Generated: {mp4_out}")
print(f"  Size: {mp4_out.stat().st_size / (1024.0 * 1024.0):.2f} MB")

# Inspect with ffprobe
ffprobe_cmd = [
    "ffprobe", "-v", "error",
    "-select_streams", "v:0",
    "-show_entries", "stream=width,height,duration,r_frame_rate,codec_name",
    "-of", "default=noprint_wrappers=1",
    str(mp4_out)
]
res = subprocess.run(ffprobe_cmd, capture_output=True, text=True, check=True)
print("ffprobe specs:\n" + res.stdout)
print("MEDIA CAPTURE AND ENCODING COMPLETE!")
