import streamlit as st
import pandas as pd
from pathlib import Path
import streamlit.components.v1 as components
import json
import subprocess
import os
import random
import re
from dotenv import load_dotenv

load_dotenv()

BASE_DIR = Path.home() / "drnemo_research"
RESULTS_DIR = BASE_DIR / "results"
SCORES_FILE = RESULTS_DIR / "master_scores.csv"
ANALYSIS_DIR = RESULTS_DIR / "analysis"
COMPLETED_DIR = BASE_DIR / "jobs" / "completed"
PENDING_DIR = BASE_DIR / "jobs" / "pending"
TARGETS_FILE = BASE_DIR / "targets.json"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def is_duplicate(peptide_seq, project, df_scores=None):
    """Check if a peptide has already been folded or is queued."""
    job_id = f"{project}_variant_{peptide_seq}"
    # Check pending queue
    if PENDING_DIR.exists():
        pending_file = PENDING_DIR / f"{job_id}.json"
        if pending_file.exists():
            return True, "already queued in pending"
    # Check historical scores
    if df_scores is not None and not df_scores.empty:
        if job_id in df_scores["Job_ID"].values:
            return True, "already folded (in master_scores.csv)"
    return False, ""

def create_job_json(peptide_seq, project, targets_data, note=""):
    """Create a pending job JSON file. Returns (success, message)."""
    if project not in targets_data:
        return False, f"Project {project} not in targets.json"
    job_id = f"{project}_variant_{peptide_seq}"
    job = {
        "job_id": job_id,
        "protein_sequence": targets_data[project].get("sequence", ""),
        "peptide_sequence": peptide_seq,
        "project": project,
        "notes": note or f"Manual submission from dashboard."
    }
    PENDING_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PENDING_DIR / f"{job_id}.json"
    with open(out_path, 'w') as f:
        json.dump(job, f, indent=4)
    return True, job_id

st.set_page_config(page_title="Dr. Nemo: In-Silico Dental Lab", layout="wide", initial_sidebar_state="expanded")

# Custom CSS for a more "Lab" look
st.markdown("""
    <style>
    .main { background-color: #0e1117; }
    stDataFrame { border-radius: 10px; }
    .stTabs [data-baseweb="tab-list"] { gap: 24px; }
    .stTabs [data-baseweb="tab"] { height: 50px; white-space: pre-wrap; background-color: #1e1e1e; border-radius: 4px 4px 0px 0px; padding-left: 10px; padding-right: 10px; }
    .workflow-card { background-color: #1e1e1e; padding: 15px; border-left: 5px solid #2ecc71; border-radius: 5px; margin-bottom: 10px; }
    .seq-box { font-family: monospace; letter-spacing: 2px; background-color: #1e1e1e; padding: 10px; border-radius: 5px; margin-bottom: 10px; border: 1px solid #454545; color: #f1c40f; }
    </style>
    """, unsafe_allow_html=True)

colA, colB = st.columns([2, 1])
with colA:
    st.title("🧬 Dr. Nemo: In-Silico Lab")
with colB:
    st.markdown("<br><br>", unsafe_allow_html=True)
    search_term = st.text_input("🔍 Global Data Search:", placeholder="Filter by Molecule ID...").strip()

# --- Load Targets ---
target_list = []
targets_data = {}
if TARGETS_FILE.exists():
    with open(TARGETS_FILE, 'r') as f:
        targets_data = json.load(f)
        target_list = list(targets_data.keys())

df = None
if SCORES_FILE.exists():
    try:
        df = pd.read_csv(SCORES_FILE)
        df = df.sort_values(by="Timestamp", ascending=False)
        
        if 'Project' not in df.columns or df['Project'].isnull().all():
            def dict_infer_project(job_id):
                job_str = str(job_id)
                if "_variant_" in job_str:
                    return job_str.split("_variant_")[0]
                for t in target_list:
                    if t in job_str: return t
                return job_str
            df.insert(1, "Project", df["Job_ID"].apply(dict_infer_project))

        legacy_mapping = {
            "S_mutans_Baseline": "S_mutans_SpaP_Baseline",
            "P_gingivalis_Protease": "P_gingivalis_RgpB_Gingipain",
            "smutans": "S_mutans_SpaP_Baseline",
            "smutans_adhesion_blocker_01": "S_mutans_SpaP_Baseline",
            "sample_smutans": "S_mutans_SpaP_Baseline",
            "S_mutans_GtfB_Glucan": "S_mutans_GtfB_Insoluble"
        }
        df["Project"] = df["Project"].map(lambda p: legacy_mapping.get(p, p))
        
        for hp in df["Project"].unique().tolist():
            if hp not in target_list and pd.notna(hp):
                target_list.append(hp)
    except Exception as e:
        df = None

def save_targets():
    with open(TARGETS_FILE, 'w') as f:
        json.dump(targets_data, f, indent=4)

# --- SIDEBAR INTERFACE ---
teacher_mode = st.sidebar.toggle("🎓 Teacher Mode", value=False)

st.sidebar.header("🔬 Project Control")
if not target_list:
    st.sidebar.warning("No targets found. Please add a new Project Target.")
    selected_project = None
else:
    selected_project = st.sidebar.selectbox("Select Active Project Target:", target_list)

with st.sidebar.expander("⚙️ Manage Projects", expanded=False):
    m_action = st.radio("Action:", ["Add New Target", "Edit Selected"], horizontal=True, label_visibility="collapsed")
    if m_action == "Add New Target":
        with st.form("add_project_form"):
            new_id = st.text_input("Target ID (No Spaces, e.g. Virus_Target):")
            new_name = st.text_input("Display Name:")
            new_desc = st.text_area("Function Description:")
            new_seq = st.text_input("Target Protein Sequence:")
            new_pep = st.text_input("Base Peptide (Optional):")
            if st.form_submit_button("Create Project Target", use_container_width=True):
                if new_id and new_seq:
                    targets_data[new_id] = {
                        "name": new_name,
                        "description": new_desc,
                        "sequence": new_seq,
                        "base_peptide": new_pep
                    }
                    save_targets()
                    st.rerun()
                else:
                    st.error("Target ID and Sequence are required.")
    elif m_action == "Edit Selected" and selected_project:
        with st.form("edit_project_form"):
            edit_name = st.text_input("Display Name:", value=targets_data[selected_project].get("name", ""))
            edit_desc = st.text_area("Function Description:", value=targets_data[selected_project].get("description", ""))
            edit_seq = st.text_input("Target Protein Sequence:", value=targets_data[selected_project].get("sequence", ""))
            edit_pep = st.text_input("Base Peptide:", value=targets_data[selected_project].get("base_peptide", ""))
            if st.form_submit_button("Save Changes", use_container_width=True):
                if edit_seq:
                    targets_data[selected_project]["name"] = edit_name
                    targets_data[selected_project]["description"] = edit_desc
                    targets_data[selected_project]["sequence"] = edit_seq
                    targets_data[selected_project]["base_peptide"] = edit_pep
                    save_targets()
                    st.rerun()
                else:
                    st.error("Target Sequence cannot be empty.")

if teacher_mode:
    with st.sidebar.expander("🎓 Teacher Insight: Targets", expanded=True):
        st.info("A 'Project Target' is usually an essential protein that a bacteria relies on. If we block it, we stop the disease.")

st.sidebar.markdown(f"**Target Info:**")
if selected_project:
    st.sidebar.info(f"{targets_data[selected_project].get('description', '')}")

st.sidebar.markdown("---")
st.sidebar.header("🧠 Research Intelligence Control")

model_choice = st.sidebar.selectbox("Active Reasoning Model:", [
    "gpt-4o", 
    "GPT-5.4 (Cloud Scale)",
    "gemini-3.1-pro", 
    "claude-haiku", 
    "llama3.1:latest", 
    "medllama2", 
    "gemma2:27b", 
    "llama3:latest"
], index=0)
os.environ["OLLAMA_MODEL"] = model_choice
st.sidebar.markdown("---")
with st.sidebar.expander("🚀 Laboratory Actions", expanded=True):
    if teacher_mode:
        st.info("👨‍🏫 *Teacher Note:* AI algorithms generate new chemical ideas. 'Mutation Scanner' acts randomly. 'Directed Evolution' mathematically combines the best prior ideas.")
        
    st.markdown("**💡 Phase 1: Broad Survey**")
    st.info("Most promising leads so far: **S. mutans SpaP** & **P. gingivalis Protease**")
    
    col_btn1, col_btn2 = st.columns(2)
    with col_btn1:
        if st.button("🚀 Start Broad Survey", use_container_width=True, help="Seeds all 12 targets and starts the autonomous fold-evolve loop"):
            try:
                subprocess.Popen(["python3", str(BASE_DIR / "run_broad_survey.py")])
                st.success("Started!")
            except Exception as e:
                st.error(f"Error: {e}")
    with col_btn2:
        if st.button("📋 Run Phase Report", use_container_width=True, help="Check progress across all targets"):
            try:
                subprocess.Popen(["x-terminal-emulator", "-e", f"python3 {BASE_DIR / 'phase_report.py'} ; read -p 'Press Enter to close...'"])
                st.success("Running!")
            except Exception as e:
                # Fallback to just running it invisibly if no terminal emulator
                subprocess.Popen(["python3", str(BASE_DIR / "phase_report.py")])
                st.success("Running in background.")
                
    st.markdown("---")
    st.markdown("**🧪 Single Target Evolution**")
    
    if st.button("Run Mutation Scanner", use_container_width=True):
        try:
            subprocess.Popen(["python3", str(BASE_DIR / "mutation_scanner.py"), selected_project])
            st.success("Dispatched!")
        except Exception as e:
            st.error(f"Error: {e}")

    if st.button("Run Directed Evolution", use_container_width=True):
        try:
            subprocess.Popen(["python3", str(BASE_DIR / "directed_evolution.py"), selected_project])
            st.success("Dispatched to Local Engine!")
        except Exception as e:
            st.error(f"Error: {e}")
            
    st.markdown("---")
    st.markdown("**🔥 High-Throughput Engine**")
    
    col_mass1, col_mass2 = st.columns(2)
    with col_mass1:
        if st.button("💥 Generate 1625 Samples", use_container_width=True, help="Queue extensive mutation sampling"):
            try:
                subprocess.Popen(["python3", str(BASE_DIR / "mutation_scanner.py"), selected_project, "1625"])
                st.success("Dispatched 1625 jobs!")
            except Exception as e:
                st.error(f"Error: {e}")
                
    with col_mass2:
        if st.button("🔄 Resume Retrieval", use_container_width=True, help="Re-ignite Orchestrator processing if the system died"):
            try:
                subprocess.Popen(["x-terminal-emulator", "-e", f"python3 {BASE_DIR / 'orchestrator.py'} ; read -p 'Finished...'"])
                st.success("Orchestrator Resumed!")
            except Exception:
                subprocess.Popen(["python3", str(BASE_DIR / "orchestrator.py")])
                st.success("Resuming quietly...")

    st.markdown("---")
    if teacher_mode:
        st.info("👨‍🏫 *Teacher Note:* DiffDock tests small chemical drugs against the protein receptor pocket.")

# --- MANUAL PEPTIDE INPUT ---
st.sidebar.markdown("---")
with st.sidebar.expander("🧪 Manual Peptide Input", expanded=False):
    st.caption("Submit your own custom peptide sequence to OpenFold3.")
    manual_peptide = st.text_input("Peptide Sequence (amino acid letters):", max_chars=30, placeholder="e.g. LKKLLKWL").strip().upper()
    manual_note = st.text_input("Notes (optional):", placeholder="My custom design")
    
    valid_input = True
    if manual_peptide:
        if len(manual_peptide) < 5:
            st.warning("Minimum 5 amino acids.")
            valid_input = False
        elif not all(c in VALID_AA for c in manual_peptide):
            bad = [c for c in manual_peptide if c not in VALID_AA]
            st.error(f"Invalid characters: {', '.join(set(bad))}. Use only: {' '.join(sorted(VALID_AA))}")
            valid_input = False
        else:
            dup, dup_reason = is_duplicate(manual_peptide, selected_project, df)
            if dup:
                st.warning(f"⚠️ Duplicate: {dup_reason}")
                valid_input = False
            else:
                st.success(f"✅ Novel sequence for {selected_project}")
    
    if st.button("🚀 Submit to OpenFold3", use_container_width=True, disabled=(not manual_peptide or not valid_input)):
        ok, msg = create_job_json(manual_peptide, selected_project, targets_data, manual_note)
        if ok:
            st.success(f"Job created: `{msg}`")
            st.balloons()
        else:
            st.error(msg)

    if st.button("Run DiffDock Prototype", use_container_width=True):
        try:
            subprocess.Popen(["python3", str(BASE_DIR / "diffdock_pipeline.py"), selected_project])
            st.success("Dispatched!")
        except Exception as e:
            st.error(f"Error: {e}")

st.sidebar.markdown("---")
st.sidebar.header("📖 Scientific Legend")
st.sidebar.markdown("""
**pLDDT (0-100)**: Predicted Local Distance Difference Test. Measures confidence in the 3D folding structure.
- **>90**: High Accuracy
- **70-90**: Good Backup
- **<50**: Low Confidence

**iPTM (0-1)**: Interface Predicted Template Modeling score. Measures confidence in the *binding affinity* between the pathogen protein and our blocker peptide. 
- **>0.8**: Very High Binding Probability
- **0.5-0.8**: Promising Interaction
""")

# --- Load and Search ---
if df is None:
    st.warning("No data found or master_scores.csv does not exist yet. Please start the 24/7 Orchestrator.")
    st.stop()

# --- LIVE QUEUE & EUREKA DETECTION ---
st.sidebar.markdown("---")
st.sidebar.header("🔄 Live Queue Status")
if PENDING_DIR.exists():
    pending_jobs = list(PENDING_DIR.glob("*.json"))
    if pending_jobs:
        st.sidebar.info(f"⏳ **{len(pending_jobs)}** jobs currently in pending queue.")
    else:
        st.sidebar.success("✅ Queue is empty.")
else:
    st.sidebar.success("✅ Queue is empty.")

st.sidebar.markdown("**Recent Processed Jobs:**")
recent_finished = df.head(5)
for _, row in recent_finished.iterrows():
    status_icon = "✅" if row['Status'] == 'SUCCESS' else "❌"
    color = "#2ecc71" if row['Status'] == 'SUCCESS' else "#e74c3c"
    st.sidebar.markdown(f"<span style='color:{color}'>{status_icon} {str(row['Job_ID'])[:22]}...</span>", unsafe_allow_html=True)

# Eureka Detection
df['pLDDT_num'] = pd.to_numeric(df['pLDDT'], errors='coerce')
df['iPTM_num'] = pd.to_numeric(df['iPTM'], errors='coerce')
df['Composite_Score_num'] = (df['pLDDT_num'] / 100) * df['iPTM_num']

eureka_hits = df[(df['pLDDT_num'] >= 99.0) | (df['iPTM_num'] >= 0.99) | (df['Composite_Score_num'] >= 0.90)]
if not eureka_hits.empty:
    st.balloons()
    best_hit = eureka_hits.iloc[0]
    st.markdown(f"""
    <div style="background: linear-gradient(90deg, #ff007f, #7928ca); padding: 20px; border-radius: 12px; border: 2px solid #fff; animation: pulse 2s infinite; text-align: center; margin-bottom: 20px;">
        <h1 style="color: white; margin: 0; font-size: 3em;">🎉 EUREKA!! 🎉</h1>
        <h3 style="color: white; margin-top: 10px;">We hit a near-perfect confidence score!</h3>
        <p style="color: #fff; font-size: 1.2em;"><strong>Target:</strong> {best_hit.get('Job_ID')} <br> <strong>Scores:</strong> pLDDT: {best_hit['pLDDT_num']} | iPTM: {best_hit['iPTM_num']} | Composite: {best_hit['Composite_Score_num']:.4f}</p>
        <a href="mailto:roisasha@gmail.com?subject=EUREKA! Perfect Sequence Found in Dr. Nemo Lab&body=We hit a perfect score!%0A%0AJob ID: {best_hit.get('Job_ID')}%0ApLDDT: {best_hit['pLDDT_num']}%0AiPTM: {best_hit['iPTM_num']}" 
           style="background-color: white; color: #7928ca; padding: 10px 20px; border-radius: 5px; text-decoration: none; font-weight: bold; display: inline-block; margin-top: 10px;">📩 Click to send Email Notification to roisasha@gmail.com</a>
    </div>
    """, unsafe_allow_html=True)



def extract_peptide(job_id):
    if "variant_" in str(job_id):
        return str(job_id).split("variant_")[-1]
    return "Base Target"

if "Peptide_Sequence" not in df.columns:
    df.insert(2, "Peptide_Sequence", df["Job_ID"].apply(extract_peptide))

def extract_peptide(job_id):
    if "variant_" in str(job_id):
        return str(job_id).split("variant_")[-1]
    return "Base Target"

if "Peptide_Sequence" not in df.columns:
    df.insert(2, "Peptide_Sequence", df["Job_ID"].apply(extract_peptide))

project_df = df[df["Project"] == selected_project]
    
if search_term:
    project_df = project_df[project_df["Job_ID"].str.contains(search_term, case=False, na=False)]

def colorize_sequence(seq):
    mapping = {
        'A': '#e67e22', 'I': '#e67e22', 'L': '#e67e22', 'M': '#e67e22', 'F': '#e67e22', 'W': '#e67e22', 'Y': '#e67e22', 'V': '#e67e22', # Hydrophobic (Orange)
        'S': '#2ecc71', 'T': '#2ecc71', 'N': '#2ecc71', 'Q': '#2ecc71', # Polar (Green)
        'R': '#3498db', 'H': '#3498db', 'K': '#3498db', # Positive (Blue)
        'D': '#e74c3c', 'E': '#e74c3c', # Negative (Red)
        'C': '#95a5a6', 'G': '#95a5a6', 'P': '#95a5a6'  # Special (Gray)
    }
    colored_html = ""
    for char in seq:
        color = mapping.get(char.upper(), "#ffffff")
        colored_html += f"<span style='color:{color}; font-weight:bold;' title='Amino Acid: {char}'>{char}</span>"
    return colored_html

def get_rating_stars(val, type="plddt"):
    try:
        v = float(val)
        if type == "plddt":
            if v >= 90: return "⭐⭐⭐⭐⭐"
            if v >= 70: return "⭐⭐⭐"
            if v >= 50: return "⭐"
        elif type == "iptm": # iptm
            if v >= 0.8: return "⭐⭐⭐⭐⭐"
            if v >= 0.5: return "⭐⭐⭐"
        elif type == "composite":
            if v >= 0.72: return "⭐⭐⭐⭐⭐"
            if v >= 0.42: return "⭐⭐⭐"
        return ""
    except: return ""

def apply_lab_styling(df_to_style):
    def score_styler(val, type="plddt"):
        try:
            v = float(val)
            if type == "plddt":
                if v >= 90: return 'background-color: #3498db; color: white;' # Blue
                if v >= 70: return 'background-color: #2ecc71; color: black;' # Green
                if v < 50: return 'background-color: #e74c3c; color: white;' # Red
            elif type == "iptm": # iptm
                if v >= 0.8: return 'background-color: #3498db; color: white;' # Blue
                if v >= 0.5: return 'background-color: #2ecc71; color: black;' # Green
                if v < 0.3: return 'background-color: #e74c3c; color: white;' # Red
            elif type == "composite":
                if v >= 0.72: return 'background-color: #3498db; color: white;' # Blue
                if v >= 0.42: return 'background-color: #2ecc71; color: black;' # Green
        except: pass
        return ''
    
    return df_to_style.style.map(lambda x: score_styler(x, "plddt"), subset=['pLDDT'])\
                            .map(lambda x: score_styler(x, "iptm"), subset=['iPTM'])\
                            .map(lambda x: score_styler(x, "composite"), subset=['Composite_Score'] if 'Composite_Score' in df_to_style.columns else [])

# ----------------- TABS -----------------
tab1, tab2, tab5, tab3, tab6, tab4 = st.tabs(["🔬 Target Monitoring", "🏆 Leaderboard", "📈 Expansion Status", "🧬 3D Structure Viewer", "🎯 Target Reference", "📚 Manual"])

with tab3:
    st.header("🧬 3D Structure Viewer")
    st.caption("Interactive macromolecular visualization and sequence mapping for top candidates.")
    
    if SCORES_FILE.exists():
        df_3d = pd.read_csv(SCORES_FILE)
        
        col_f1, col_f2 = st.columns(2)
        with col_f1:
            default_proj_idx = target_list.index(selected_project) if selected_project in target_list else 0
            tab3_proj = st.selectbox("🎯 Target Project:", target_list, index=default_proj_idx, key="tab3_proj")
        with col_f2:
            tab3_metric = st.radio("🏆 Sort Top 10 by:", ["Composite Score", "pLDDT", "iPTM"], horizontal=True)
            
        df_3d = df_3d[(df_3d['Status'] == 'SUCCESS') & (df_3d['Project'] == tab3_proj)].copy()
        
        if not df_3d.empty:
            if tab3_metric == "Composite Score" and 'Composite_Score' in df_3d.columns:
                df_3d['Composite_Score'] = pd.to_numeric(df_3d['Composite_Score'], errors='coerce')
                df_3d = df_3d.sort_values(by="Composite_Score", ascending=False)
            elif tab3_metric == "pLDDT" and 'pLDDT' in df_3d.columns:
                df_3d['pLDDT'] = pd.to_numeric(df_3d['pLDDT'], errors='coerce')
                df_3d = df_3d.sort_values(by="pLDDT", ascending=False)
            elif tab3_metric == "iPTM" and 'iPTM' in df_3d.columns:
                df_3d['iPTM'] = pd.to_numeric(df_3d['iPTM'], errors='coerce')
                df_3d = df_3d.sort_values(by="iPTM", ascending=False)
                
            top_3d = df_3d.head(10)
            
            options = []
            job_mapping = {}
            for idx, row in top_3d.reset_index(drop=True).iterrows():
                if tab3_metric == "Composite Score": val = row.get("Composite_Score", 0)
                elif tab3_metric == "pLDDT": val = row.get("pLDDT", 0)
                else: val = row.get("iPTM", 0)
                
                lbl = f"#{idx+1} ⭐ {val} | {row.get('Peptide_Sequence', 'N/A')} ({row['Job_ID']})"
                options.append(lbl)
                job_mapping[lbl] = row['Job_ID']
                
            st.markdown("---")
            selected_lbl = st.selectbox("🔬 Select Candidate from Top 10:", options)
            
            if selected_lbl:
                sel_job = job_mapping[selected_lbl]
                pdb_path = RESULTS_DIR / f"{sel_job}.pdb"
                
                col_3d, col_seq = st.columns([1.3, 1])
                
                with col_3d:
                    if pdb_path.exists():
                        st.success(f"Rendering 3D Fold: `{sel_job}.pdb`")
                        import json
                        with open(pdb_path, "r") as f:
                            safe_pdb_code = json.dumps(f.read())

                        html_code = f"""
                        <div id="container_{sel_job}" style="height: 550px; width: 100%; position: relative;"></div>
                        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                        <script>
                            function render3D_{sel_job}() {{
                                if (typeof $3Dmol === 'undefined') {{ setTimeout(render3D_{sel_job}, 100); return; }}
                                let element = document.getElementById('container_{sel_job}');
                                let viewer = $3Dmol.createViewer(element, {{backgroundColor: '#111111'}});
                                let pdbData = {safe_pdb_code};
                                viewer.addModel(pdbData, "pdb");
                                viewer.setStyle({{}}, {{stick: {{colorscheme: 'greenCarbon', radius: 0.15}} }});
                                viewer.addStyle({{}}, {{cartoon: {{color: 'spectrum', ribbon:true}} }});
                                viewer.zoomTo();
                                viewer.spin("y", 0.5);
                                viewer.render();
                            }}
                            render3D_{sel_job}();
                        </script>
                        """
                        st.components.v1.html(html_code, height=580)
                    else:
                        st.error(f"PDB file not found: {pdb_path}")
                        
                with col_seq:
                    st.subheader(f"🎯 Visualizing Candidate")
                    st.caption(f"Job ID: `{sel_job}`")
                    
                    json_file = COMPLETED_DIR / f"{sel_job}.json"
                    seq_target = targets_data[tab3_proj].get("sequence", "SEQUENCE_DATA_UNAVAILABLE")
                    seq_peptide = extract_peptide(sel_job)
                    
                    if json_file.exists():
                        with open(json_file, 'r') as f:
                            import json
                            j_data = json.load(f)
                            seq_target = j_data.get("protein_sequence", seq_target)
                            seq_peptide = j_data.get("peptide_sequence", seq_peptide)
                            
                    st.markdown("**Sequence Legend:** 🟠 Hydrophobic | 🟢 Polar | 🔵 Positive (+) | 🔴 Negative (-)")
                    st.markdown(f"**Target Protein:**<div class='seq-box'>{colorize_sequence(seq_target)}</div>", unsafe_allow_html=True)
                    st.markdown(f"**Blocker Molecule:**<div class='seq-box'>{colorize_sequence(seq_peptide)}</div>", unsafe_allow_html=True)
                    
                    st.info("The 3D model on the left shows the physical geometry of this interaction. Refer to the electrostatic colors above to understand why the blocker adheres tightly to the target.")
        else:
            st.warning(f"No successful simulations found yet for {tab3_proj}.")
    else:
        st.warning("Master scores file not found.")


with tab6:
    st.header("🧬 Target Protein Reference Library")
    st.caption("These are fixed, well-characterized UniProt protein sequences. They NEVER change — only the blocker peptide varies.")
    
    if target_list:
        ref_data = []
        for t_id in target_list:
            t_info = targets_data.get(t_id, {})
            ref_data.append({
                "Project": t_id,
                "Protein Name": t_info.get("name", "N/A"),
                "Category": t_info.get("category", "N/A"),
                "UniProt ID": t_info.get("uniprot", "N/A"),
                "Seq Length": len(t_info.get("sequence", "")),
                "Base Peptide": t_info.get("base_peptide", "N/A")
            })
        
        df_ref = pd.DataFrame(ref_data)
        st.dataframe(df_ref, hide_index=True, use_container_width=True)
        
        st.subheader("🔎 Sequence Viewer")
        view_target = st.selectbox("Select a target to view its full sequence:", target_list, key="ref_viewer")
        if view_target:
            t_info = targets_data.get(view_target, {})
            st.markdown(f"**{t_info.get('name', view_target)}** — UniProt: `{t_info.get('uniprot', 'N/A')}` — {len(t_info.get('sequence', ''))} AA")
            st.markdown(f"<div class='seq-box' style='font-size: 12px; overflow-wrap: break-word; white-space: normal;'>{colorize_sequence(t_info.get('sequence', ''))}</div>", unsafe_allow_html=True)
            st.caption(f"📝 {t_info.get('description', '')}")
    else:
        st.warning("No targets loaded.")

with tab5:
    st.markdown("## 📈 Expansion Status")
    STATS_FILE = BASE_DIR / "expansion_stats.json"
    
    # Live auto-refresh every 30 seconds
    import time as _time
    if 'last_refresh' not in st.session_state:
        st.session_state.last_refresh = _time.time()
    
    # Read live sample count directly from CSV
    live_sample_count = 0
    if SCORES_FILE.exists():
        try:
            live_sample_count = sum(1 for _ in open(SCORES_FILE)) - 1  # subtract header
        except:
            pass
    
    num_targets = len(target_list) if target_list else 13
    total_target = 1625 * num_targets
    live_progress = live_sample_count / total_target if total_target > 0 else 0
    
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("📊 Live Samples", f"{live_sample_count:,}")
    with c2:
        st.metric("🎯 Target", f"{total_target:,}")
    with c3:
        st.metric("📈 Completion", f"{live_progress*100:.1f}%")
    with c4:
        pending_count = len(list(PENDING_DIR.glob('*.json'))) if PENDING_DIR.exists() else 0
        st.metric("⏳ Pending Jobs", pending_count)
    
    st.progress(min(live_progress, 1.0), text=f"Expansion Progress: {live_progress*100:.1f}% ({live_sample_count:,} / {total_target:,})")
    
    if STATS_FILE.exists():
        with open(STATS_FILE, 'r') as f:
            stats = json.load(f)
        st.info(f"🧬 **Lab Status:** {stats.get('status', 'UNKNOWN')} | Avg Composite: {stats.get('avg_score', 0):.3f} | Last Sync: {stats.get('last_update', 'N/A')}")
    
    st.markdown("---")
    col_a, col_b = st.columns(2)
    with col_a:
        if st.button("🚀 Start 1625-Sample Expansion", use_container_width=True, type="primary", help="Launches autonomous_lab.py targeting 1625 samples per project"):
            try:
                subprocess.Popen(["python3", str(BASE_DIR / "autonomous_lab.py")])
                st.success("✅ Autonomous Lab launched! Monitor progress here.")
                st.balloons()
            except Exception as e:
                st.error(f"Error: {e}")
    with col_b:
        if st.button("♻️ Refresh Status", use_container_width=True):
            st.rerun()
    
    # Auto-refresh every 30 seconds
    if _time.time() - st.session_state.last_refresh > 30:
        st.session_state.last_refresh = _time.time()
        st.rerun()

with tab1:
    # ── Live Telemetry Bar ──
    pending_count_live = len(list(PENDING_DIR.glob('*.json'))) if PENDING_DIR.exists() else 0
    completed_count_live = len(list(COMPLETED_DIR.glob('*.json'))) if COMPLETED_DIR.exists() else 0
    results_count_live = len(list(RESULTS_DIR.glob('*.pdb'))) if RESULTS_DIR.exists() else 0
    success_count = len(df[df['Status'] == 'SUCCESS']) if df is not None else 0
    total_processed = len(df) if df is not None else 0
    success_rate = (success_count / total_processed * 100) if total_processed > 0 else 0
    
    tm1, tm2, tm3, tm4 = st.columns(4)
    with tm1:
        st.metric("⏳ Pending", f"{pending_count_live:,}", help="Jobs waiting in queue")
    with tm2:
        st.metric("✅ PDB Results", f"{results_count_live:,}", help="Successfully folded structures")
    with tm3:
        st.metric("📊 Success Rate", f"{success_rate:.1f}%", help="Successful / Total processed")
    with tm4:
        # Estimate throughput from timestamps
        if total_processed > 1 and df is not None:
            try:
                timestamps = pd.to_datetime(df['Timestamp'], errors='coerce').dropna()
                if len(timestamps) > 1:
                    elapsed_hrs = (timestamps.max() - timestamps.min()).total_seconds() / 3600
                    throughput = total_processed / elapsed_hrs if elapsed_hrs > 0 else 0
                    st.metric("⚡ Throughput", f"{throughput:.0f}/hr")
                else:
                    st.metric("⚡ Throughput", "—")
            except:
                st.metric("⚡ Throughput", "—")
        else:
            st.metric("⚡ Throughput", "—")
    
    if pending_count_live > 0:
        live_pct = pending_count_live / (pending_count_live + results_count_live) if (pending_count_live + results_count_live) > 0 else 0
        st.progress(1.0 - live_pct, text=f"🔬 Active Processing: {pending_count_live} pending / {results_count_live} completed")
    
    st.markdown("---")
    st.write(f"### 📊 Processed Generation Log ({len(project_df)} results for {selected_project})")
    
    if teacher_mode:
         st.info("👨‍🏫 **Teacher Note:** This table is your 'Lab Notebook'. It records every experiment. Status=SUCCESS means the physics simulation completed. iPTM is the score you want to maximize.")
    else:
         st.info("💡 Select one or MULTIPLE rows below to visualize or export the combined data.")
    
    # Add star rating columns
    project_df = project_df.copy()
    project_df['pLDDT'] = pd.to_numeric(project_df['pLDDT'], errors='coerce')
    project_df['iPTM'] = pd.to_numeric(project_df['iPTM'], errors='coerce')
    
    # Calculate Composite Score for the log too
    project_df['Composite_Score'] = (project_df['pLDDT'] / 100) * project_df['iPTM']
    
    project_df["pLDDT_Rating"] = project_df["pLDDT"].apply(lambda x: get_rating_stars(x, "plddt"))
    project_df["iPTM_Rating"] = project_df["iPTM"].apply(lambda x: get_rating_stars(x, "iptm"))
    project_df["Composite_Rating"] = project_df["Composite_Score"].apply(lambda x: get_rating_stars(x, "composite"))
    
    # Define columns to display - REMOVING 'Project' from view
    view_cols = ["Timestamp", "Job_ID", "Peptide_Sequence", "pLDDT", "pLDDT_Rating", "iPTM", "iPTM_Rating", "Composite_Score", "Composite_Rating", "Status"]
    
    event = st.dataframe(
        apply_lab_styling(project_df[view_cols]),
        width="stretch",
        hide_index=True,
        on_select="rerun",
        selection_mode="multi-row"
    )

    selected_rows = event.selection.rows
    valid_jobs = project_df[project_df["Status"] == "SUCCESS"]

    if valid_jobs.empty:
        st.warning(f"No successful folds yet for project: {selected_project}")
    else:
        st.markdown("---")
        
        selected_job_ids = []
        if selected_rows:
            selected_job_ids = project_df.iloc[selected_rows]["Job_ID"].tolist()
            selected_df = project_df.iloc[selected_rows].copy()
            
            if len(selected_job_ids) >= 2:
                st.markdown("---")
                st.subheader(f"🔬 Peptide Comparison ({len(selected_job_ids)} selected)")
                
                # Extract peptide sequences
                seqs = [extract_peptide(jid) for jid in selected_job_ids]
                labels = [jid.split('_variant_')[-1] if '_variant_' in jid else jid for jid in selected_job_ids]
                
                # Alignment grid — pad to equal length
                max_len = max(len(s) for s in seqs)
                padded = [s.ljust(max_len, '-') for s in seqs]
                
                # Property classification
                aa_property = {
                    'A': 'Hydrophobic', 'I': 'Hydrophobic', 'L': 'Hydrophobic', 'M': 'Hydrophobic',
                    'F': 'Hydrophobic', 'W': 'Hydrophobic', 'Y': 'Hydrophobic', 'V': 'Hydrophobic',
                    'S': 'Polar', 'T': 'Polar', 'N': 'Polar', 'Q': 'Polar',
                    'R': 'Positive', 'H': 'Positive', 'K': 'Positive',
                    'D': 'Negative', 'E': 'Negative',
                    'C': 'Special', 'G': 'Special', 'P': 'Special', '-': 'Gap'
                }
                prop_color = {
                    'Hydrophobic': '#e67e22', 'Polar': '#2ecc71', 'Positive': '#3498db',
                    'Negative': '#e74c3c', 'Special': '#95a5a6', 'Gap': '#555'
                }
                
                # Build alignment HTML
                positions = list(range(1, max_len + 1))
                
                # Find consensus at each position
                consensus = []
                for pos_i in range(max_len):
                    residues = [s[pos_i] for s in padded if s[pos_i] != '-']
                    if residues:
                        from collections import Counter
                        consensus.append(Counter(residues).most_common(1)[0][0])
                    else:
                        consensus.append('-')
                
                # Build table
                header_row = '<th style="padding:6px 10px; background:#1a1a2e; color:#00ffcc; font-size:12px;">Position</th>'
                for p in positions:
                    header_row += f'<th style="padding:6px 8px; background:#1a1a2e; color:#aaa; font-size:11px; min-width:30px; text-align:center;">{p}</th>'
                header_row += '<th style="padding:6px 10px; background:#1a1a2e; color:#00ffcc;">pLDDT</th><th style="padding:6px 10px; background:#1a1a2e; color:#00ffcc;">iPTM</th><th style="padding:6px 10px; background:#1a1a2e; color:#00ffcc;">Composite</th>'
                
                body_rows = ""
                for idx, (seq, label) in enumerate(zip(padded, labels)):
                    row_data = selected_df.iloc[idx]
                    row = f'<td style="padding:6px 10px; color:#fff; font-weight:600; font-size:11px; white-space:nowrap;">{label}</td>'
                    for pos_i, aa in enumerate(seq):
                        prop = aa_property.get(aa, 'Gap')
                        bg = prop_color.get(prop, '#555')
                        # Highlight if matches consensus
                        match = aa == consensus[pos_i]
                        border = 'border: 2px solid #00ffcc;' if match and len(padded) > 1 else ''
                        row += f'<td style="text-align:center; padding:4px; background:{bg}33; color:{bg}; font-weight:bold; font-size:14px; {border}" title="{prop}">{aa}</td>'
                    plddt_val = row_data.get('pLDDT', 0)
                    iptm_val = row_data.get('iPTM', 0)
                    comp_val = row_data.get('Composite_Score', 0)
                    row += f'<td style="text-align:center; padding:6px; color:#fff;">{plddt_val:.1f}</td>'
                    row += f'<td style="text-align:center; padding:6px; color:#fff;">{iptm_val:.4f}</td>'
                    row += f'<td style="text-align:center; padding:6px; color:#fff; font-weight:bold;">{comp_val:.4f}</td>'
                    body_rows += f'<tr style="border-bottom:1px solid #333;">{row}</tr>'
                
                # Consensus row
                cons_row = '<td style="padding:6px 10px; color:#00ffcc; font-weight:800; font-size:12px;">CONSENSUS</td>'
                for pos_i, aa in enumerate(consensus):
                    prop = aa_property.get(aa, 'Gap')
                    bg = prop_color.get(prop, '#555')
                    cons_row += f'<td style="text-align:center; padding:4px; background:#00ffcc22; color:#00ffcc; font-weight:bold; font-size:14px; border:2px solid #00ffcc;">{aa}</td>'
                cons_row += '<td colspan="3" style="text-align:center; color:#00ffcc; font-style:italic;">Best-of pattern</td>'
                body_rows += f'<tr style="border-top: 2px solid #00ffcc;">{cons_row}</tr>'
                
                alignment_html = f'''
                <div style="overflow-x:auto; border-radius:12px; border:1px solid #333; background:#0e1117;">
                    <table style="border-collapse:collapse; width:100%; font-family: monospace;">
                        <thead><tr>{header_row}</tr></thead>
                        <tbody>{body_rows}</tbody>
                    </table>
                </div>
                '''
                st.markdown(alignment_html, unsafe_allow_html=True)
                st.caption("🟠 Hydrophobic | 🟢 Polar | 🔵 Positive | 🔴 Negative — Green borders = matches consensus")
                
                # ============ BREEDING ALGORITHM ============
                st.markdown("---")
                st.subheader("🧬 Breed Hybrid Candidates")
                st.caption("Generate new peptides by combining the best features of your selected sequences.")
                
                if st.button("🧬 Generate 5 Hybrid Candidates", use_container_width=True, type="primary"):
                    composites = [selected_df.iloc[i].get('Composite_Score', 0) for i in range(len(seqs))]
                    hybrids = set()
                    
                    # 1. Positional Best — pick AA from the parent with highest composite at each position
                    positional_best = []
                    for pos_i in range(max_len):
                        best_idx = 0
                        best_score = -1
                        for s_idx, seq in enumerate(padded):
                            if seq[pos_i] != '-' and composites[s_idx] > best_score:
                                best_score = composites[s_idx]
                                best_idx = s_idx
                        positional_best.append(padded[best_idx][pos_i])
                    hybrids.add(''.join(positional_best).replace('-', ''))
                    
                    # 2. Consensus motif directly
                    hybrids.add(''.join(consensus).replace('-', ''))
                    
                    # 3-5. Random crossovers
                    attempts = 0
                    while len(hybrids) < 5 and attempts < 50:
                        attempts += 1
                        p1, p2 = random.sample(range(len(padded)), min(2, len(padded)))
                        crossover_pt = random.randint(1, max_len - 1)
                        child = padded[p1][:crossover_pt] + padded[p2][crossover_pt:]
                        # Introduce 1 random mutation at 20% chance
                        child_list = list(child)
                        if random.random() < 0.2:
                            mut_pos = random.randint(0, len(child_list) - 1)
                            child_list[mut_pos] = random.choice(list(VALID_AA))
                        child_seq = ''.join(child_list).replace('-', '')
                        if child_seq and child_seq not in [s.replace('-', '') for s in padded]:  # Not same as any parent
                            hybrids.add(child_seq)
                    
                    # Display and submit
                    hybrids_list = list(hybrids)
                    submitted = 0
                    skipped_dup = 0
                    results_html = ""
                    for h_seq in hybrids_list[:5]:
                        dup, dup_reason = is_duplicate(h_seq, selected_project, df)
                        if dup:
                            results_html += f'<div style="padding:8px; margin:4px 0; background:#e74c3c22; border-left:3px solid #e74c3c; border-radius:4px;">⚠️ <code>{h_seq}</code> — SKIPPED ({dup_reason})</div>'
                            skipped_dup += 1
                        else:
                            ok, msg = create_job_json(h_seq, selected_project, targets_data, f"Bred from: {', '.join(labels)}")
                            if ok:
                                results_html += f'<div style="padding:8px; margin:4px 0; background:#2ecc7122; border-left:3px solid #2ecc71; border-radius:4px;">✅ <code>{h_seq}</code> — Queued as <code>{msg}</code></div>'
                                submitted += 1
                            else:
                                results_html += f'<div style="padding:8px; margin:4px 0; background:#e74c3c22; border-left:3px solid #e74c3c; border-radius:4px;">❌ <code>{h_seq}</code> — Error: {msg}</div>'
                    
                    st.markdown(results_html, unsafe_allow_html=True)
                    if submitted > 0:
                        st.success(f"🧬 {submitted} hybrid(s) queued! {skipped_dup} duplicate(s) skipped.")
                        st.balloons()
                    elif skipped_dup > 0:
                        st.warning(f"All {skipped_dup} candidates were duplicates. Select different parents or modify manually.")
                


        # Top item visualization (Default or Selected)
        job_to_inspect = selected_job_ids[0] if selected_job_ids else valid_jobs.iloc[0]["Job_ID"]
        
        st.subheader(f"🎯 Visualizing Candidate: {job_to_inspect}")
        
        # Load Sequence Data from JSON
        json_file = COMPLETED_DIR / f"{job_to_inspect}.json"
        seq_target = targets_data[selected_project].get("sequence", "SEQUENCE_DATA_UNAVAILABLE") # Fallback to prevent empty headers
        seq_peptide = extract_peptide(job_to_inspect)
        
        if json_file.exists():
            with open(json_file, 'r') as f:
                j_data = json.load(f)
                seq_target = j_data.get("protein_sequence", seq_target)
                seq_peptide = j_data.get("peptide_sequence", seq_peptide)
                
        col1, col2 = st.columns([1.5, 1])
        
        with col1:
            st.markdown("**Sequence Legend:** 🟠 Hydrophobic | 🟢 Polar | 🔵 Positive (+) | 🔴 Negative (-)")
            if teacher_mode:
                st.info("👨‍🏫 **Teacher Note:** Amino acids have different chemical properties. Notice the colors below! These magnetic and chemical charges determine how the molecules stick together in 3D space.")
            
            st.markdown(f"**Target Protein:**<div class='seq-box'>{colorize_sequence(seq_target)}</div>", unsafe_allow_html=True)
            st.markdown(f"**Blocker Molecule:**<div class='seq-box'>{colorize_sequence(seq_peptide)}</div>", unsafe_allow_html=True)
            
            pdb_path = RESULTS_DIR / f"{job_to_inspect}.pdb"
            if pdb_path.exists():
                with open(pdb_path, 'r') as f:
                    pdb_data = f.read().replace("'", "\\'").replace("\\n", "\\\\n").replace("\n", "\\n")
                
                html_code = f"""
                <div id="container" style="height: 550px; width: 100%; position: relative;"></div>
                <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                <script>
                    function render3D() {{
                        if (typeof $3Dmol === 'undefined') {{ setTimeout(render3D, 100); return; }}
                        let element = document.getElementById('container');
                        let view = $3Dmol.createViewer(element, {{backgroundColor: "#111111"}});
                        let data = '{pdb_data}';
                        view.addModel(data, "pdb");
                        
                        let colorfunc = function(atom) {{
                            if(atom.chain === 'B') return "#f1c40f"; // Peptide
                            
                            let b = atom.b;
                            if (b < 50) return "#e74c3c"; // Red
                            if (b < 70) return "#f39c12"; // Orange/Yellow
                            if (b < 90) return "#2ecc71"; // Green (High Confidence)
                            return "#3498db"; // Cyan/Blue (Very High Confidence)
                        }};
                        
                        view.setStyle({{chain: 'A'}}, {{cartoon: {{colorfunc: colorfunc}}}});
                        view.setStyle({{chain: 'B'}}, {{stick: {{color: "#f1c40f", radius: 0.2}}, sphere: {{color: "#f1c40f", radius: 0.8}}}});
                        
                        view.zoomTo();
                        view.render();
                    }}
                    render3D();
                </script>
                """
                components.html(html_code, height=550)
                st.caption("🟢 Green/Blue: High Confidence Target Fold | 🔴 Red: Low Confidence Frame | 🟡 Designed Peptide Interactor")
            else:
                placeholder_css = f'''
                <div style="height: 550px; width: 100%; border-radius: 15px; border: 2px dashed #00ffcc; display: flex; align-items: center; justify-content: center; background: radial-gradient(circle, #0e1117 0%, #000000 100%); overflow: hidden; position: relative;">
                    <div style="position: absolute; width: 200%; height: 200%; background: conic-gradient(from 0deg, transparent 0%, rgba(0,255,204,0.1) 30%, transparent 40%); animation: radar 4s linear infinite;"></div>
                    <div style="text-align: center; font-family: sans-serif; position: relative; z-index: 10; animation: pulse 2s infinite;">
                        <span style="font-size: 80px; display: block; margin-bottom: 20px;">🧊</span>
                        <h3 style="color: #00ffcc; margin: 0; letter-spacing: 2px; text-transform: uppercase;">In-Silico Render Engine Ready</h3>
                        <p style="color: #95a5a6; margin-top: 10px; font-size: 14px;">Awaiting 3D fold geometry for<br><strong style="color: #fff;">{job_to_inspect}</strong></p>
                    </div>
                </div>
                <style>
                @keyframes pulse {{ 0% {{ opacity: 0.8; transform: scale(0.98); }} 50% {{ opacity: 1; transform: scale(1.02); }} 100% {{ opacity: 0.8; transform: scale(0.98); }} }}
                @keyframes radar {{ 100% {{ transform: rotate(360deg); }} }}
                </style>
                '''
                components.html(placeholder_css, height=560)
                
        with col2:
            st.subheader("MedGemma 4B Analysis", help="Our local LLM reads the raw physics scores and turns them into a clinical summary.")
            
            if teacher_mode:
                st.info("👨‍🏫 **Teacher Note:** MedGemma is an Artificial Intelligence trained on medical textbooks. It acts as our Lab Diagnostician, verifying the physics results.")
                
            analysis_path = ANALYSIS_DIR / f"{job_to_inspect}.md"
            if analysis_path.exists():
                with open(analysis_path, 'r', encoding="utf-8") as f:
                    content = f.read()
                    st.markdown(f'<div style="background-color: #262730; padding: 20px; border-radius: 10px; border: 1px solid #454545;">{content}</div>', unsafe_allow_html=True)
                
                st.markdown("<br>", unsafe_allow_html=True)
                pitch_report = f"# Institutional Research Pitch: {job_to_inspect}\\n\\n## Computational Finding\\nWe have successfully resolved a high-confidence structural interaction for targeting the pathogen mechanism.\\n\\n### Molecular Construction\\n- **Target Sequence:** {seq_target}\\n- **Inhibitor Sequence:** {seq_peptide}\\n\\n{content}\\n\\n## Protocol Requirement\\nRequesting immediate validation of this candidate via structural biochemistry and in-vitro bacterial adhesion assays.\\n"
                st.download_button("📥 Export IRB Pitch Deck (.md)", data=pitch_report, file_name=f"IRB_Pitch_{job_to_inspect}.md", mime="text/markdown", type="primary", use_container_width=True)
                
            else:
                st.info("Analysis pending. The local MedGemma model may still be generating this report.")

with tab2:
    st.markdown(f"## 🏆 Leaderboard: {selected_project}")
    st.write("Tracking the highest-confidence molecule variants achieved so far in this project.")
    
    if teacher_mode:
        st.info("👨‍🏫 **Teacher Note:** This chart is the proof that Directed Evolution works. As we breed new chemicals, you should see the line trend upwards towards 1.0.")
        
    if not valid_jobs.empty:
        valid_jobs = valid_jobs.copy()
        valid_jobs['iPTM'] = pd.to_numeric(valid_jobs['iPTM'], errors='coerce')
        valid_jobs['pLDDT'] = pd.to_numeric(valid_jobs['pLDDT'], errors='coerce')
        
        # Calculate Composite Score (Logical AND gate filtering)
        valid_jobs['Composite_Score'] = (valid_jobs['pLDDT'] / 100) * valid_jobs['iPTM']
        import numpy as np
        valid_jobs['Composite_Score'] = valid_jobs['Composite_Score'].replace([np.inf, -np.inf, np.nan], 0)
        
        # Add ratings to leaderboard
        valid_jobs["pLDDT_Rating"] = valid_jobs["pLDDT"].apply(lambda x: get_rating_stars(x, "plddt"))
        valid_jobs["iPTM_Rating"] = valid_jobs["iPTM"].apply(lambda x: get_rating_stars(x, "iptm"))
        valid_jobs["Composite_Rating"] = valid_jobs["Composite_Score"].apply(lambda x: get_rating_stars(x, "composite"))
        
        leaderboard = valid_jobs.sort_values(by="Composite_Score", ascending=False).head(10)
        
        leaderboard = leaderboard.copy()
        
        l_cols = ["Job_ID", "Peptide_Sequence", "pLDDT", "pLDDT_Rating", "iPTM", "iPTM_Rating", "Composite_Score", "Composite_Rating", "Timestamp"]
        
        st.dataframe(
            apply_lab_styling(leaderboard[l_cols]),
            width="stretch",
            hide_index=True
        )
        
        st.markdown("### 📈 Composite Score Progression Over Time")
        trend_data = valid_jobs.sort_values(by="Timestamp")
        st.line_chart(trend_data.set_index("Timestamp")["Composite_Score"])
    else:
        st.write("No successful jobs plotted yet to rank.")
        
    # ============ BULK EXPORT ENGINE ============
    st.markdown("---")
    st.subheader("📚 Synthesis-Ready Export")
    export_top_n = st.slider("Select Top N Candidates to Export:", min_value=1, max_value=min(200, len(valid_jobs) if not valid_jobs.empty else 1), value=min(10, len(valid_jobs) if not valid_jobs.empty else 1))
    
    if not valid_jobs.empty:
        export_df = valid_jobs.sort_values(by="Composite_Score", ascending=False).head(export_top_n)
        
        import io
        import zipfile
        
        csv_data = export_df.to_csv(index=False).encode('utf-8')
        
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
            for _, row in export_df.iterrows():
                jid = row['Job_ID']
                pdb_path = RESULTS_DIR / f"{jid}.pdb"
                if pdb_path.exists():
                    zip_file.write(pdb_path, f"{jid}.pdb")
        zip_data = zip_buffer.getvalue()
        
        ecol1, ecol2 = st.columns(2)
        with ecol1:
            st.download_button("📊 Export Target Data (.csv)", data=csv_data, file_name=f"Top_{export_top_n}_Candidates_{selected_project}.csv", mime="text/csv", type="primary", use_container_width=True)
        with ecol2:
            st.download_button("🧬 Download PDB Structures (.zip)", data=zip_data, file_name=f"Top_{export_top_n}_Structures_{selected_project}.zip", mime="application/zip", type="primary", use_container_width=True)
    else:
        st.info("No candidates available to export.")
    
    # ============ PAN-BIOFILM CROSS-TEST ============
    st.markdown("---")
    with st.expander("🌐 Pan-Biofilm Cross-Test", expanded=False):
        st.caption("Test a top peptide against ALL other target proteins to discover universal multi-pathogen inhibitors.")
        
        if teacher_mode:
            st.info("👨‍🏫 **Teacher Note:** A peptide built to block SpaP might also block GtfB or even C. albicans. Cross-testing reveals 'broad-spectrum' candidates — the holy grail of antimicrobial peptides.")
        
        if not valid_jobs.empty:
            cross_candidates = valid_jobs.sort_values(by="Composite_Score", ascending=False).head(10)
            cross_options = [f"{row['Job_ID']} (Composite: {row['Composite_Score']:.4f})" for _, row in cross_candidates.iterrows()]
            cross_mapping = {opt: row['Job_ID'] for opt, (_, row) in zip(cross_options, cross_candidates.iterrows())}
            
            selected_cross = st.selectbox("Select candidate to cross-test:", cross_options, key="cross_test_sel")
            
            if st.button("🧬 Run Pan-Biofilm Cross-Test", use_container_width=True, type="primary"):
                cross_job_id = cross_mapping[selected_cross]
                # Extract peptide sequence
                cross_peptide = extract_peptide(cross_job_id)
                # Get the source project
                source_project = cross_job_id.rsplit('_variant_', 1)[0] if '_variant_' in cross_job_id else selected_project
                
                cross_queued = 0
                cross_skipped = 0
                for t_id in target_list:
                    if t_id == source_project:
                        continue  # skip the original project
                    new_job_id = f"{t_id}_variant_{cross_peptide}"
                    new_job_path = PENDING_DIR / f"{new_job_id}.json"
                    if new_job_path.exists():
                        cross_skipped += 1
                        continue
                    # Check if already in scores
                    if df is not None and new_job_id in df['Job_ID'].values:
                        cross_skipped += 1
                        continue
                    
                    protein_seq = targets_data[t_id].get('sequence', '')
                    if not protein_seq:
                        continue
                    
                    job_data = {
                        "job_id": new_job_id,
                        "protein_sequence": protein_seq,
                        "peptide_sequence": cross_peptide,
                        "project": t_id,
                        "notes": f"Pan-Biofilm Cross-Test from {source_project}"
                    }
                    PENDING_DIR.mkdir(parents=True, exist_ok=True)
                    with open(new_job_path, 'w') as f:
                        json.dump(job_data, f, indent=4)
                    cross_queued += 1
                
                if cross_queued > 0:
                    st.success(f"🌐 Queued {cross_queued} cross-tests across {len(target_list) - 1} targets! ({cross_skipped} duplicates skipped)")
                    st.balloons()
                else:
                    st.warning(f"All cross-tests already exist or were skipped ({cross_skipped} duplicates).")
        else:
            st.info("No candidates available for cross-testing.")

with tab4:
    st.markdown("## 📚 Dr. Nemo Pipeline Guide")
    st.write("A visual guide to the exact assembly line process of your virtual laboratory.")
    
    html_operator = """
    <style>
        .flow-container {
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 20px;
            padding: 40px;
            font-family: 'Inter', sans-serif;
            background: transparent;
        }
        .flow-item {
            width: 320px;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            color: white;
            font-weight: 600;
            box-shadow: 0 4px 15px rgba(0,0,0,0.3);
            border: 2px solid rgba(255,255,255,0.1);
        }
        .flow-item.step { background: #2980b9; border-color: #3498db; }
        .flow-item.queue { background: #8e44ad; border-color: #9b59b6; border-radius: 50px; }
        .flow-item.auto { background: #27ae60; border-color: #2ecc71; }
        .flow-item.review { background: #d35400; border-color: #e67e22; }
        .arrow { font-size: 24px; color: rgba(255,255,255,0.4); margin: -5px 0; }
        .sub-text { display: block; font-size: 13px; font-weight: 400; margin-top: 8px; opacity: 0.9; }
        .edge-label { font-size: 12px; color: #bdc3c7; background: #1e1e1e; padding: 2px 8px; border-radius: 10px; }
    </style>
    <div class="flow-container">
        <div class="flow-item step">👤 STEP 1: Select Target Project<span class="sub-text">Choose your biological target via the sidebar.</span></div>
        <div class="arrow">↓</div>
        <div class="flow-item step">🧬 STEP 2: Generate Blocker Candidates<span class="sub-text">Directed Evolution breeds new peptide variants from top performers.</span></div>
        <div class="edge-label">Creates 10+ Jobs</div>
        <div class="arrow">↓</div>
        <div class="flow-item queue">📂 jobs/pending/<span class="sub-text">Queue fills with folding tasks.</span></div>
        <div class="edge-label">4 Parallel Workers</div>
        <div class="arrow">↓</div>
        <div class="flow-item auto">⚙️ STEP 3: OpenFold3 Pipeline<span class="sub-text">4 NVIDIA API keys fold structures in parallel. Each key = 1 dedicated worker.</span></div>
        <div class="edge-label">Saves PDB + Scores</div>
        <div class="arrow">↓</div>
        <div class="flow-item queue">📊 master_scores.csv<span class="sub-text">Central database of all molecule fitness data.</span></div>
        <div class="arrow">↓</div>
        <div class="flow-item review">👤 STEP 4: Review Leaderboard<span class="sub-text">Inspect Composite Scores. Blue/Green = high binding confidence.</span></div>
        <div class="edge-label">If no perfect hits → loop back</div>
        <div class="arrow">↓</div>
        <div class="flow-item step">🔁 STEP 5: Evolve Best Hits<span class="sub-text">Select top candidates → Breed hybrids → Repeat until Eureka.</span></div>
    </div>
    """
    components.html(html_operator, height=1000, scrolling=True)
    
    st.markdown("---")
    st.markdown("""    
**Scientific Legend:**
- **pLDDT (0-100)**: Confidence in 3D fold accuracy. >90 = High, 70-90 = Good, <50 = Low.
- **iPTM (0-1)**: Confidence in binding affinity. >0.8 = Excellent, 0.5-0.8 = Promising.
- **Composite Score**: `(pLDDT/100) × iPTM`. Our unified fitness metric.
- **Target Protein**: The fixed microbial protein we want to block (UniProt verified, never changes).
- **Blocker Peptide**: The variable molecule we are evolving to bind and inhibit the target.
    """)

