"""Pin every entry in manuscript/references.bib to its verified PMID, DOI, and title."""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BIB_PATH = REPO_ROOT / "manuscript" / "references.bib"

# Expected metadata for all 22 entries in manuscript/references.bib:
# (pmid, doi, exact_title_or_pubmed_title)
# Entries without a PubMed record (pre-MEDLINE, unindexed, book, software) have pmid=None.
EXPECTED_REFERENCES = {
    "Trott2010": (
        "19499576",
        "10.1002/jcc.21334",
        "AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading",
    ),
    "Eberhardt2021": (
        "34278794",
        "10.1021/acs.jcim.1c00203",
        "AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings",
    ),
    "Quiroga2016": (
        "27171006",
        "10.1371/journal.pone.0155183",
        "Vinardo: A Scoring Function Based on Autodock Vina Improves Scoring, Docking, and Virtual Screening",
    ),
    "Ito2011": (
        "21354427",
        "10.1016/j.jmb.2011.02.028",
        "Crystal structure of glucansucrase from the dental caries pathogen Streptococcus mutans",
    ),
    "Ramasubbu1996": (
        "15299664",
        "10.1107/S0907444995014119",
        "Structure of human salivary alpha-amylase at 1.6 A resolution: implications for its role in the oral cavity",
    ),
    "Stephan1944": (
        None,
        "10.1177/00220345440230040401",
        "Intra-Oral Hydrogen-Ion Concentrations Associated with Dental Caries Activity",
    ),
    "Popper1959": (
        None,
        None,
        "The Logic of Scientific Discovery",
    ),
    "Bowen2018": (
        "29097091",
        "10.1016/j.tim.2017.09.008",
        "Oral Biofilms: Pathogens, Matrix, and Polymicrobial Interactions in Microenvironments",
    ),
    "Koo2013": (
        "24045647",
        "10.1177/0022034513504218",
        "The exopolysaccharide matrix: a virulence determinant of cariogenic biofilm",
    ),
    "Eastman2017": (
        "28746339",
        "10.1371/journal.pcbi.1005659",
        "OpenMM 7: Rapid development of high performance algorithms for molecular dynamics",
    ),
    "Maier2015": (
        "26574453",
        "10.1021/acs.jctc.5b00255",
        "ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB",
    ),
    "Nguyen2013": (
        "25788871",
        "10.1021/ct3010485",
        "Improved Generalized Born Solvent Model Parameters for Protein Simulations",
    ),
    "Halgren1996": (
        None,
        "10.1002/(SICI)1096-987X(199604)17:5/6<490::AID-JCC1>3.0.CO;2-P",
        "Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94",
    ),
    "OBoyle2011": (
        "21982300",
        "10.1186/1758-2946-3-33",
        "Open Babel: An open chemical toolbox",
    ),
    "Hopkins2004": (
        "15109945",
        "10.1016/S1359-6446(04)03069-7",
        "Ligand efficiency: a useful metric for lead selection",
    ),
    "Craik2013": (
        "23253135",
        "10.1111/cbdd.12055",
        "The future of peptide-based drugs",
    ),
    "Marsh2006": (
        "15153690",
        "10.1159/000077756",
        "Dental plaque as a microbial biofilm",
    ),
    "Zero1999": (
        "15153701",
        "10.1159/000077767",
        "Sugars - the arch criminal?",
    ),
    "Burne2000": (
        "11094270",
        "10.1111/j.1574-6968.2000.tb09393.x",
        "Alkali production by oral bacteria and protection against dental caries",
    ),
    "Xiao2012": (
        "22496649",
        "10.1371/journal.ppat.1002623",
        "The exopolysaccharide matrix modulates the interaction between 3D architecture and virulence of a mixed-species oral biofilm",
    ),
    "Scharnow2019": (
        "31391878",
        "10.1039/c9md00015a",
        "Targeting S. mutans biofilms: a perspective on preventing dental caries",
    ),
    "Landrum2020": (
        None,
        None,
        "RDKit: Open-source cheminformatics software",
    ),
}


def _parse_bib(text: str) -> dict[str, dict[str, str]]:
    entries = re.split(r"\n(?=@)", text.strip())
    parsed: dict[str, dict[str, str]] = {}
    for raw in entries:
        m = re.match(r"@(\w+)\{([^,]+),", raw)
        assert m is not None, f"Could not parse BibTeX header: {raw[:60]}"
        key = m.group(2).strip()
        fields: dict[str, str] = {}
        for fm in re.finditer(r"(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}", raw):
            fields[fm.group(1).lower()] = re.sub(r"\s+", " ", fm.group(2)).strip()
        parsed[key] = fields
    return parsed


def _normalize_title(title: str) -> str:
    # Strip simple LaTeX formatting commands like \textit{...} while preserving text.
    cleaned = re.sub(r"\\textit\{([^}]*)\}", r"\1", title)
    return cleaned.strip().lower()


def test_references_bib_has_exact_22_entries():
    bib = _parse_bib(BIB_PATH.read_text(encoding="utf-8"))
    assert set(bib.keys()) == set(EXPECTED_REFERENCES.keys())
    assert len(bib) == 22


def test_references_bib_pinned_to_pmid_doi_and_title():
    bib = _parse_bib(BIB_PATH.read_text(encoding="utf-8"))
    for key, (exp_pmid, exp_doi, exp_title) in EXPECTED_REFERENCES.items():
        entry = bib[key]
        if exp_pmid is not None and "pmid" in entry:
            assert entry["pmid"] == exp_pmid, f"{key}: pmid {entry.get('pmid')} != {exp_pmid}"
        assert entry.get("doi") == exp_doi, f"{key}: doi {entry.get('doi')} != {exp_doi}"
        assert _normalize_title(entry.get("title", "")) == _normalize_title(exp_title), (
            f"{key}: title {entry.get('title')!r} != {exp_title!r}"
        )
