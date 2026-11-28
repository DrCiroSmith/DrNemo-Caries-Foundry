"""Guards over the two bibliography files.

Twelve of the twenty-two entries in manuscript/references.bib were wrong when
this file was written, and four of them named a journal, volume, issue, page
range and DOI that belong to no paper that could be located, while carrying the
author list and title of a real paper published somewhere else. The two
bibliography files also disagreed with each other, and both were wrong.

The table below is a transcription of a field-by-field check against PubMed
carried out on 2026-11-27, with the PMID that was read recorded beside each
entry. These tests do not call PubMed: they are offline and deterministic, and
they fail if the bibliography drifts away from the state that was checked. A
new entry has to be added to the table by hand, which is the point: adding one
forces someone to go and look it up.
"""

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
REFERENCES = REPO_ROOT / "manuscript" / "references.bib"
PAPER_BIB = REPO_ROOT / "paper" / "paper.bib"
PREPRINT_TEX = REPO_ROOT / "manuscript" / "DrNemo_Preprint.tex"
PAPER_MD = REPO_ROOT / "paper" / "paper.md"

# key -> (pmid or None, journal, volume, number, pages, doi)
VERIFIED = {
    "Trott2010": ("19499576", "Journal of Computational Chemistry", "31", "2", "455--461", "10.1002/jcc.21334"),
    "Eberhardt2021": ("34278794", "Journal of Chemical Information and Modeling", "61", "8", "3891--3898", "10.1021/acs.jcim.1c00203"),
    "Quiroga2016": ("27171006", "PLOS ONE", "11", "5", "e0155183", "10.1371/journal.pone.0155183"),
    "Ito2011": ("21354427", "Journal of Molecular Biology", "408", "2", "177--186", "10.1016/j.jmb.2011.02.028"),
    "Ramasubbu1996": ("15299664", "Acta Crystallographica Section D: Biological Crystallography", "52", "3", "435--446", "10.1107/S0907444995014119"),
    "Stephan1944": (None, "Journal of Dental Research", "23", "4", "257--266", "10.1177/00220345440230040401"),
    "Bowen2018": ("29097091", "Trends in Microbiology", "26", "3", "229--242", "10.1016/j.tim.2017.09.008"),
    "Koo2013": ("24045647", "Journal of Dental Research", "92", "12", "1065--1073", "10.1177/0022034513504218"),
    "Eastman2017": ("28746339", "PLOS Computational Biology", "13", "7", "e1005659", "10.1371/journal.pcbi.1005659"),
    "Maier2015": ("26574453", "Journal of Chemical Theory and Computation", "11", "8", "3696--3713", "10.1021/acs.jctc.5b00255"),
    "Nguyen2013": ("25788871", "Journal of Chemical Theory and Computation", "9", "4", "2020--2034", "10.1021/ct3010485"),
    "Halgren1996": (None, "Journal of Computational Chemistry", "17", "5-6", "490--519", "10.1002/(SICI)1096-987X(199604)17:5/6<490::AID-JCC1>3.0.CO;2-P"),
    "OBoyle2011": ("21982300", "Journal of Cheminformatics", "3", "1", "33", "10.1186/1758-2946-3-33"),
    "Hopkins2004": ("15109945", "Drug Discovery Today", "9", "10", "430--431", "10.1016/S1359-6446(04)03069-7"),
    "Craik2013": ("23253135", "Chemical Biology \\& Drug Design", "81", "1", "136--147", "10.1111/cbdd.12055"),
    "Marsh2004": ("15153690", "Caries Research", "38", "3", "204--211", "10.1159/000077756"),
    "Zero2004": ("15153701", "Caries Research", "38", "3", "277--285", "10.1159/000077767"),
    "Burne2000": ("11094270", "FEMS Microbiology Letters", "193", "1", "1--6", "10.1111/j.1574-6968.2000.tb09393.x"),
    "Xiao2012": ("22496649", "PLOS Pathogens", "8", "4", "e1002623", "10.1371/journal.ppat.1002623"),
    "Scharnow2019": ("31391878", "MedChemComm", "10", "7", "1057--1067", "10.1039/c9md00015a"),
}

# Entries that are not journal articles, so they carry no journal or DOI.
NON_ARTICLE_KEYS = {"Popper1959", "Landrum2020"}

# Author names that were in the bibliography and are not the names of the
# authors of the cited papers. None of them may come back.
WITHDRAWN_AUTHOR_NAMES = [
    "Wu, Henny",          # PMID 29097091 gives Wu, Hui
    "Lu, Bob",            # PMID 22496649 gives Lu, Bingwen
    "Hey, Christopher",   # PMID 22496649 gives Heydorn, Arne
    "Ito, Kazuo",         # PMID 21354427 gives Ito, Keisuke
    "Ito, Shingo",        # PMID 21354427 gives Ito, Sohei
    "Scharnow, Anna M.",  # PMID 31391878 gives Scharnow, Amber M
    "Ragunath, Chandran", # not an author of PMID 15299664
    "Bhandary",           # not an author of PMID 15299664
]

# DOIs that were in the bibliography and that no located paper carries. None of
# them may come back, in any file.
WITHDRAWN_DOIS = [
    "10.1006/jmbi.1996.0603",      # claimed for Ramasubbu1996
    "10.1007/s40496-018-0168-3",   # claimed for Bowen2018
    "10.1038/nrmicro3115",         # claimed for Koo2013
    "10.1038/nrd1546",             # claimed for Hopkins2004
    "10.1159/000095697",           # claimed for Marsh, whose paper is 2004
    "10.1159/000016501",           # claimed for Zero, whose paper is 2004
    "10.1021/ct300962b",           # claimed for Nguyen2013
]


def comment_blocks(path):
    """Return {key: the % comment lines immediately above that entry}.

    A PMID has to be found in the block attached to its own entry. An earlier
    version of this check searched the whole file, and so passed when an entry's
    own PMID line was deleted, because the same PMID was quoted in another
    entry's comment as a corroborating record. That is the same defect as a
    guard pinned wider than the claim it protects.
    """
    blocks, pending = {}, []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.lstrip()
        if stripped.startswith("%"):
            pending.append(stripped)
        elif stripped.startswith("@"):
            key = re.match(r"@\w+\s*\{\s*([A-Za-z0-9_:-]+)", stripped)
            if key:
                blocks[key.group(1)] = "\n".join(pending)
            pending = []
        elif not stripped:
            continue
        else:
            pending = []
    return blocks


# The two files that name the withdrawn DOIs in order to record that they are
# withdrawn: this test's own table, and section 3.12 of the negative-result
# paper, which is the register of what each wrong entry said. Nothing else in
# the repository has a reason to contain one, and
# test_no_bibliography_entry_carries_a_withdrawn_doi checks the doi fields
# themselves, so neither file can hide a withdrawn DOI used as a citation.
SWEEP_EXEMPT = {
    Path(__file__).resolve(),
    (REPO_ROOT / "paper" / "negative_result_GTF-SI.md").resolve(),
}


def parse_bib(path):
    """Return {key: {field: value}} for a flat .bib file, ignoring % comments."""
    text = "\n".join(
        line for line in path.read_text(encoding="utf-8").splitlines()
        if not line.lstrip().startswith("%")
    )
    entries = {}
    for match in re.finditer(r"@(\w+)\s*\{\s*([A-Za-z0-9_:-]+)\s*,(.*?)\n\}", text, re.S):
        kind, key, body = match.group(1), match.group(2), match.group(3)
        fields = {"__kind__": kind.lower()}
        for fmatch in re.finditer(r"(\w+)\s*=\s*\{(.*?)\}\s*(?:,\s*\n|\n*$)", body, re.S):
            fields[fmatch.group(1).lower()] = fmatch.group(2).strip()
        entries[key] = fields
    return entries


@pytest.fixture(scope="module")
def references():
    return parse_bib(REFERENCES)


def test_every_entry_is_either_verified_or_declared_a_non_article(references):
    assert references, "no entries parsed out of references.bib"
    unaccounted = set(references) - set(VERIFIED) - NON_ARTICLE_KEYS
    assert not unaccounted, (
        "these entries are in references.bib but not in this file's verification "
        f"table: {sorted(unaccounted)}. Check each one against PubMed and add it."
    )
    absent = set(VERIFIED) - set(references)
    assert not absent, f"verified entries missing from references.bib: {sorted(absent)}"


def test_each_entry_matches_the_fields_that_were_checked(references):
    checked = 0
    for key, (_pmid, journal, volume, number, pages, doi) in VERIFIED.items():
        entry = references[key]
        assert entry["journal"] == journal, f"{key} journal"
        assert entry["volume"] == volume, f"{key} volume"
        assert entry["number"] == number, f"{key} number"
        assert entry["pages"] == pages, f"{key} pages"
        assert entry["doi"] == doi, f"{key} doi"
        checked += 1
    assert checked == 20, f"expected 20 journal entries, checked {checked}"


def test_each_pubmed_checked_entry_records_the_pmid_that_was_read():
    blocks = comment_blocks(REFERENCES)
    checked = 0
    for key, (pmid, *_rest) in VERIFIED.items():
        if pmid is None:
            continue
        assert key in blocks, f"{key} has no comment block above it"
        assert f"PMID {pmid}" in blocks[key], (
            f"the comment block above {key} no longer records PMID {pmid}, the "
            f"record that was read to check it"
        )
        checked += 1
    assert checked == 18, f"expected 18 PubMed-checked entries, checked {checked}"


def test_the_citation_key_year_matches_the_entry_year(references):
    """Marsh2006 pointed at a 2004 paper and Zero1999 at a 2004 paper.

    A key whose year disagrees with the entry invites anyone writing prose to
    cite the wrong year, so the two have to agree.
    """
    for key, entry in references.items():
        match = re.search(r"(1[89]\d\d|20\d\d)$", key)
        if not match:
            continue
        assert match.group(1) == entry["year"], (
            f"citation key {key} says {match.group(1)} but its year field says "
            f"{entry['year']}"
        )


def test_no_withdrawn_doi_appears_anywhere_in_the_repository():
    searched = 0
    offenders = []
    for path in REPO_ROOT.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        if path.suffix.lower() not in {".bib", ".tex", ".md", ".py", ".html", ".json", ".txt", ".yml", ".yaml", ".ipynb"}:
            continue
        if path.resolve() in SWEEP_EXEMPT:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        searched += 1
        for doi in WITHDRAWN_DOIS:
            if doi in text:
                offenders.append(f"{path.relative_to(REPO_ROOT)}: {doi}")
    assert searched > 20, f"only searched {searched} files, the sweep is not working"
    assert not offenders, "withdrawn DOIs are back: " + "; ".join(offenders)


def test_no_withdrawn_author_name_appears_in_either_bibliography():
    offenders = []
    for path in (REFERENCES, PAPER_BIB):
        text = path.read_text(encoding="utf-8")
        for name in WITHDRAWN_AUTHOR_NAMES:
            # The comment blocks in references.bib name what was withdrawn, so
            # only the author field counts.
            for line in text.splitlines():
                if line.lstrip().startswith("author") and name in line:
                    offenders.append(f"{path.name}: {name}")
    assert not offenders, "withdrawn author names are back: " + "; ".join(offenders)


def test_the_two_bibliographies_agree_on_every_entry_they_share():
    """paper/paper.bib held a different, also wrong, copy of three entries."""
    a, b = parse_bib(REFERENCES), parse_bib(PAPER_BIB)
    shared = set(a) & set(b)
    assert len(shared) == 9, f"expected 9 shared entries, found {len(shared)}"
    for key in sorted(shared):
        assert a[key] == b[key], (
            f"{key} differs between manuscript/references.bib and "
            f"paper/paper.bib; edit it in references.bib and copy it across"
        )
    only_in_paper = set(b) - set(a)
    assert not only_in_paper, (
        f"paper/paper.bib has entries references.bib does not: {sorted(only_in_paper)}"
    )


def test_every_cited_key_resolves_in_the_bibliography_that_document_declares():
    tex_keys = set()
    for group in re.findall(r"\\cite\{([^}]*)\}", PREPRINT_TEX.read_text(encoding="utf-8")):
        tex_keys.update(k.strip() for k in group.split(","))
    md_keys = set()
    for group in re.findall(r"\[@([^\]]+)\]", PAPER_MD.read_text(encoding="utf-8")):
        md_keys.update(k.strip().lstrip("@") for k in group.split(";"))

    assert tex_keys, "no \\cite keys found in the preprint source"
    assert md_keys, "no [@key] citations found in paper.md"

    missing_tex = tex_keys - set(parse_bib(REFERENCES))
    assert not missing_tex, (
        f"the preprint cites keys that references.bib does not define: {sorted(missing_tex)}"
    )
    missing_md = md_keys - set(parse_bib(PAPER_BIB))
    assert not missing_md, (
        f"paper.md cites keys that paper.bib does not define: {sorted(missing_md)}"
    )


def test_no_bibliography_entry_carries_a_withdrawn_doi():
    """Field-level check, so a withdrawn DOI cannot return as a citation even if
    the comment blocks are rewritten and the raw sweep above is loosened."""
    offenders = []
    for path in (REFERENCES, PAPER_BIB):
        for key, entry in parse_bib(path).items():
            if entry.get("doi") in WITHDRAWN_DOIS:
                offenders.append(f"{path.name}:{key}")
    assert not offenders, "withdrawn DOIs used as citations: " + "; ".join(offenders)
