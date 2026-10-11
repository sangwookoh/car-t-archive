#!/usr/bin/env python3
"""
CAR-T Publication Archive — PubMed collector.

Append-only by design: a PMID already present in the archive is never
rewritten, so past records stay byte-identical across runs.

Usage
-----
  python scripts/fetch_pubmed.py weekly
  python scripts/fetch_pubmed.py weekly --lookback 45
  python scripts/fetch_pubmed.py backfill --from 1987 --to 2014
  python scripts/fetch_pubmed.py backfill --from 2025 --to 2025

Environment
-----------
  NCBI_API_KEY   optional; raises the NCBI rate limit from 3/s to 10/s.
"""

import argparse
import datetime as dt
import json
import os
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
DELAY = 0.11 if API_KEY else 0.36          # stay under the NCBI rate limit
TOOL = "cart-publication-archive"
EMAIL = os.environ.get("NCBI_EMAIL", "").strip()

# NOTE on precision: do NOT use "CAR-T"[All Fields]. PubMed tokenises the
# hyphen, so any paper containing "car" and "t" matches — car-seat injury
# studies, fMRI papers and the like. Anchor on the MeSH descriptor plus exact
# Title/Abstract phrases instead.
# Bare "CAR-T"/"CAR-NK" matter: many abstracts write "CAR-T therapy" with no
# "cell", and requiring "cell" drops ~1,200 real papers.
# "CAR-M" is deliberately absent: it also denotes cancer-cell-membrane coatings
# in nanomedicine, so only unambiguous macrophage spellings are listed.
QUERY = (
    '("Receptors, Chimeric Antigen"[Mesh] '
    'OR "chimeric antigen receptor"[tiab] '
    'OR "chimeric antigen receptors"[tiab] '
    'OR "chimeric autoantibody receptor"[tiab] '
    'OR "chimeric autoantibody receptors"[tiab] '
    'OR "CAR-T"[tiab] OR "CAR T cell"[tiab] OR "CAR T cells"[tiab] '
    'OR "CAR-NK"[tiab] OR "CAR NK cell"[tiab] OR "CAR NK cells"[tiab] '
    'OR "CAAR-T"[tiab] OR "CAAR T cell"[tiab] OR "CAAR T cells"[tiab] '
    'OR "CAR-macrophage"[tiab] OR "CAR-macrophages"[tiab] OR "CAR-Ms"[tiab] '
    'OR "CAR-Treg"[tiab] OR "CAR-Tregs"[tiab] '
    'OR "CAR-NKT"[tiab] OR "CAR-iNKT"[tiab] '
    'OR "CAR-neutrophil"[tiab] OR "CAR-neutrophils"[tiab] '
    'OR "CAR-HSC"[tiab] OR "CAR gamma delta"[tiab] '
    # Pre-2010 terminology. "Chimeric antigen receptor" only became standard
    # around 2009-2010; the founding literature says T-body, chimeric
    # immunoreceptor or chimeric T cell receptor instead. Bare "chimeric
    # receptor" is deliberately excluded — it pulls in >1,200 unrelated
    # chimeric GPCR and cytokine-receptor papers.
    'OR "T-body"[tiab] OR "T-bodies"[tiab] '
    'OR "chimeric immunoreceptor"[tiab] OR "chimeric immunoreceptors"[tiab] '
    'OR "chimeric immune receptor"[tiab] OR "chimeric immune receptors"[tiab] '
    'OR "chimeric T cell receptor"[tiab] OR "chimeric T cell receptors"[tiab] '
    # Clinical papers often name only the product.
    'OR "tisagenlecleucel"[tiab] OR "axicabtagene"[tiab] '
    'OR "brexucabtagene"[tiab] OR "lisocabtagene"[tiab] '
    'OR "idecabtagene"[tiab] OR "ciltacabtagene"[tiab] '
    'OR "obecabtagene"[tiab] OR "equecabtagene"[tiab] '
    'OR "relmacabtagene"[tiab] OR "zevorcabtagene"[tiab] '
    'OR "anitocabtagene"[tiab] '
    'OR "axi-cel"[tiab] OR "tisa-cel"[tiab] OR "liso-cel"[tiab] '
    'OR "ide-cel"[tiab] OR "cilta-cel"[tiab] OR "brexu-cel"[tiab] '
    'OR "obe-cel"[tiab] '
    'OR "Kymriah"[tiab] OR "Yescarta"[tiab] OR "Tecartus"[tiab] '
    'OR "Breyanzi"[tiab] OR "Abecma"[tiab] OR "Carvykti"[tiab] '
    'OR "Aucatzyl"[tiab] '
    # Platform variants and target-prefixed spellings.
    'OR "UniCAR"[tiab] OR "switchable CAR"[tiab] OR "adapter CAR"[tiab] '
    'OR "universal chimeric antigen receptor"[tiab] '
    'OR "CAR construct"[tiab] OR "CAR constructs"[tiab] '
    'OR "CAR therapy"[tiab] OR "CAR therapies"[tiab] '
    'OR "CAR cell therapy"[tiab] '
    'OR "CD19 CAR"[tiab] OR "BCMA CAR"[tiab] '
    # --- Non-CAR synthetic receptor platforms -------------------------------
    # Same engineering premise (install a synthetic receptor so a cell becomes
    # the drug), different recognition logic. Kept in scope; tagged separately
    # so the CAR view stays clean.
    #
    # Tolerance platforms that present an autoantigen instead of an scFv:
    # "BAR-T" is overwhelmingly the Barrett's oesophagus cell line, so the bare
    # abbreviation is NOT searched. The spelled-out name is unambiguous, and
    # real BAR papers say "chimeric antigen/autoantibody receptor" too.
    'OR "B-cell antibody receptor"[tiab] OR "BAR Treg"[tiab] '
    'OR "BAR Tregs"[tiab] '
    # AIR: inflammation biosensor, not an antigen binder. PubMed has no phrase
    # index entry for "artificial immune receptor", so proximity is required —
    # a plain quoted search returns zero.
    'OR "artificial immune receptor"[tiab:~1] '
    'OR "AIR-Treg"[tiab] OR "AIR Treg"[tiab] OR "AIR Tregs"[tiab] '
    # Receptor architectures outside the CAR backbone:
    'OR "synNotch"[tiab] OR "synthetic Notch"[tiab] OR "SNIPR"[tiab] '
    'OR "TRuC"[tiab] OR "TRuCs"[tiab] '
    'OR "T cell receptor fusion construct"[tiab] '
    'OR "T cell receptor fusion protein"[tiab] '
    'OR "chimeric engulfment receptor"[tiab] '
    'OR "antibody-coupled T cell receptor"[tiab] OR "ACTR T cells"[tiab] '
    'OR "HIT receptor"[tiab] OR "synthetic immune receptor"[tiab:~1] '
    # CAR-adjacent control modules. Bare "iCAR" is omitted — it also denotes
    # unrelated constructs; the spelled-out form is unambiguous.
    'OR "inhibitory chimeric antigen receptor"[tiab] '
    'OR "chimeric switch receptor"[tiab] '
    'OR "chimeric costimulatory receptor"[tiab]) '
    'NOT (Editorial[pt] OR Comment[pt] OR Letter[pt] OR "Case Reports"[pt] '
    'OR Retracted Publication[pt])'
)

CATEGORY_RULES = [
    ("CAAR-T",     [r"chimeric autoantibody receptor", r"\bCAAR[-\s]?T\b", r"\bCAAR\b"]),
    ("CAR-NK",     [r"CAR[-\s]?NK", r"chimeric antigen receptor[^.]{0,40}natural killer"]),
    ("CAR-M",      [r"CAR[-\s]macrophage", r"\bCAR-Ms\b",
                    r"chimeric antigen receptor[^.]{0,40}macrophage",
                    r"\bCAR-M\b[^.]{0,60}macrophage", r"macrophage[^.]{0,60}\bCAR-M\b"]),
    ("CAR-Treg",   [r"CAR[-\s]?Tregs?\b", r"chimeric antigen receptor[^.]{0,40}regulatory T"]),
    ("CAR-NKT",    [r"CAR[-\s]?i?NKT", r"chimeric antigen receptor[^.]{0,40}NKT"]),
    ("CAR-γδT",    [r"CAR[-\s]?gamma[-\s]?delta", r"CAR[-\s]?γδ",
                    r"chimeric antigen receptor[^.]{0,40}(gamma delta|γδ)"]),
    ("CAR-neutrophil", [r"CAR[-\s]?neutrophil"]),
    # Includes pre-2010 names and spellings that never write out "CAR-T", so
    # the founding literature is findable under the modern term.
    ("CAR-T",      [r"CAR[-\s]?T\b", r"chimeric antigen receptor",
                    r"\bT-bod(y|ies)\b", r"chimeric immunoreceptor",
                    r"chimeric immune receptor", r"chimeric T[-\s]cell receptor",
                    r"\bUniCAR\b", r"\bCAR (construct|therap|platform|cell)",
                    r"\b(CD19|BCMA|CD22|GD2|HER2|CD7|CD20|CD30)[-\s]CAR\b",
                    r"tisagenlecleucel", r"cabtagene", r"\b(axi|tisa|liso|ide|cilta|brexu|obe)-cel\b",
                    r"Kymriah|Yescarta|Tecartus|Breyanzi|Abecma|Carvykti|Aucatzyl"]),
    ("in vivo",    [r"\bin vivo\b[^.]{0,50}\bCAR\b", r"lipid nanoparticle", r"\bLNP\b"]),
    ("Autoimmune", [r"autoimmun", r"\bSLE\b", r"lupus", r"myasthenia", r"pemphigus",
                    r"multiple sclerosis", r"myositis", r"scleroderma", r"arthritis",
                    r"NMOSD", r"neuromyelitis", r"type 1 diabetes"]),
    ("SLE",        [r"\bSLE\b", r"lupus erythematosus"]),
    # "acetylcholine receptor" alone is not enough: fetal AChR is a CAR target
    # in rhabdomyosarcoma, which has nothing to do with myasthenia gravis.
    ("Myasthenia", [r"myasthenia", r"\bMuSK\b", r"\bAChR\b[^.]{0,60}myasthen",
                    r"myasthen[^.]{0,60}\bAChR\b"]),
    ("Pemphigus",  [r"pemphigus", r"desmoglein", r"\bDsg[13]\b"]),
    ("NMOSD",      [r"neuromyelitis", r"NMOSD", r"aquaporin[-\s]?4", r"\bAQP4\b"]),
    ("Solid tumor",[r"solid tumou?r", r"glioblastoma", r"\bGBM\b", r"pancreatic cancer",
                    r"hepatocellular", r"ovarian cancer", r"sarcoma"]),
    ("Hematologic",[r"lymphoma", r"leukemi|leukaemi", r"myeloma", r"\bB-ALL\b", r"\bT-ALL\b"]),
    ("Manufacturing", [r"manufactur", r"\bGMP\b", r"vector production", r"transduction efficiency",
                       r"scale[-\s]?up", r"cryopreserv"]),
    ("Toxicity",   [r"cytokine release syndrome", r"\bCRS\b", r"\bICANS\b", r"neurotoxicity",
                    r"on-target off-tumor", r"adverse event"]),
    ("Approved product", [r"tisagenlecleucel", r"cabtagene", r"\baxi-cel\b", r"\btisa-cel\b",
                          r"\bliso-cel\b", r"\bide-cel\b", r"\bcilta-cel\b", r"\bbrexu-cel\b",
                          r"\bobe-cel\b", r"Kymriah", r"Yescarta", r"Tecartus", r"Breyanzi",
                          r"Abecma", r"Carvykti", r"Aucatzyl"]),
    ("Early era",  [r"\bT-bod(y|ies)\b", r"chimeric immunoreceptor",
                    r"chimeric immune receptor", r"chimeric T[-\s]cell receptor"]),
    # --- non-CAR synthetic receptor platforms ---
    ("BAR",        [r"B[-\s]cell antibody receptor", r"\bBAR[-\s]Tregs?\b",
                    r"\bBAR[-\s]T cells?\b"]),
    ("AIR",        [r"artificial immune receptor", r"\bAIR[-\s]Tregs?\b",
                    r"\bCD40[-\s]AIR\b"]),
    ("synNotch",   [r"syn[-\s]?Notch", r"synthetic Notch", r"\bSNIPR\b"]),
    ("TRuC/TFP",   [r"\bTRuCs?\b", r"T cell receptor fusion (construct|protein)",
                    r"\bTCR fusion\b"]),
    ("CER",        [r"chimeric engulfment receptor", r"\bCER[-\s]?T\b"]),
    ("iCAR/switch",[r"inhibitory chimeric antigen receptor", r"\biCARs?\b",
                    r"chimeric switch receptor", r"chimeric costimulatory receptor"]),
]
# ------------------------------------------------- cell/gene therapy gate
# A keyword can be right and the paper still wrong: "BAR-T" is usually the
# Barrett's oesophagus cell line, "CAR-M" a cancer-cell-membrane nanoparticle,
# "AIR" an airway abbreviation. Papers whose ONLY hook is such an abbreviation
# must show, in the abstract, that an engineered receptor is being put into
# cells. Papers matching an unambiguous term skip the gate entirely.

# Terms specific enough that no corroboration is needed.
UNAMBIGUOUS = [
    r"chimeric antigen receptor", r"chimeric autoantibody receptor",
    # MeSH inverts the word order, and a record with no abstract may carry
    # nothing else. Without this, indexed papers get thrown away.
    r"Receptors, Chimeric Antigen", r"Immunotherapy, Adoptive",
    r"\bCAR[-\s]?T\b", r"\bCAR[-\s]?NK\b", r"\bCAAR[-\s]?T\b",
    r"CAR[-\s]macrophage", r"CAR[-\s]?Tregs?\b", r"CAR[-\s]?i?NKT\b",
    r"chimeric immunoreceptor", r"chimeric immune receptor",
    r"chimeric T[-\s]cell receptor", r"\bT-bod(y|ies)\b",
    r"B[-\s]cell antibody receptor", r"artificial immune receptor",
    r"syn[-\s]?Notch", r"\bSNIPR\b", r"chimeric engulfment receptor",
    r"T cell receptor fusion (construct|protein)",
    r"tisagenlecleucel", r"cabtagene", r"\b(axi|tisa|liso|ide|cilta|brexu|obe)-cel\b",
    r"Kymriah|Yescarta|Tecartus|Breyanzi|Abecma|Carvykti|Aucatzyl",
]

# Evidence that an engineered receptor is being introduced into cells.
EVIDENCE = [
    # receptor architecture
    r"\bscFv\b", r"single[-\s]chain variable", r"\bCD3\s?[ζz]\b", r"CD3[-\s]?zeta",
    r"costimulatory domain", r"co[-\s]stimulatory domain", r"\b4-1BB\b", r"\bCD28\b",
    r"hinge", r"transmembrane domain", r"ectodomain", r"signaling domain",
    r"binding domain", r"\bbinder\b",
    # gene delivery
    r"lentivir", r"retrovir", r"transduc", r"transfect", r"electroporat",
    r"\bCRISPR\b", r"gene[-\s]edit", r"\bvector\b", r"piggyBac", r"sleeping beauty",
    r"mRNA[-\s](delivery|encoding)", r"lipid nanoparticle",
    # cell therapy practice
    r"adoptive (transfer|cell|immunotherapy)", r"cell therapy", r"cellular therapy",
    r"engineered (T|NK|immune|Treg)", r"infusion", r"leukapheres",
    r"lymphodeplet", r"cytotoxicity assay", r"tumor (cell )?(killing|lysis)",
    r"effector (cell|function)", r"immunotherap",
]

# Strong signals the match is something else entirely.
EXCLUDE = [
    r"Barrett", r"oesophag|esophag", r"bile acid receptor", r"\bFXR\b", r"\bTGR5\b",
    r"cell membrane[-\s]coat", r"membrane[-\s]coated nanoparticle",
    r"airway (inflammation|resistance|hyperrespons)", r"\bCAR[-\s]?M\b.{0,30}nano",
]


def is_cell_gene_therapy(blob):
    """True when the abstract supports reading this as cell/gene therapy.

    Returns (verdict, reason) so the audit command can explain itself.
    """
    if any(re.search(p, blob, re.I) for p in UNAMBIGUOUS):
        return True, "unambiguous term"
    hits = [p for p in EVIDENCE if re.search(p, blob, re.I)]
    bad = [p for p in EXCLUDE if re.search(p, blob, re.I)]
    if len(hits) >= 2 and not bad:
        return True, f"{len(hits)} evidence signals"
    if bad:
        return False, "excluded context"
    return False, f"only {len(hits)} evidence signal(s)"


# --------------------------------------------------------------- disease axis
# Deliberately finer than the broad category tags: the point is to be able to
# ask "what is being tried in myasthenia gravis" and get an answer, which a
# single "Autoimmune" bucket cannot give.
DISEASE_RULES = [
    # --- haematologic malignancy ---
    ("B-ALL",            [r"\bB[-\s]?ALL\b", r"B[-\s]cell acute lymphoblastic",
                          r"acute lymphoblastic leukemi|leukaemi"]),
    ("T-ALL",            [r"\bT[-\s]?ALL\b", r"T[-\s]cell acute lymphoblastic"]),
    ("AML",              [r"\bAML\b", r"acute myeloid leukemi|leukaemi"]),
    ("CLL",              [r"\bCLL\b", r"chronic lymphocytic leukemi|leukaemi"]),
    ("DLBCL/NHL",        [r"\bDLBCL\b", r"diffuse large B[-\s]cell lymphoma",
                          r"non[-\s]Hodgkin", r"follicular lymphoma",
                          r"mantle cell lymphoma", r"\bMCL\b"]),
    ("Hodgkin lymphoma", [r"Hodgkin lymphoma", r"Hodgkin's disease"]),
    ("Multiple myeloma", [r"multiple myeloma", r"\bMM\b(?![-\s]?\d)", r"plasma cell myeloma"]),
    ("T-cell lymphoma",  [r"T[-\s]cell lymphoma", r"cutaneous T[-\s]cell",
                          r"peripheral T[-\s]cell lymphoma"]),
    # --- solid tumour ---
    ("Glioma/GBM",       [r"glioblastoma", r"\bGBM\b", r"glioma", r"\bDIPG\b",
                          r"diffuse midline glioma"]),
    ("Neuroblastoma",    [r"neuroblastoma"]),
    ("Sarcoma",          [r"sarcoma", r"osteosarcoma", r"Ewing"]),
    ("Lung cancer",      [r"lung cancer", r"\bNSCLC\b", r"\bSCLC\b",
                          r"lung adenocarcinoma"]),
    ("Breast cancer",    [r"breast cancer", r"breast carcinoma", r"\bTNBC\b"]),
    ("Pancreatic cancer",[r"pancreatic cancer", r"pancreatic ductal adenocarcinoma", r"\bPDAC\b"]),
    ("Liver cancer",     [r"hepatocellular carcinoma", r"\bHCC\b", r"liver cancer"]),
    ("Gastric/CRC",      [r"gastric cancer", r"colorectal cancer", r"\bCRC\b",
                          r"colon cancer", r"gastroesophageal"]),
    ("Ovarian cancer",   [r"ovarian cancer", r"ovarian carcinoma"]),
    ("Prostate cancer",  [r"prostate cancer", r"prostate carcinoma"]),
    ("Renal cancer",     [r"renal cell carcinoma", r"\bRCC\b", r"kidney cancer"]),
    ("Melanoma",         [r"melanoma"]),
    ("Mesothelioma",     [r"mesothelioma"]),
    ("Head and neck",    [r"head and neck (squamous|cancer|carcinoma)", r"\bHNSCC\b"]),
    # --- autoimmune / immune-mediated ---
    ("SLE",              [r"\bSLE\b", r"systemic lupus erythematosus", r"lupus nephritis"]),
    ("Myasthenia gravis",[r"myasthenia gravis", r"\bgMG\b",
                          r"\bMuSK\b[^.]{0,60}myasthen", r"myasthen[^.]{0,60}\bMuSK\b"]),
    ("Pemphigus",        [r"pemphigus", r"\bDsg[13]\b", r"desmoglein"]),
    ("NMOSD",            [r"neuromyelitis optica", r"\bNMOSD?\b", r"aquaporin[-\s]?4", r"\bAQP4\b"]),
    ("Multiple sclerosis",[r"multiple sclerosis", r"\bMS\b(?=[^.]{0,40}(patient|disease|treat))"]),
    ("Rheumatoid arthritis",[r"rheumatoid arthritis", r"\bRA\b(?=[^.]{0,40}patient)"]),
    ("Systemic sclerosis",[r"systemic sclerosis", r"scleroderma", r"\bSSc\b"]),
    ("Myositis",         [r"myositis", r"dermatomyositis", r"\bIIM\b"]),
    ("Type 1 diabetes",  [r"type 1 diabetes", r"\bT1D\b", r"autoimmune diabetes"]),
    ("IBD",              [r"inflammatory bowel disease", r"Crohn", r"ulcerative colitis"]),
    ("Sjogren",          [r"Sj[öo]gren"]),
    ("ANCA vasculitis",  [r"\bANCA\b", r"vasculitis"]),
    ("ITP/AIHA",         [r"immune thrombocytopeni", r"\bITP\b",
                          r"autoimmune hemolytic|haemolytic anemia|anaemia"]),
    ("Hemophilia inhibitor",[r"h(a)?emophilia", r"\bFVIII\b", r"factor VIII inhibitor"]),
    # --- transplantation and other ---
    ("GvHD/transplant",  [r"graft[-\s]versus[-\s]host", r"\bGvHD\b",
                          r"transplant rejection", r"allograft tolerance"]),
    ("Infection",        [r"\bHIV\b", r"hepatitis B virus", r"\bHBV\b",
                          r"fungal infection", r"\bCMV\b"]),
    ("Fibrosis",         [r"cardiac fibrosis", r"\bFAP\b[^.]{0,40}fibro", r"liver fibrosis"]),
    ("Senescence",       [r"senescen", r"senolytic", r"\buPAR\b"]),
]

# --------------------------------------------------------------- target axis
TARGET_RULES = [
    (t, [rf"\b{re.escape(t)}\b"]) for t in
    ["CD19", "CD20", "CD22", "CD30", "CD33", "CD38", "CD123", "CD7", "CD5",
     "CD70", "BCMA", "GPRC5D", "SLAMF7", "FCRL5", "GD2", "HER2", "EGFRvIII",
     "EGFR", "IL13Ra2", "B7-H3", "MUC1", "MUC16", "PSMA", "ROR1", "CEA",
     "NKG2D", "FAP", "uPAR", "GPC3", "Claudin18.2", "CLDN18.2", "AQP4",
     "MuSK", "Dsg3", "Dsg1", "PLA2R", "FVIII", "CD40L", "TL1A"]
] + [
    ("Mesothelin", [r"mesothelin"]),
    ("AChR",       [r"\bAChR\b", r"acetylcholine receptor"]),
    ("CD40",       [r"\bCD40\b(?![Ll])"]),
]

REVIEW_TYPES = {"Review", "Systematic Review", "Meta-Analysis"}
TRIAL_TYPES = {"Clinical Trial", "Randomized Controlled Trial", "Clinical Trial, Phase I",
               "Clinical Trial, Phase II", "Clinical Trial, Phase III", "Observational Study"}


# ----------------------------------------------------------------- utilities
def _url(endpoint, **params):
    params.setdefault("db", "pubmed")
    params["tool"] = TOOL
    if EMAIL:
        params["email"] = EMAIL
    if API_KEY:
        params["api_key"] = API_KEY
    return f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"


def _get(url, retries=5):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": TOOL})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as exc:                      # noqa: BLE001
            last = exc
            wait = min(60, 2 ** attempt * 2)
            print(f"  retry {attempt + 1}/{retries} in {wait}s — {exc}", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"request failed after {retries} attempts: {last}")


def text_of(node):
    """Flatten an element's text, including nested markup such as <i>/<sup>."""
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip() if node is not None else ""


def classify(blob, rules):
    out = []
    for name, pats in rules:
        if any(re.search(p, blob, re.I) for p in pats):
            out.append(name)
    # A CAR-NK or CAAR-T paper says "chimeric antigen receptor" too, which
    # would otherwise tag it CAR-T. Keep CAR-T only when the T-cell form is
    # named explicitly.
    others = {"CAR-NK", "CAAR-T", "CAR-M", "CAR-Treg", "CAR-NKT",
              "CAR-γδT", "CAR-neutrophil"}
    if "CAR-T" in out and others & set(out):
        if not re.search(r"CAR[-\s]?T\b|chimeric antigen receptor[^.]{0,40}T cell", blob, re.I):
            out.remove("CAR-T")
    # Umbrella tag: lets the reader separate CAR work from everything else in
    # one click, without having to know each platform's acronym.
    NON_CAR = {"BAR", "AIR", "synNotch", "TRuC/TFP", "CER"}
    if NON_CAR & set(out):
        out.append("Non-CAR platform")
    return out


def iso_week(y, m, d):
    try:
        return "%d-W%02d" % dt.date(int(y), int(m or 1), int(d or 1)).isocalendar()[:2]
    except Exception:                                  # noqa: BLE001
        return None


# ----------------------------------------------------------------- esearch
def esearch(term, mindate=None, maxdate=None):
    """Return every PMID for a query, paging through the History server."""
    params = dict(term=term, usehistory="y", retmax=0, retmode="json")
    if mindate:
        params.update(datetype="pdat", mindate=mindate, maxdate=maxdate or mindate)
    meta = json.loads(_get(_url("esearch.fcgi", **params)))["esearchresult"]
    total = int(meta["count"])
    print(f"  esearch: {total} records")
    if total == 0:
        return []
    # NCBI caps retstart at 9,999 for esearch. Anything larger has to be split
    # into narrower date windows by the caller.
    if total > 9999:
        raise RuntimeError(
            f"esearch returned {total} records; NCBI only serves the first 9,999. "
            "Split the query into smaller date windows.")
    pmids, step = [], 9000
    for start in range(0, total, step):
        page = dict(term=term, WebEnv=meta["webenv"], query_key=meta["querykey"],
                    retstart=start, retmax=step, retmode="json")
        if mindate:
            page.update(datetype="pdat", mindate=mindate, maxdate=maxdate or mindate)
        pmids += json.loads(_get(_url("esearch.fcgi", **page)))["esearchresult"]["idlist"]
        time.sleep(DELAY)
    return pmids


# ----------------------------------------------------------------- efetch
def parse_article(art):
    cit = art.find("MedlineCitation")
    if cit is None:
        return None
    pmid = text_of(cit.find("PMID"))
    if not pmid:
        return None
    a = cit.find("Article")
    title = text_of(a.find("ArticleTitle")).rstrip(".") if a is not None else ""

    parts = []
    for ab in a.findall("./Abstract/AbstractText") if a is not None else []:
        label = ab.get("Label")
        body = text_of(ab)
        parts.append(f"{label}: {body}" if label else body)
    abstract = " ".join(p for p in parts if p)

    authors, full = [], []
    for au in a.findall("./AuthorList/Author") if a is not None else []:
        last = text_of(au.find("LastName"))
        if not last:
            coll = text_of(au.find("CollectiveName"))
            if coll:
                authors.append(coll)
                full.append(coll)
            continue
        ini = text_of(au.find("Initials"))
        fore = text_of(au.find("ForeName"))
        authors.append(f"{last} {ini}".strip())
        full.append(f"{fore} {last}".strip())

    j = a.find("Journal") if a is not None else None
    journal = text_of(j.find("ISOAbbreviation")) if j is not None else ""
    journal_full = text_of(j.find("Title")) if j is not None else ""

    # publication date: prefer ArticleDate (electronic), fall back to JournalIssue
    y = m = d = ""
    ad = a.find("ArticleDate") if a is not None else None
    if ad is not None:
        y, m, d = text_of(ad.find("Year")), text_of(ad.find("Month")), text_of(ad.find("Day"))
    if not y and j is not None:
        pd = j.find("./JournalIssue/PubDate")
        if pd is not None:
            y = text_of(pd.find("Year"))
            m = text_of(pd.find("Month"))
            d = text_of(pd.find("Day"))
            if not y:
                mo = re.search(r"\b(19|20)\d{2}\b", text_of(pd.find("MedlineDate")))
                y = mo.group(0) if mo else ""
    months = {"jan": "01", "feb": "02", "mar": "03", "apr": "04", "may": "05", "jun": "06",
              "jul": "07", "aug": "08", "sep": "09", "oct": "10", "nov": "11", "dec": "12"}
    if m and not m.isdigit():
        m = months.get(m[:3].lower(), "")
    m = m.zfill(2) if m else ""
    d = d.zfill(2) if d else ""

    doi = ""
    for aid in art.findall("./PubmedData/ArticleIdList/ArticleId"):
        if aid.get("IdType") == "doi":
            doi = text_of(aid)
    if not doi and a is not None:
        for eid in a.findall("ELocationID"):
            if eid.get("EIdType") == "doi":
                doi = text_of(eid)

    kw = [text_of(k) for k in cit.findall("./KeywordList/Keyword")]
    mesh = [text_of(h.find("DescriptorName")) for h in cit.findall("./MeshHeadingList/MeshHeading")]
    ptypes = [text_of(p) for p in a.findall("./PublicationTypeList/PublicationType")] if a is not None else []

    blob = " ".join([title, abstract, " ".join(kw), " ".join(mesh)])

    # Gate before anything else: a keyword hit is not a reason to keep a paper
    # about a cell line that happens to share an abbreviation.
    ok, why = is_cell_gene_therapy(blob)
    if not ok:
        return {"_rejected": True, "pmid": pmid, "title": title, "reason": why}

    if set(ptypes) & REVIEW_TYPES:
        pub_type = "Review"
    elif set(ptypes) & TRIAL_TYPES:
        pub_type = "Clinical"
    else:
        pub_type = "Original"

    return {
        "pmid": pmid,
        "doi": doi or None,
        "title": title,
        "abstract": abstract,
        "authors": authors,
        "authors_full": full,
        "first_author": authors[0] if authors else None,
        "last_author": authors[-1] if len(authors) > 1 else None,
        "journal": journal or journal_full,
        "journal_full": journal_full or journal,
        "year": int(y) if y.isdigit() else None,
        "date": "-".join(x for x in (y, m, d) if x),
        "week": iso_week(y, m, d),
        "keywords": kw,
        "mesh": mesh,
        "category": classify(blob, CATEGORY_RULES),
        "disease": classify(blob, DISEASE_RULES),
        "target": classify(blob, TARGET_RULES),
        "pub_type": pub_type,
        "article_types": ptypes,
        "added_at": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def efetch(pmids, batch=200, show_rejects=8):
    out, rejected = [], []
    for i in range(0, len(pmids), batch):
        chunk = pmids[i:i + batch]
        xml = _get(_url("efetch.fcgi", id=",".join(chunk), retmode="xml"))
        try:
            root = ET.fromstring(xml)
        except ET.ParseError as exc:
            print(f"  ! XML parse error on batch {i // batch}: {exc}", file=sys.stderr)
            continue
        for art in root.findall("PubmedArticle"):
            rec = parse_article(art)
            if not rec:
                continue
            if rec.get("_rejected"):
                rejected.append(rec)
            else:
                out.append(rec)
        print(f"  efetch {min(i + batch, len(pmids))}/{len(pmids)}")
        time.sleep(DELAY)
    if rejected:
        print(f"  gate rejected {len(rejected)} not-cell/gene-therapy papers; e.g.")
        for r in rejected[:show_rejects]:
            print(f"    - {r['pmid']} ({r['reason']}): {r['title'][:72]}")
    return out


# ----------------------------------------------------------------- storage
def shard_path(kind, year):
    return DATA / f"{kind}-{year}.json"


def load_json(path, default):
    if path.exists():
        with path.open(encoding="utf-8") as f:
            return json.load(f)
    return default


def save_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, separators=(",", ":"))


def known_pmids():
    """Every PMID already stored, across all year shards."""
    seen = set()
    for p in DATA.glob("idx-*.json"):
        for r in load_json(p, {}).get("records", []):
            seen.add(r["pmid"])
    return seen


def slim(rec):
    """Index record: everything the list view and search need except the abstract."""
    return {
        "pmid": rec["pmid"], "doi": rec["doi"], "title": rec["title"],
        "authors": rec["authors"], "first_author": rec["first_author"],
        "last_author": rec["last_author"], "journal": rec["journal"],
        "year": rec["year"], "date": rec["date"], "week": rec["week"],
        "category": rec["category"], "disease": rec.get("disease", []),
        "target": rec.get("target", []), "pub_type": rec["pub_type"],
        "keywords": rec["keywords"][:12], "mesh": rec["mesh"][:14],
    }


def store(records):
    """Append new records into their year shards. Existing PMIDs are left untouched."""
    by_year = {}
    for r in records:
        by_year.setdefault(r["year"] or 0, []).append(r)

    added_total = 0
    for year, recs in sorted(by_year.items()):
        ip, ap = shard_path("idx", year), shard_path("abs", year)
        idx = load_json(ip, {"year": year, "records": []})
        abs_ = load_json(ap, {})
        have = {r["pmid"] for r in idx["records"]}
        added = 0
        for r in recs:
            if r["pmid"] in have:
                continue
            idx["records"].append(slim(r))
            if r["abstract"]:
                abs_[r["pmid"]] = r["abstract"]
            have.add(r["pmid"])
            added += 1
        if added:
            idx["records"].sort(key=lambda r: (r["date"] or "", r["pmid"]), reverse=True)
            idx["count"] = len(idx["records"])
            idx["updated_at"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            save_json(ip, idx)
            save_json(ap, abs_)
            added_total += added
        print(f"  {year}: +{added} (total {len(idx['records'])})")
    return added_total


def write_manifest():
    years, total = [], 0
    for p in sorted(DATA.glob("idx-*.json")):
        d = load_json(p, {})
        n = len(d.get("records", []))
        years.append({"year": d.get("year"), "count": n})
        total += n
    years.sort(key=lambda y: y["year"] or 0, reverse=True)
    save_json(DATA / "manifest.json", {
        "schema_version": 2,
        "source": "PubMed (NCBI E-utilities)",
        "query": QUERY,
        "total": total,
        "years": years,
        "generated_at": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    })
    print(f"manifest: {total} records across {len(years)} years")


# ----------------------------------------------------------------- commands
def run(pmids, label):
    print(f"[{label}] candidates: {len(pmids)}")
    fresh = [p for p in pmids if p not in known_pmids()]
    print(f"[{label}] new: {len(fresh)}")
    if not fresh:
        write_manifest()
        return 0
    added = store(efetch(fresh))
    write_manifest()
    print(f"[{label}] added {added}")
    return added


def cmd_weekly(args):
    # A 45-day window, not 7: PubMed indexes some papers late, and the PMID
    # dedup above means re-seeing a paper costs nothing. A 7-day window drops
    # anything indexed after its publication week, permanently.
    end = dt.date.today()
    start = end - dt.timedelta(days=args.lookback)
    print(f"weekly window {start} .. {end}")
    run(esearch(QUERY, start.strftime("%Y/%m/%d"), end.strftime("%Y/%m/%d")), "weekly")


def cmd_prune(args):
    """Remove stored records that the current QUERY no longer matches.

    Append-only keeps history stable, but it also preserves mistakes: records
    collected under an earlier, looser query stay forever. This rebuilds the
    valid PMID set from one esearch pass (cheap — no efetch) and drops
    anything outside it. Korean summaries are never touched.
    """
    # Year by year, because one esearch can only serve 9,999 records.
    stored_years = sorted(
        y for y in (load_json(p, {}).get("year") for p in DATA.glob("idx-*.json")) if y
    )
    if not stored_years:
        print("nothing stored; nothing to prune")
        return
    lo, hi = stored_years[0] - 1, stored_years[-1] + 1
    valid = set()
    for year in range(lo, hi + 1):
        print(f"--- esearch {year} ---")
        valid |= set(esearch(QUERY, f"{year}/01/01", f"{year}/12/31"))
    print(f"valid PMIDs under current query: {len(valid)}")
    if len(valid) < 1000:
        print("refusing to prune: suspiciously small result set", file=sys.stderr)
        sys.exit(1)

    removed_total = 0
    for ip in sorted(DATA.glob("idx-*.json")):
        idx = load_json(ip, {"records": []})
        year = idx.get("year")
        abs_map = load_json(shard_path("abs", year), {})

        def survives(r):
            if r["pmid"] not in valid:
                return False, "outside query"
            blob = " ".join([r.get("title") or "", abs_map.get(r["pmid"], ""),
                             " ".join(r.get("keywords") or []),
                             " ".join(r.get("mesh") or [])])
            ok, why = is_cell_gene_therapy(blob)
            return ok, why

        verdicts = {r["pmid"]: survives(r) for r in idx["records"]}
        keep = [r for r in idx["records"] if verdicts[r["pmid"]][0]]
        dropped = [r["pmid"] for r in idx["records"] if not verdicts[r["pmid"]][0]]
        if not dropped:
            continue
        for r in idx["records"]:
            if not verdicts[r["pmid"]][0]:
                print(f"  drop {year} {r['pmid']} ({verdicts[r['pmid']][1]}): {r['title'][:66]}")
        idx["records"] = keep
        idx["count"] = len(keep)
        idx["updated_at"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        if keep:
            save_json(ip, idx)
        else:
            ip.unlink()
        ap = shard_path("abs", year)
        if ap.exists():
            for pmid in dropped:
                abs_map.pop(pmid, None)
            save_json(ap, abs_map) if abs_map else ap.unlink()
        removed_total += len(dropped)
    write_manifest()
    print(f"pruned {removed_total} records")


def cmd_retag(args):
    """Recompute category and pub_type for stored records, offline.

    Classification rules change as the field's vocabulary changes; the
    bibliographic facts do not. This rewrites only the derived fields, leaving
    titles, authors, dates and added_at untouched, and needs no PubMed call.
    """
    changed_files = changed_records = 0
    for ip in sorted(DATA.glob("idx-*.json")):
        idx = load_json(ip, {"records": []})
        year = idx.get("year")
        abs_ = load_json(shard_path("abs", year), {})
        dirty = False
        for r in idx["records"]:
            blob = " ".join([r.get("title") or "", abs_.get(r["pmid"], ""),
                             " ".join(r.get("keywords") or []),
                             " ".join(r.get("mesh") or [])])
            new = {
                "category": classify(blob, CATEGORY_RULES),
                "disease": classify(blob, DISEASE_RULES),
                "target": classify(blob, TARGET_RULES),
            }
            if any(new[k] != r.get(k) for k in new):
                r.update(new)
                dirty = True
                changed_records += 1
        if dirty:
            idx["updated_at"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            save_json(ip, idx)
            changed_files += 1
            print(f"  {year}: retagged")
    write_manifest()
    print(f"retagged {changed_records} records across {changed_files} years")


def cmd_audit(args):
    """Dry run of the cell/gene-therapy gate over stored records.

    Changes nothing. Use it to see what prune would remove, and to sanity-check
    the gate before trusting it with the archive.
    """
    total = flagged = 0
    reasons, samples = {}, []
    for ip in sorted(DATA.glob("idx-*.json")):
        idx = load_json(ip, {"records": []})
        year = idx.get("year")
        abs_map = load_json(shard_path("abs", year), {})
        for r in idx["records"]:
            total += 1
            blob = " ".join([r.get("title") or "", abs_map.get(r["pmid"], ""),
                             " ".join(r.get("keywords") or []),
                             " ".join(r.get("mesh") or [])])
            ok, why = is_cell_gene_therapy(blob)
            if ok:
                continue
            flagged += 1
            reasons[why] = reasons.get(why, 0) + 1
            if len(samples) < 40:
                has_abs = "abstract" if abs_map.get(r["pmid"]) else "NO ABSTRACT"
                samples.append(f"  {year} {r['pmid']} [{why}; {has_abs}] {r['title'][:66]}")
    print(f"audited {total} records; gate would reject {flagged} "
          f"({100 * flagged / max(total, 1):.1f}%)")
    for why, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"  {n:6d}  {why}")
    print("\nsamples:")
    print("\n".join(samples))


def cmd_backfill(args):
    for year in range(args.from_year, args.to_year + 1):
        print(f"--- {year} ---")
        run(esearch(QUERY, f"{year}/01/01", f"{year}/12/31"), str(year))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("weekly");   w.add_argument("--lookback", type=int, default=45)
    b = sub.add_parser("backfill")
    b.add_argument("--from", dest="from_year", type=int, required=True)
    b.add_argument("--to", dest="to_year", type=int, required=True)
    sub.add_parser("prune")
    sub.add_parser("retag")
    sub.add_parser("audit")
    args = ap.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)
    {"weekly": cmd_weekly, "backfill": cmd_backfill, "prune": cmd_prune,
     "retag": cmd_retag, "audit": cmd_audit}[args.cmd](args)


if __name__ == "__main__":
    main()
