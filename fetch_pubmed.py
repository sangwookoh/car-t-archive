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
QUERY = (
    '("Receptors, Chimeric Antigen"[Mesh] '
    'OR "chimeric antigen receptor"[tiab] '
    'OR "chimeric antigen receptors"[tiab] '
    'OR "chimeric autoantibody receptor"[tiab] '
    'OR "chimeric autoantibody receptors"[tiab] '
    'OR "CAR-T cell"[tiab] OR "CAR-T cells"[tiab] '
    'OR "CAR T cell"[tiab] OR "CAR T cells"[tiab] '
    'OR "CAR-NK cell"[tiab] OR "CAR-NK cells"[tiab] '
    'OR "CAR NK cell"[tiab] OR "CAR NK cells"[tiab] '
    'OR "CAAR-T"[tiab] OR "CAAR T cell"[tiab]) '
    'NOT (Editorial[pt] OR Comment[pt] OR Letter[pt] OR "Case Reports"[pt] '
    'OR Retracted Publication[pt])'
)

CATEGORY_RULES = [
    ("CAAR-T",     [r"chimeric autoantibody receptor", r"\bCAAR[-\s]?T\b", r"\bCAAR\b"]),
    ("CAR-NK",     [r"CAR[-\s]?NK", r"chimeric antigen receptor[^.]{0,40}natural killer"]),
    ("CAR-T",      [r"CAR[-\s]?T\b", r"chimeric antigen receptor"]),
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
    if "CAR-T" in out and ("CAR-NK" in out or "CAAR-T" in out):
        if not re.search(r"CAR[-\s]?T\b|chimeric antigen receptor[^.]{0,40}T cell", blob, re.I):
            out.remove("CAR-T")
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
    pmids, step = [], 10000
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
        "pub_type": pub_type,
        "article_types": ptypes,
        "added_at": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def efetch(pmids, batch=200):
    out = []
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
            if rec:
                out.append(rec)
        print(f"  efetch {min(i + batch, len(pmids))}/{len(pmids)}")
        time.sleep(DELAY)
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
        "category": rec["category"], "pub_type": rec["pub_type"],
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
    valid = set(esearch(QUERY))
    print(f"valid PMIDs under current query: {len(valid)}")
    if len(valid) < 1000:
        print("refusing to prune: suspiciously small result set", file=sys.stderr)
        sys.exit(1)

    removed_total = 0
    for ip in sorted(DATA.glob("idx-*.json")):
        idx = load_json(ip, {"records": []})
        keep = [r for r in idx["records"] if r["pmid"] in valid]
        dropped = [r["pmid"] for r in idx["records"] if r["pmid"] not in valid]
        if not dropped:
            continue
        year = idx.get("year")
        for r in idx["records"]:
            if r["pmid"] not in valid:
                print(f"  drop {year} {r['pmid']}: {r['title'][:70]}")
        idx["records"] = keep
        idx["count"] = len(keep)
        idx["updated_at"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        if keep:
            save_json(ip, idx)
        else:
            ip.unlink()
        ap = shard_path("abs", year)
        if ap.exists():
            abs_ = load_json(ap, {})
            for pmid in dropped:
                abs_.pop(pmid, None)
            save_json(ap, abs_) if abs_ else ap.unlink()
        removed_total += len(dropped)
    write_manifest()
    print(f"pruned {removed_total} records")


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
    args = ap.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)
    {"weekly": cmd_weekly, "backfill": cmd_backfill, "prune": cmd_prune}[args.cmd](args)


if __name__ == "__main__":
    main()
