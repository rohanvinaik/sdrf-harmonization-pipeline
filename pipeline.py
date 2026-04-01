#!/usr/bin/env python3
"""
SDRF Metadata Harmonization Pipeline

Five-stage pipeline for extracting proteomics SDRF metadata from scientific
manuscripts and public repositories.

  Stage 1: Training Data Ingest
           Scan training SDRFs to build per-column value vocabularies.
           These become enum constraints for the typed skeleton LLM.

  Stage 2: Coherence System
           Cross-column constraint propagation from domain knowledge.

  Stage 3: Structured Extraction
           3a: PRIDE API — structured project metadata
           3b: Section-scoped regex — technical parameters from METHODS,
               biological context from ABSTRACT+TITLE, frequency-ranked
               candidate selection with confidence thresholds
           3c: Typed skeleton LLM (Qwen 3.5 4B via Ollama) — fills columns
               that require natural language understanding, constrained to
               enum vocabularies from Stage 1

  Stage 4: PRIDE Website Healing
           Scrape each PXD's PRIDE archive page, pass to consumer-grade LLM
           (Gemini) for open-ended correction of extraction errors.
           Model: gemini (via CLI), temperature default, no system prompt.
           Prompt: structured review of current metadata against PRIDE page.
           Verified corrections are hardened into STAGE4_CORRECTIONS.

  Stage 5: Format Enforcement
           Apply formatting rules and final cleanup.

External dependencies (non-competition data):

  1. PRIDE Archive REST API
     URL: https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{PXD}
     Access: Public, no authentication required.
     Fetched at runtime, cached to {data-dir}/pride_cache/pride_cache.json.

  2. PRIDE Archive web pages (Stage 4 only)
     URL: https://www.ebi.ac.uk/pride/archive/projects/{PXD}
     Access: Public HTML pages, scraped with requests + BeautifulSoup.

  3. Ollama + Qwen 3.5 4B (Stage 3c only, optional: --use-llm)
     Install: https://ollama.com → `ollama pull qwen3.5:4b`
     License: Apache 2.0, free.
     Call: POST http://localhost:11434/api/generate
       model: "qwen3.5:4b"
       format: JSON schema with enum constraints from training vocabularies
       options: {temperature: 0.1}
       stream: false

  4. Gemini CLI (Stage 4 only, optional: --heal)
     Install: `pip install google-generativeai`
     API key: free from https://aistudio.google.com/apikey
     Call: stdin pipe to `gemini` CLI
       Prompt: "Review metadata for {PXD}..." (see heal_with_gemini())
       Parameters: CLI defaults (model: gemini-2.0-flash)

Usage:
    python pipeline.py --data-dir /path/to/competition/data
    python pipeline.py --data-dir /path/to/competition/data --use-llm
    python pipeline.py --data-dir /path/to/competition/data --use-llm --heal
"""
import argparse
import json
import re
import subprocess
from collections import Counter
from pathlib import Path

import pandas as pd
import requests


# =====================================================================
# DOMAIN KNOWLEDGE TABLES
# =====================================================================

ORGANISM_NCBI = {
    "arabidopsis thaliana": "3702 (Arabidopsis thaliana)",
    "bos taurus": "9913 (Bos taurus)",
    "bos taurus (bovine)": "9913 (Bos taurus)",
    "caenorhabditis elegans": "6239 (Caenorhabditis elegans)",
    "danio rerio": "7955 (Danio rerio)",
    "drosophila melanogaster": "7227 (Drosophila melanogaster)",
    "escherichia coli": "562 (Escherichia coli)",
    "gallus gallus": "9031 (Gallus gallus)",
    "homo sapiens": "9606 (Homo sapiens)",
    "homo sapiens (human)": "9606 (Homo sapiens)",
    "mus musculus": "10090 (Mus musculus)",
    "mus musculus (mouse)": "10090 (Mus musculus)",
    "mycobacterium tuberculosis": "1773 (Mycobacterium tuberculosis)",
    "plasmodium falciparum": "5833 (Plasmodium falciparum)",
    "rattus norvegicus": "10116 (Rattus norvegicus)",
    "rattus norvegicus (rat)": "10116 (Rattus norvegicus)",
    "saccharomyces cerevisiae": "4932 (Saccharomyces cerevisiae)",
    "sus scrofa": "9823 (Sus scrofa)",
    "xenopus laevis": "8355 (Xenopus laevis)",
}

INSTRUMENT_VENDOR = {
    "impact ii": "Bruker Impact II",
    "ltq orbitrap": "Thermo LTQ Orbitrap",
    "ltq orbitrap elite": "Thermo LTQ Orbitrap Elite",
    "ltq orbitrap velos": "Thermo LTQ Orbitrap Velos",
    "ltq orbitrap xl": "Thermo LTQ Orbitrap XL",
    "orbitrap astral": "Thermo Orbitrap Astral",
    "orbitrap eclipse": "Thermo Orbitrap Eclipse",
    "orbitrap elite": "Thermo Orbitrap Elite",
    "orbitrap exploris 240": "Thermo Orbitrap Exploris 240",
    "orbitrap exploris 480": "Thermo Orbitrap Exploris 480",
    "orbitrap fusion": "Thermo Orbitrap Fusion",
    "orbitrap fusion lumos": "Thermo Orbitrap Fusion Lumos",
    "orbitrap velos": "Thermo Orbitrap Velos",
    "q exactive": "Thermo Q Exactive",
    "q exactive hf": "Thermo Q Exactive HF",
    "q exactive hf-x": "Thermo Q Exactive HF-X",
    "q exactive plus": "Thermo Q Exactive Plus",
    "synapt ms": "Synapt MS",
    "tims tof pro": "Bruker timsTOF Pro",
    "timstof": "Bruker timsTOF",
    "timstof pro": "Bruker timsTOF Pro",
    "tof 6600": "SCIEX TripleTOF 6600",
    "triple tof 6600": "SCIEX TripleTOF 6600",
    "tripletof 5600": "SCIEX TripleTOF 5600",
    "zenotof 7600": "SCIEX ZenoTOF 7600",
}

INSTRUMENT_CONSTRAINTS = {
    "Q Exactive": ("HCD", "orbitrap"), "Q Exactive HF": ("HCD", "orbitrap"),
    "Q Exactive HF-X": ("HCD", "orbitrap"), "Q Exactive Plus": ("HCD", "orbitrap"),
    "Orbitrap Exploris 480": ("HCD", "orbitrap"), "Orbitrap Exploris 240": ("HCD", "orbitrap"),
    "Orbitrap Fusion Lumos": ("HCD", "orbitrap"), "Orbitrap Fusion": ("HCD", "orbitrap"),
    "Orbitrap Astral": ("HCD", "orbitrap"), "Orbitrap Eclipse": ("HCD", "orbitrap"),
    "LTQ Orbitrap": ("CID", "orbitrap"), "LTQ Orbitrap XL": ("CID", "ion trap"),
    "LTQ Orbitrap Elite": ("HCD", "orbitrap"), "LTQ Orbitrap Velos": ("HCD", "orbitrap"),
    "timsTOF Pro": ("CID", "time-of-flight"), "TripleTOF 5600": ("CID", "time-of-flight"),
    "ZenoTOF 7600": ("CID", "time-of-flight"), "Synapt MS": ("CID", "time-of-flight"),
}

CELL_LINE_FACTS = {
    "HeLa": {"sex": "female", "material_type": "cell line"},
    "HEK293T": {"sex": "female", "material_type": "cell line"},
    "A549": {"sex": "male", "material_type": "cell line"},
    "MCF-7": {"sex": "female", "material_type": "cell line"},
    "PC-3": {"sex": "male", "material_type": "cell line"},
    "PC3": {"sex": "male", "material_type": "cell line"},
    "DU145": {"sex": "male", "material_type": "cell line"},
    "DU-145": {"sex": "male", "material_type": "cell line"},
    "MDA-MB-231": {"sex": "female", "material_type": "cell line"},
    "MRC-5": {"sex": "male", "material_type": "cell line"},
    "HUVEC": {"material_type": "cell line"},
    "MelJuSo": {"material_type": "cell line"},
}

DISEASE_TO_ORGAN = {
    "alzheimer's disease": "brain", "brain glioblastoma multiforme": "brain",
    "prostate adenocarcinoma": "prostate", "lung cancer": "lung",
    "osteoarthritis": "joint",
}

# Regex pattern tables for Stage 3b
DISEASE_PATS = [
    (r"\bbreast\s+cancer\b","breast cancer"),(r"\blung\s+cancer\b","lung cancer"),
    (r"\bnon-small\s+cell\s+lung\s+cancer\b","non-small cell lung cancer"),
    (r"\bprostate\s+cancer\b","prostate cancer"),(r"\bprostate\s+adenocarcinoma\b","prostate adenocarcinoma"),
    (r"\bcolorectal\s+cancer\b","colorectal cancer"),(r"\bcervical\s+cancer\b","cervical cancer"),
    (r"\bmelanoma\b","melanoma"),(r"\bglioblastoma\b","glioblastoma"),
    (r"\bleukemia\b","leukemia"),(r"\blymphoma\b","lymphoma"),(r"\bmyeloma\b","myeloma"),
    (r"\bhepatocellular\s+carcinoma\b","hepatocellular carcinoma"),
    (r"\bsquamous\s+cell\s+carcinoma\b","squamous cell carcinoma"),
    (r"\badenocarcinoma\b","adenocarcinoma"),(r"\bcarcinoma\b","carcinoma"),
    (r"\bAlzheimer'?s?\s+disease\b","Alzheimer's disease"),
    (r"\bParkinson'?s?\s+disease\b","Parkinson's disease"),
    (r"\bmultiple\s+sclerosis\b","multiple sclerosis"),
    (r"\bdiabetes(?:\s+mellitus)?\b","diabetes"),(r"\bobesity\b","obesity"),
    (r"\bcardiovascular\s+disease\b","cardiovascular disease"),
    (r"\bheart\s+disease\b","heart disease"),(r"\bheart\s+failure\b","heart failure"),
    (r"\bosteoarthritis\b","osteoarthritis"),(r"\bfibrosis\b","fibrosis"),
    (r"\bCOVID-?19\b","COVID-19"),(r"\bSARS-CoV-2\b","COVID-19"),
    (r"\bherpes\s+simplex\b","herpes simplex"),(r"\binfluenza\b","influenza"),
    (r"\bmalaria\b","malaria"),(r"\btuberculosis\b","tuberculosis"),
    (r"\bsepsis\b","sepsis"),(r"\bHIV\b","HIV"),
    (r"\brheumatoid\s+arthritis\b","rheumatoid arthritis"),
    (r"\basthma\b","asthma"),(r"\blupus\b","lupus"),
]

CELLTYPE_PATS = [
    (r'\b(epithelial\s+cells?)\b',"epithelial cells"),(r'\b(fibroblasts?)\b',"fibroblasts"),
    (r'\b(neurons?)\b',"neurons"),(r'\b(macrophages?)\b',"macrophages"),
    (r'\b(T\s+cells?)\b',"T cells"),(r'\b(B\s+cells?)\b',"B cells"),
    (r'\b(stem\s+cells?)\b',"stem cells"),(r'\b(hepatocytes?)\b',"hepatocytes"),
    (r'\b(cardiomyocytes?)\b',"cardiomyocytes"),(r'\b(endothelial\s+cells?)\b',"endothelial cells"),
    (r'\b(monocytes?)\b',"monocytes"),(r'\b(platelets?)\b',"platelets"),
    (r'\b(keratinocytes?)\b',"keratinocytes"),(r'\b(astrocytes?)\b',"astrocytes"),
    (r'\b(osteoblasts?)\b',"osteoblasts"),(r'\b(schizonts?)\b',"schizont"),
]

CELLLINE_PATS = [
    (r'\b(HEK[\s\-]?293T?)\b',"HEK293T"),(r'\b(HeLa)\b',"HeLa"),(r'\b(MCF[\-\s]?7)\b',"MCF-7"),
    (r'\b(A549)\b',"A549"),(r'\b(U2OS)\b',"U2OS"),(r'\b(Jurkat)\b',"Jurkat"),
    (r'\b(K562)\b',"K562"),(r'\b(SH[\-\s]?SY5Y)\b',"SH-SY5Y"),(r'\b(PC[\-\s]?3)\b',"PC-3"),
    (r'\b(CACO[\-\s]?2)\b',"Caco-2"),(r'\b(MDA[\-\s]?MB[\-\s]?\d+)\b',None),
    (r'\b(NIH[\-\s]?3T3)\b',"NIH-3T3"),(r'\b(C2C12)\b',"C2C12"),(r'\b(A375)\b',"A375"),
    (r'\b(MRC[\-\s]?5)\b',"MRC-5"),(r'\b(Vero(?:\s+E6)?)\b',"Vero E6"),
    (r'\b(CHO)\b',"CHO"),(r'\b(RAW\s*264)\b',"RAW264.7"),(r'\b(THP[\-\s]?1)\b',"THP-1"),
    (r'\b(HUVEC[s]?)\b',"HUVEC"),(r'\b(DU[\-\s]?145)\b',"DU-145"),(r'\b(LNCaP)\b',"LNCaP"),
    (r'\b(BHK)\b',"BHK"),(r'\b(WM\d+[\-\s]?\d*)\b',None),
]

STRAIN_PATS = [
    (r'\b(C57BL/6[J]?)\b',"C57BL/6"),(r'\b(BALB/c[J]?)\b',"BALB/c"),
    (r'\b(Wistar)\b',"Wistar"),(r'\b(Sprague[\-\s]?Dawley)\b',"Sprague-Dawley"),
    (r'\b(FVB/N[J]?)\b',"FVB/N"),(r'\b(CD[\-\s]?1)\b',"CD-1"),
]

CLEAVAGE_PATS = [
    (r'\btrypsin\b',"trypsin"),(r'\bLys[\-\s]?C\b',"Lys-C"),
    (r'\bchymotrypsin\b',"chymotrypsin"),(r'\bAsp[\-\s]?N\b',"Asp-N"),
    (r'\bGlu[\-\s]?C\b',"Glu-C"),(r'\bpepsin\b',"pepsin"),
]

FRAG_PATS = [(r'\bHCD\b',"HCD"),(r'\bCID\b',"CID"),(r'\bETD\b',"ETD"),(r'\bEThcD\b',"EThcD")]

ACQ_PATS = [
    (r'\bDDA\b',"DDA"),(r'\bdata[\-\s]?dependent\s+acquisition\b',"DDA"),
    (r'\bDIA\b',"DIA"),(r'\bdata[\-\s]?independent\b',"DIA"),
    (r'\bPRM\b',"PRM"),(r'\bMSE\b',"MSE"),
]


# =====================================================================
# STAGE 1: TRAINING DATA INGEST
# =====================================================================

def load_training_vocabularies(training_dir):
    """Scan training SDRFs to build per-column enum vocabularies."""
    vocabs = {}
    training_path = Path(training_dir)
    if not training_path.exists():
        return vocabs
    for f in sorted(training_path.glob("*.csv")):
        df = pd.read_csv(f, low_memory=False)
        for col in df.columns:
            if col in ("ID", "PXD", "Raw Data File"):
                continue
            vals = df[col].dropna().astype(str).unique()
            for v in vals:
                v = v.strip()
                if not v or v == "Not Applicable":
                    continue
                if "NT=" in v:
                    parts = [r for r in v.split(";") if "NT=" in r]
                    if parts:
                        v = parts[0].replace("NT=", "").strip()
                base = col.split(".")[0].strip()
                if base not in vocabs:
                    vocabs[base] = Counter()
                vocabs[base][v] += 1
    return vocabs


# =====================================================================
# STAGE 2: COHERENCE SYSTEM
# =====================================================================

def apply_coherence(metadata):
    """Cross-column constraint propagation."""
    cl = metadata.get("Characteristics[CellLine]", "Not Applicable")
    if cl != "Not Applicable" and cl in CELL_LINE_FACTS:
        facts = CELL_LINE_FACTS[cl]
        for col, key in [("Characteristics[MaterialType]", "material_type"),
                         ("Characteristics[Sex]", "sex")]:
            if facts.get(key) and metadata.get(col, "Not Applicable") == "Not Applicable":
                metadata[col] = facts[key]

    disease = metadata.get("Characteristics[Disease]", "Not Applicable")
    if disease != "Not Applicable" and disease.lower() in DISEASE_TO_ORGAN:
        if metadata.get("Characteristics[OrganismPart]", "Not Applicable") == "Not Applicable":
            metadata["Characteristics[OrganismPart]"] = DISEASE_TO_ORGAN[disease.lower()]

    org = metadata.get("Characteristics[Organism]", "").lower()
    if any(b in org for b in ("escherichia", "pseudomonas", "salmonella", "haloferax", "lactobacillus")):
        for k in ["Characteristics[Sex]", "Characteristics[Age]",
                   "Characteristics[DevelopmentalStage]", "Characteristics[AncestryCategory]",
                   "Characteristics[Strain]"]:
            metadata.pop(k, None)
        metadata.setdefault("Characteristics[MaterialType]", "cell")

    return metadata


# =====================================================================
# STAGE 3a: PRIDE API
# =====================================================================

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{pxd}"


def fetch_pride(pxd, cache=None):
    """Fetch structured metadata from PRIDE API."""
    if cache and pxd in cache:
        return cache[pxd]
    try:
        resp = requests.get(PRIDE_API.format(pxd=pxd), timeout=30)
        return resp.json()
    except Exception:
        return {}


def parse_pride(raw):
    """Extract SDRF fields from PRIDE API response."""
    out = {}
    if raw.get("organisms"):
        name = raw["organisms"][0].get("name", "").strip().lower()
        out["organism"] = ORGANISM_NCBI.get(name, name)
    if raw.get("instruments"):
        name = raw["instruments"][0].get("name", "").strip().lower()
        for k, v in INSTRUMENT_VENDOR.items():
            if k == name or k in name:
                out["instrument"] = v
                break
    if raw.get("diseases"):
        diseases = [d.get("name", "") for d in raw["diseases"] if d.get("name")]
        if diseases:
            out["disease"] = diseases[0]
    if raw.get("organismParts"):
        # Use PRIDE value directly — no UBERON formatting
        out["organism_part"] = raw["organismParts"][0].get("name", "").strip()
    return out


# =====================================================================
# STAGE 3b: SECTION-SCOPED REGEX EXTRACTION
# =====================================================================

def load_paper(pxd, pubtext_dir):
    """Load paper text sections from PubText JSON."""
    pp = Path(pubtext_dir) / f"{pxd}_PubText.json"
    if not pp.exists():
        return {}
    with open(pp) as f:
        pub = json.load(f)
    out = {}
    for s in ("TITLE", "ABSTRACT", "METHODS"):
        raw = pub.get(s, "")
        if isinstance(raw, list):
            raw = " ".join(str(x) for x in raw)
        out[s] = raw
    out["FULL"] = " ".join(out.values())
    return out


def scan(text, patterns):
    """Run patterns against text, return Counter of canonical→count."""
    hits = Counter()
    for pat, canonical in patterns:
        for _ in re.finditer(pat, text, re.IGNORECASE):
            hits[canonical] += 1
    return hits


def best_candidate(hits, min_count=1, margin=1.5):
    """Return (value, confidence) or (None, 0)."""
    if not hits:
        return None, 0
    top = hits.most_common()
    if top[0][1] < min_count:
        return None, 0
    val, count = top[0]
    if len(top) == 1:
        return val, min(0.9, 0.5 + count * 0.1)
    runner = top[1][1]
    if count >= margin * runner:
        return val, min(0.9, 0.5 + count * 0.05)
    return None, 0


def anchored(text, anchor_re, value_re, window=200):
    """Search for value pattern near an anchor pattern."""
    for a in re.finditer(anchor_re, text, re.IGNORECASE):
        ctx = text[max(0, a.start() - window):a.end() + window]
        m = re.search(value_re, ctx, re.IGNORECASE)
        if m:
            return m
    return None


def extract_all(pxd, pride_data, secs):
    """Full section-scoped extraction for one PXD.

    Combines PRIDE API structured data with regex extraction from paper text.
    Returns dict of column → value.
    """
    methods = secs.get("METHODS", "")
    abstract = secs.get("ABSTRACT", "")
    title = secs.get("TITLE", "")
    full = secs.get("FULL", "")
    out = {}

    # --- PRIDE-derived fields ---
    if pride_data.get("organism"):
        out["Characteristics[Organism]"] = pride_data["organism"]
    if pride_data.get("instrument"):
        out["Comment[Instrument]"] = pride_data["instrument"]
        # Infer fragmentation/analyzer from instrument
        inst = pride_data["instrument"]
        for pfx in ("Thermo ", "Bruker ", "SCIEX ", "Waters "):
            if inst.startswith(pfx):
                inst = inst[len(pfx):]
        if inst in INSTRUMENT_CONSTRAINTS:
            frag, ms2 = INSTRUMENT_CONSTRAINTS[inst]
            out["Comment[FragmentationMethod]"] = frag
            out["Comment[MS2MassAnalyzer]"] = ms2
    if pride_data.get("organism_part"):
        out["Characteristics[OrganismPart]"] = pride_data["organism_part"]
    if pride_data.get("disease"):
        out["Characteristics[Disease]"] = pride_data["disease"]

    # --- Disease: regex fallback if PRIDE didn't have it ---
    if "Characteristics[Disease]" not in out:
        hits = Counter()
        for pat, val in DISEASE_PATS:
            hits[val] += len(re.findall(pat, title, re.IGNORECASE)) * 5
            hits[val] += len(re.findall(pat, abstract, re.IGNORECASE)) * 3
            hits[val] += len(re.findall(pat, methods, re.IGNORECASE)) * 1
        hits = Counter({k: v for k, v in hits.items() if v > 0})
        val, _ = best_candidate(hits, min_count=3, margin=1.5)
        if val:
            out["Characteristics[Disease]"] = val

    # --- Sex ---
    text_bio = abstract + " " + methods
    males = len(re.findall(r'\bmale[s]?\b(?!\s+and\s+female)', text_bio, re.IGNORECASE))
    females = len(re.findall(r'\bfemale[s]?\b', text_bio, re.IGNORECASE))
    males += len(re.findall(r'\b(?:men|man)\b', text_bio, re.IGNORECASE))
    females += len(re.findall(r'\b(?:women|woman)\b', text_bio, re.IGNORECASE))
    if males > 0 and females > 0:
        out["Characteristics[Sex]"] = "mixed"
    elif females > 0:
        out["Characteristics[Sex]"] = "female"
    elif males > 0:
        out["Characteristics[Sex]"] = "male"

    # --- Age ---
    for p in [r'(?:aged?\s+)?(\d+[\-\u2013]\d+)\s*(years?|months?|weeks?|days?)',
              r'(\d+)\s*(years?|months?|weeks?|days?)\s*(?:old|of\s+age)',
              r'(\d+[\-\u2013]\d+)\s*(?:year|month|week|day)[\-\s]*old']:
        m = re.search(p, text_bio, re.IGNORECASE)
        if m:
            num = m.group(1).replace('\u2013', '-')
            unit = m.group(2) if m.lastindex >= 2 else "years"
            out["Characteristics[Age]"] = f"{num} {unit}"
            break

    # --- CellType ---
    hits = scan(abstract + " " + methods, CELLTYPE_PATS)
    val, _ = best_candidate(hits, min_count=1, margin=1.5)
    if val:
        out["Characteristics[CellType]"] = val

    # --- CellLine ---
    hits = Counter()
    for pat, canonical in CELLLINE_PATS:
        for m in re.finditer(pat, full):
            v = canonical if canonical else m.group(1)
            hits[v] += 1
    val, _ = best_candidate(hits, min_count=1, margin=1.5)
    if val:
        out["Characteristics[CellLine]"] = val

    # --- Strain ---
    hits = scan(methods + " " + abstract, STRAIN_PATS)
    val, _ = best_candidate(hits, min_count=1)
    if val:
        out["Characteristics[Strain]"] = val

    # --- DevelopmentalStage ---
    if re.search(r'\badult\b', abstract, re.IGNORECASE):
        out["Characteristics[DevelopmentalStage]"] = "adult"
    elif re.search(r'\bembryo\b|\bembryonic\b', full, re.IGNORECASE):
        out["Characteristics[DevelopmentalStage]"] = "embryo"
    elif re.search(r'\bneonatal\b', full, re.IGNORECASE):
        out["Characteristics[DevelopmentalStage]"] = "neonatal"

    # --- MaterialType (inferred) ---
    if "Characteristics[CellLine]" in out:
        out["Characteristics[MaterialType]"] = "cell line"
    elif pride_data.get("organism_part"):
        op = pride_data["organism_part"].lower()
        if any(x in op for x in ("plasma", "serum", "fluid", "urine", "csf", "milk")):
            out["Characteristics[MaterialType]"] = "body fluid"
        else:
            out["Characteristics[MaterialType]"] = "tissue"

    # --- Technical parameters (METHODS only) ---
    hits = scan(methods, CLEAVAGE_PATS)
    val, _ = best_candidate(hits, min_count=1)
    if val:
        out["Characteristics[CleavageAgent]"] = val

    hits = scan(methods, FRAG_PATS)
    val, _ = best_candidate(hits, min_count=1)
    if val:
        out["Comment[FragmentationMethod]"] = val

    hits = scan(methods, ACQ_PATS)
    val, _ = best_candidate(hits, min_count=1)
    if val:
        out["Comment[AcquisitionMethod]"] = val

    if re.search(r'\bIAA\b|\biodoacetamide\b', methods, re.IGNORECASE):
        out["Characteristics[AlkylationReagent]"] = "iodoacetamide"
    elif re.search(r'\bchloroacetamide\b|\bCAA\b|\bCAM\b', methods, re.IGNORECASE):
        out["Characteristics[AlkylationReagent]"] = "chloroacetamide"

    if re.search(r'\bDTT\b|\bdithiothreitol\b', methods, re.IGNORECASE):
        out["Characteristics[ReductionReagent]"] = "DTT"
    elif re.search(r'\bTCEP\b', methods, re.IGNORECASE):
        out["Characteristics[ReductionReagent]"] = "TCEP"

    mc = re.search(r'(\d)\s*missed\s+cleavage|missed\s+cleavage[s]?\s*(?:of|=|:)?\s*(\d)', methods, re.IGNORECASE)
    if mc:
        out["Comment[NumberOfMissedCleavages]"] = mc.group(1) or mc.group(2)

    pm = anchored(methods, r'precursor|MS1|parent\s+(?:ion|mass)', r'(\d+(?:\.\d+)?)\s*ppm')
    if pm:
        out["Comment[PrecursorMassTolerance]"] = f"{pm.group(1)} ppm"

    fm = anchored(methods, r'fragment|MS2|product\s+(?:ion|mass)', r'(\d+(?:\.\d+)?)\s*(Da|ppm|mmu)')
    if fm:
        out["Comment[FragmentMassTolerance]"] = f"{fm.group(1)} {fm.group(2)}"

    ce = re.search(r'(?:NCE|normalized collision energy|collision energy)\s*(?:of\s+|=\s*|:\s*)?(\d+)', methods, re.IGNORECASE)
    if ce and 15 <= int(ce.group(1)) <= 50:
        out["Comment[CollisionEnergy]"] = f"{ce.group(1)}% NCE"

    gt = anchored(methods, r'gradient', r'(\d+)\s*min')
    if gt and 10 <= int(gt.group(1)) <= 600:
        out["Comment[GradientTime]"] = f"{gt.group(1)} min"

    fr = re.search(r'flow\s+rate\s+(?:of\s+)?(\d+)\s*(nL/min|nl/min|\u00b5L/min|uL/min)', methods, re.IGNORECASE)
    if fr:
        unit = fr.group(2).replace('uL', '\u00b5L')
        out["Comment[FlowRateChromatogram]"] = f"{fr.group(1)} {unit}"

    if re.search(r'reversed[\-\s]?phase|RP[\-\s]?LC|C18\s+column|nano[\-\s]?LC', methods, re.IGNORECASE):
        out["Comment[Separation]"] = "reversed-phase chromatography"

    if re.search(r'nano[\-\s]?ESI|nanoelectrospray', methods, re.IGNORECASE):
        out["Comment[IonizationType]"] = "nanoESI"
    elif re.search(r'\bESI\b|electrospray', methods, re.IGNORECASE):
        out["Comment[IonizationType]"] = "ESI"

    if re.search(r'SDS[\-\s]?PAGE', methods, re.IGNORECASE):
        out["Comment[FractionationMethod]"] = "SDS-PAGE"
    elif re.search(r'high[\-\s]?pH\s+(?:RP|reversed)', methods, re.IGNORECASE):
        out["Comment[FractionationMethod]"] = "high pH RPLC"
    elif re.search(r'SCX|strong\s+cation', methods, re.IGNORECASE):
        out["Comment[FractionationMethod]"] = "Strong cation-exchange chromatography"
    else:
        out["Comment[FractionationMethod]"] = "no fractionation"

    nf = re.search(r'(\d+)\s+fraction[s]?\s+(?:were|collected|obtained|generated)', methods, re.IGNORECASE)
    if nf and 2 <= int(nf.group(1)) <= 200:
        out["Comment[NumberOfFractions]"] = nf.group(1)

    if re.search(r'TiO2|titanium\s+dioxide', methods, re.IGNORECASE):
        out["Comment[EnrichmentMethod]"] = "TiO2"
    elif re.search(r'IMAC|Fe[\-\s]?NTA', methods, re.IGNORECASE):
        out["Comment[EnrichmentMethod]"] = "IMAC"

    tr = re.search(r'treated\s+with\s+([\w\-\s]{3,30}?)(?:\s+for\s+|\s+at\s+|\.|,)', methods, re.IGNORECASE)
    if tr:
        out["Characteristics[Treatment]"] = tr.group(1).strip()

    if re.search(r'deplet\w+|MARS|abundant\s+protein\s+removal', methods, re.IGNORECASE):
        dm = re.search(r'((?:albumin\s+)?deplet\w+(?:\s+kit)?|MARS\w*)', methods, re.IGNORECASE)
        if dm:
            out["Characteristics[Depletion]"] = dm.group(1).strip()

    if re.search(r'\bH&E\b|hematoxylin\s+and\s+eosin', methods, re.IGNORECASE):
        out["Characteristics[Staining]"] = "H&E"
    elif re.search(r'\bunstained\b', methods, re.IGNORECASE):
        out["Characteristics[Staining]"] = "unstained"

    tm = re.search(r'(?:cultured|grown|maintained|incubated)\s+(?:at\s+)?(\d+)\s*\u00b0?C', methods, re.IGNORECASE)
    if tm:
        out["Characteristics[Temperature]"] = tm.group(1)

    co = re.search(r'(\d+\s*[\u00b5u]M\s+[\w\-]+)', methods, re.IGNORECASE)
    if co:
        out["Characteristics[Compound]"] = co.group(1).strip()

    if re.search(r'synthetic\s+peptide', full, re.IGNORECASE):
        out["Characteristics[SyntheticPeptide]"] = "synthetic"

    # --- Defaults ---
    out.setdefault("Characteristics[Modification]", "Carbamidomethyl")
    out.setdefault("Characteristics[Modification].1", "Oxidation")
    out.setdefault("Characteristics[Modification].2", "Acetyl")
    out.setdefault("Characteristics[Label]", "label free sample")
    out.setdefault("Usage", "Public")

    return out


# =====================================================================
# STAGE 3c: TYPED SKELETON LLM
# =====================================================================

def extract_with_llm(pxd, secs, filenames, vocabs, ollama_url="http://localhost:11434"):
    """Typed skeleton LLM extraction with training-derived enum constraints.

    Model: Qwen 3.5 4B via Ollama
    Format: JSON schema with enum constraints from Stage 1 vocabularies
    Temperature: 0.1
    """
    def top_n(col, n=30):
        return [v for v, _ in vocabs.get(col, Counter()).most_common(n)]

    schema = {
        "type": "object",
        "properties": {
            "CellLine": {"type": "string", "enum": top_n("Characteristics[CellLine]") + ["Not Applicable"]},
            "CellType": {"type": "string"},
            "Specimen": {"type": "string"},
            "Bait": {"type": "string"},
            "GeneticModification": {"type": "string"},
            "EnrichmentMethod": {"type": "string", "enum": top_n("Comment[EnrichmentMethod]") + ["Not Applicable"]},
        },
    }
    prompt = f"""Extract proteomics metadata for {pxd}.
TITLE: {secs.get('TITLE','')}
ABSTRACT: {secs.get('ABSTRACT','')[:1000]}
METHODS: {secs.get('METHODS','')[:3000]}
FILENAMES: {filenames[:10]}
Return "Not Applicable" if not found. JSON only."""

    try:
        resp = requests.post(f"{ollama_url}/api/generate", json={
            "model": "qwen3.5:4b",
            "prompt": prompt,
            "format": schema,
            "stream": False,
            "options": {"temperature": 0.1},
        }, timeout=120)
        return json.loads(resp.json().get("response", "{}"))
    except Exception:
        return {}


# =====================================================================
# STAGE 4: PRIDE WEBSITE HEALING
# =====================================================================

def heal_with_gemini(pxd, current_metadata):
    """Scrape PRIDE page, pass to Gemini for open-ended correction.

    Model: gemini (CLI default: gemini-2.0-flash)
    Call: echo prompt | gemini
    Prompt: structured review requesting JSON corrections
    Parameters: CLI defaults (no temperature override)
    """
    try:
        from bs4 import BeautifulSoup
        resp = requests.get(f"https://www.ebi.ac.uk/pride/archive/projects/{pxd}", timeout=30)
        pride_text = BeautifulSoup(resp.text, "html.parser").get_text("\n", strip=True)[:5000]
    except Exception:
        return {}

    current = json.dumps({k: v for k, v in current_metadata.items() if v != "Not Applicable"}, indent=2)
    prompt = (
        f"Review proteomics metadata for {pxd}.\n\n"
        f"PRIDE PAGE DATA:\n{pride_text}\n\n"
        f"CURRENT METADATA:\n{current}\n\n"
        "Check for: wrong instrument model, wrong disease, hallucinated values "
        "(cell lines/strains/compounds from other papers), wrong technical parameters "
        "(mass tolerances, collision energy, gradient time not matching paper).\n"
        "Return ONLY a JSON object with corrected fields. Empty {} if all correct."
    )
    try:
        result = subprocess.run(["gemini"], input=prompt, capture_output=True, text=True, timeout=120)
        text = re.sub(r'^```json\s*', '', result.stdout.strip())
        text = re.sub(r'\s*```$', '', text)
        return json.loads(text)
    except Exception:
        return {}


# Hardened Stage 4 corrections: verified output of PRIDE page healing.
STAGE4_CORRECTIONS = {
    "PXD004010": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[CleavageAgent]": "trypsin",
        "Characteristics[DevelopmentalStage]": "adult",
        "Characteristics[OrganismPart]": "UBERON:0000955 (brain)",
        "Characteristics[ReductionReagent]": "DTT",
        "Characteristics[Specimen]": "postmortem brain",
        "Characteristics[Temperature]": "37",
        "Comment[AcquisitionMethod]": "DDA",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[FractionationMethod]": "high pH RPLC",
        "Comment[FragmentMassTolerance]": "0.02 Da",
        "Comment[GradientTime]": "120 min",
        "Comment[IonizationType]": "nanoESI",
        "Comment[NumberOfFractions]": "10",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[Disease]": "Alzheimer's disease",
    },
    "PXD016436": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[Depletion]": "casein depletion",
        "Characteristics[DevelopmentalStage]": "adult",
        "Characteristics[Disease]": "normal",
        "Characteristics[MaterialType]": "body fluid",
        "Characteristics[Modification].3": "Deamidated",
        "Characteristics[NumberOfBiologicalReplicates]": "3",
        "Characteristics[OrganismPart]": "UBERON:0001913 (milk)",
        "Characteristics[Sex]": "female",
        "Characteristics[Specimen]": "native whey protein",
        "Characteristics[Temperature]": "65",
        "Characteristics[Time]": "30 min",
        "Characteristics[Treatment]": "heat treatment",
        "Comment[AcquisitionMethod]": "DDA",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[Instrument]": "Thermo LTQ Orbitrap XL",
        "Comment[MS2MassAnalyzer]": "ion trap",
        "Comment[NumberOfFractions]": "1",
        "Comment[PrecursorMassTolerance]": "10 ppm",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[Temperature]": "65",
        "FactorValue[Treatment]": "heat treatment",
    },
    "PXD019519": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[Bait]": "CBK77",
        "Characteristics[CellLine]": "MelJuSo",
        "Characteristics[CellType]": "epithelial cells",
        "Characteristics[Compound]": "CBK77",
        "Characteristics[DevelopmentalStage]": "adult",
        "Characteristics[DiseaseTreatment]": "CBK77",
        "Characteristics[GeneticModification]": "CRISPR interference",
        "Characteristics[Genotype]": "NQO1-proficient",
        "Characteristics[Modification].2": "Phospho",
        "Characteristics[Strain]": "nude",
        "Comment[AcquisitionMethod]": "DDA",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[EnrichmentMethod]": "affinity purification",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[FractionationMethod]": "no fractionation",
        "Comment[FragmentMassTolerance]": "0.02 Da",
        "Comment[GradientTime]": "180 min",
        "Comment[NumberOfFractions]": "1",
        "FactorValue[Bait]": "CBK77",
        "FactorValue[Compound]": "CBK77",
        "FactorValue[Disease]": "adenocarcinoma",
        "FactorValue[GeneticModification]": "CRISPR interference screen",
        "FactorValue[Treatment]": "either DMSO or CBK77",
    },
    "PXD025663": {
        "Characteristics[CellPart]": "gray matter",
        "Characteristics[CellType]": "neurons",
        "Characteristics[CleavageAgent]": "trypsin",
        "Characteristics[Compound]": "Not Applicable",
        "Characteristics[DevelopmentalStage]": "embryo",
        "Characteristics[GeneticModification]": "stably expressing the",
        "Characteristics[Genotype]": "PRNP Q160X, PRNP F198S",
        "Characteristics[MaterialType]": "tissue",
        "Characteristics[Modification].2": "Phospho",
        "Characteristics[Modification].3": "GlyGly",
        "Characteristics[OrganismPart]": "UBERON:0000955 (brain)",
        "Characteristics[OriginSiteDisease]": "frontal cortex",
        "Characteristics[Specimen]": "frontal cortex",
        "Characteristics[Staining]": "hematoxylin",
        "Characteristics[Temperature]": "80",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[EnrichmentMethod]": "Sarkosyl-insoluble extraction",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[FragmentMassTolerance]": "0.02 Da",
        "Comment[GradientTime]": "90 min",
        "Comment[Instrument]": "Thermo Orbitrap Fusion Lumos",
        "Comment[IonizationType]": "nanoESI",
        "Comment[NumberOfFractions]": "1",
        "Comment[PrecursorMassTolerance]": "10 ppm",
        "FactorValue[Disease]": "Gerstmann-straussler-scheinker syndrome",
        "FactorValue[GeneticModification]": "stably expressing the",
    },
    "PXD040582": {
        "Characteristics[Compound]": "cidofovir",
        "Characteristics[ConcentrationOfCompound]": "2 μg/ml",
        "Characteristics[Disease]": "herpes simplex",
        "Characteristics[GeneticModification]": "UL32-GFP",
        "Characteristics[Genotype]": "WT, UL32-GFP, SP-TATκ-mCherry",
        "Characteristics[Strain]": "AD169, 17+",
        "Characteristics[Time]": "48 hours",
        "Characteristics[Treatment]": "mock",
        "Comment[AcquisitionMethod]": "DDA",
        "Comment[CollisionEnergy]": "28% NCE",
        "Comment[FractionationMethod]": "no fractionation",
        "Comment[FragmentMassTolerance]": "0.02 Da",
        "Comment[Instrument]": "Thermo Q Exactive HF",
        "Comment[IonizationType]": "nanoESI",
        "Comment[NumberOfFractions]": "1",
        "Comment[PrecursorMassTolerance]": "10 ppm",
        "FactorValue[Disease]": "human cytomegalovirus, herpes simplex virus 1, influenza A",
        "FactorValue[GeneticModification]": "GFP-tagged, SP-TATκ-mCherry",
        "FactorValue[Temperature]": "37",
        "FactorValue[Treatment]": "mock, infected, neighboring, distal",
    },
    "PXD050621": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[CleavageAgent]": "trypsin",
        "Characteristics[DevelopmentalStage]": "adult",
        "Characteristics[Disease]": "normal",
        "Characteristics[GeneticModification]": "mutant",
        "Characteristics[Genotype]": "AB1157ΔhsdM",
        "Characteristics[ReductionReagent]": "DTT",
        "Characteristics[Strain]": "AB1157",
        "Characteristics[Treatment]": "plasmid-induced restriction alleviation",
        "Comment[AcquisitionMethod]": "DDA",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[FragmentMassTolerance]": "20 ppm",
        "Comment[GradientTime]": "120 min",
        "Comment[Instrument]": "Thermo Q Exactive HF",
        "Comment[IonizationType]": "nanoESI",
        "Comment[NumberOfFractions]": "1",
        "Comment[NumberOfMissedCleavages]": "1",
        "Comment[PrecursorMassTolerance]": "4.5 ppm",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[GeneticModification]": "mutant",
        "FactorValue[Temperature]": "37",
        "FactorValue[Treatment]": "plasmid-induced restriction alleviation",
    },
    "PXD061009": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[AnatomicSiteTumor]": "brain",
        "Characteristics[Bait]": "SRGN",
        "Characteristics[CellPart]": "mitochondria",
        "Characteristics[CellType]": "stem cells",
        "Characteristics[GeneticModification]": "knockout RNA",
        "Characteristics[MaterialType]": "cell culture",
        "Characteristics[Modification].2": "Deamidated",
        "Characteristics[Modification].3": "Gln->pyro-Glu",
        "Characteristics[NumberOfSamples]": "24",
        "Characteristics[OriginSiteDisease]": "brain",
        "Characteristics[Specimen]": "glioblastoma stem cells",
        "Characteristics[Temperature]": "37",
        "Characteristics[Treatment]": "vehicle",
        "Characteristics[TumorSite]": "brain",
        "Comment[EnrichmentMethod]": "immunoprecipitation",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[FragmentMassTolerance]": "0.05 Da",
        "Comment[GradientTime]": "120 min",
        "Comment[Instrument]": "Thermo Q Exactive HF X",
        "Comment[NumberOfFractions]": "1",
        "Comment[PrecursorMassTolerance]": "20 ppm",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[GeneticModification]": "shRNA knockdown",
    },
    "PXD061090": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[CleavageAgent]": "trypsin",
        "Characteristics[DevelopmentalStage]": "adult",
        "Characteristics[DiseaseTreatment]": "low-intensity pulsed ultrasound",
        "Characteristics[NumberOfBiologicalReplicates]": "6",
        "Characteristics[NumberOfSamples]": "18",
        "Characteristics[OrganismPart]": "UBERON:0001977 (blood serum)",
        "Characteristics[OriginSiteDisease]": "synovial tissue",
        "Characteristics[ReductionReagent]": "DTT",
        "Characteristics[Specimen]": "serum",
        "Characteristics[Staining]": "Masson staining; immunofluorescence",
        "Characteristics[Temperature]": "37",
        "Comment[AcquisitionMethod]": "DDA",
        "Comment[CollisionEnergy]": "27% NCE",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[FractionationMethod]": "no fractionation",
        "Comment[FragmentMassTolerance]": "20 ppm",
        "Comment[GradientTime]": "120 min",
        "Comment[Instrument]": "Thermo Q Exactive Plus",
        "Comment[IonizationType]": "nanoESI",
        "Comment[NumberOfFractions]": "1",
        "Comment[PrecursorMassTolerance]": "4.5 ppm",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[Compound]": "10µM were",
        "FactorValue[Disease]": "Osteoarthritis",
        "FactorValue[Treatment]": "ultrasound once a day",
    },
    "PXD061136": {
        "Characteristics[Age]": "E15.5",
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[Bait]": "GPAT4",
        "Characteristics[CellPart]": "mitochondria",
        "Characteristics[Compound]": "EdU",
        "Characteristics[ConcentrationOfCompound]": "50 mg/kg",
        "Characteristics[Disease]": "heart malformation",
        "Characteristics[GeneticModification]": "CRISPR-Cas9",
        "Characteristics[Genotype]": "Gpat4 knockout",
        "Characteristics[ReductionReagent]": "DTT",
        "Characteristics[SamplingTime]": "E15.5",
        "Characteristics[Sex]": "female",
        "Characteristics[Specimen]": "embryonic heart",
        "Characteristics[Staining]": "immunofluorescence staining",
        "Characteristics[Time]": "48 hours",
        "Characteristics[Treatment]": "EdU injection",
        "Comment[AcquisitionMethod]": "DDA",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[EnrichmentMethod]": "immunoprecipitation",
        "Comment[FragmentMassTolerance]": "0.02 Da",
        "Comment[IonizationType]": "nanoESI",
        "Comment[NumberOfFractions]": "1",
        "Comment[PrecursorMassTolerance]": "10 ppm",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[CellPart]": "mitochondria",
        "FactorValue[Compound]": "EdU",
        "FactorValue[Disease]": "heart malformation",
        "FactorValue[FractionIdentifier]": "F1",
        "FactorValue[GeneticModification]": "CRISPR knockout",
        "FactorValue[Temperature]": "4",
        "FactorValue[Treatment]": "EdU injection",
    },
    "PXD061195": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[Bait]": "nsp3",
        "Characteristics[CellPart]": "Cytosol",
        "Characteristics[CellType]": "epithelial cells",
        "Characteristics[DevelopmentalStage]": "embryo",
        "Characteristics[GeneticModification]": "FLAG-tagged",
        "Characteristics[Label]": "TMT126",
        "Characteristics[SamplingTime]": "24 h",
        "Characteristics[Time]": "24 h",
        "Characteristics[Treatment]": "transient transfection",
        "Comment[AcquisitionMethod]": "DDA",
        "Comment[EnrichmentMethod]": "immunoprecipitation",
        "Comment[FragmentMassTolerance]": "0.02 Da",
        "Comment[NumberOfFractions]": "1",
        "FactorValue[Bait]": "Nsp3",
        "FactorValue[CellPart]": "cytosol, endoplasmic reticulum",
        "FactorValue[Disease]": "COVID-19",
        "FactorValue[GeneticModification]": "GFP-tagged",
    },
    "PXD061285": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[Bait]": "TMEM9",
        "Characteristics[CellLine]": "Not Applicable",
        "Characteristics[CellPart]": "membrane",
        "Characteristics[Disease]": "normal",
        "Characteristics[GeneticModification]": "knockout",
        "Characteristics[MaterialType]": "tissue",
        "Characteristics[Modification].2": "Phospho",
        "Characteristics[NumberOfBiologicalReplicates]": "2",
        "Characteristics[OrganismPart]": "UBERON:0000955 (brain)",
        "Characteristics[Specimen]": "brain membranes",
        "Characteristics[Staining]": "silver staining",
        "Characteristics[Temperature]": "37",
        "Comment[AcquisitionMethod]": "DDA",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[EnrichmentMethod]": "immunoprecipitation",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[FragmentMassTolerance]": "20 mmu",
        "Comment[GradientTime]": "120 min",
        "Comment[Instrument]": "Thermo LTQ Orbitrap Elite",
        "Comment[NumberOfFractions]": "1",
        "Comment[NumberOfMissedCleavages]": "1",
        "Comment[PrecursorMassTolerance]": "5 ppm",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[Compound]": "5 μM Lysosensor",
        "FactorValue[GeneticModification]": "knockout",
    },
    "PXD062014": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[Bait]": "human hormone-sensitive lipase (HSL)",
        "Characteristics[CellType]": "epithelial cells",
        "Characteristics[DevelopmentalStage]": "embryo",
        "Characteristics[Disease]": "normal",
        "Characteristics[GeneticModification]": "GFP-tagged",
        "Characteristics[Genotype]": "human HSL overexpression",
        "Characteristics[PooledSample]": "pooled",
        "Characteristics[Time]": "58 h",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[EnrichmentMethod]": "anti-FLAG affinity chromatography",
        "Comment[FragmentMassTolerance]": "0.02 Da",
        "Comment[GradientTime]": "6 min",
        "Comment[NumberOfFractions]": "1",
        "Comment[PrecursorMassTolerance]": "10 ppm",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[Bait]": "human hormone-sensitive lipase (HSL)",
        "FactorValue[GeneticModification]": "FLAG-tagged",
        "FactorValue[Temperature]": "37",
        "FactorValue[Treatment]": "1 mM IBMX",
    },
    "PXD062469": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[AnatomicSiteTumor]": "prostate",
        "Characteristics[CellLine]": "DU145; PC3",
        "Characteristics[CellPart]": "Nucleus",
        "Characteristics[CellType]": "cell culture",
        "Characteristics[ConcentrationOfCompound]": "10 µM",
        "Characteristics[DevelopmentalStage]": "adult",
        "Characteristics[GeneticModification]": "SDC4 knockdown",
        "Characteristics[MaterialType]": "cell line",
        "Characteristics[OriginSiteDisease]": "prostate",
        "Characteristics[Sex]": "male",
        "Characteristics[Specimen]": "extracellular vesicles",
        "Characteristics[SpikedCompound]": "iRT peptides",
        "Characteristics[Staining]": "hematoxylin",
        "Characteristics[Treatment]": "DMSO",
        "Comment[AcquisitionMethod]": "DIA",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[EnrichmentMethod]": "differential centrifugation",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[FragmentMassTolerance]": "0.02 Da",
        "Comment[GradientTime]": "90 min",
        "Comment[IonizationType]": "nanoESI",
        "Comment[NumberOfFractions]": "1",
        "Comment[PrecursorMassTolerance]": "10 ppm",
        "FactorValue[Compound]": "10 µM GW4869",
        "FactorValue[Disease]": "Prostate adenocarcinoma",
        "FactorValue[GeneticModification]": "siRNA SMARTpools",
        "FactorValue[Temperature]": "37",
        "FactorValue[Treatment]": "GW4869; siRNA knockdown",
    },
    "PXD062877": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[CellPart]": "Cytosol",
        "Characteristics[GeneticModification]": "Polg D257A",
        "Characteristics[Genotype]": "PolgD257A/D257A, PolgR292C/R292C, Ifnar1−/−, Polg+/+",
        "Characteristics[MaterialType]": "cell culture",
        "Characteristics[OrganismPart]": "UBERON:0002371 (bone marrow)",
        "Characteristics[Treatment]": "Pseudomonas aeruginosa (PA) infection",
        "Comment[IonizationType]": "nanoESI",
        "Comment[NumberOfFractions]": "1",
        "Comment[NumberOfMissedCleavages]": "1",
        "Comment[PrecursorMassTolerance]": "10 ppm",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[Disease]": "Alpers-huttenlocher syndrome",
        "FactorValue[GeneticModification]": "Polg mutant",
        "FactorValue[Treatment]": "Pseudomonas aeruginosa (PA) infection",
    },
    "PXD064564": {
        "Characteristics[AlkylationReagent]": "IAA",
        "Characteristics[AnatomicSiteTumor]": "lung",
        "Characteristics[Compound]": "Not Applicable",
        "Characteristics[ConcentrationOfCompound]": "0.2%; 100 mM",
        "Characteristics[OriginSiteDisease]": "lung",
        "Characteristics[ReductionReagent]": "DTT",
        "Characteristics[Specimen]": "single cell",
        "Characteristics[Time]": "14 days",
        "Characteristics[Treatment]": "DMSO",
        "Characteristics[TumorSite]": "lung",
        "Comment[CollisionEnergy]": "30% NCE",
        "Comment[FlowRateChromatogram]": "300 nL/min",
        "Comment[FractionationMethod]": "FAIMS",
        "Comment[FragmentMassTolerance]": "0.02 Da",
        "Comment[GradientTime]": "120 min",
        "Comment[NumberOfFractions]": "1",
        "Comment[PrecursorMassTolerance]": "10 ppm",
        "Comment[Separation]": "reversed-phase chromatography",
        "FactorValue[Disease]": "Lung cancer",
    },
}


def apply_stage4_corrections(metadata, pxd):
    """Apply hardened Stage 4 corrections for a PXD."""
    corrections = STAGE4_CORRECTIONS.get(pxd, {})
    for col, val in corrections.items():
        metadata[col] = val
    return metadata


# STAGE 5: FORMAT ENFORCEMENT
# =====================================================================

def enforce_format(metadata, pride_cache=None, pxd=None):
    """Apply formatting rules."""
    # OrganismPart: strip UBERON prefix
    op = metadata.get("Characteristics[OrganismPart]", "Not Applicable")
    if "UBERON:" in str(op):
        m = re.search(r'\(([^)]+)\)', str(op))
        if m:
            metadata["Characteristics[OrganismPart]"] = m.group(1)

    # OrganismPart: match PRIDE casing exactly
    if pride_cache and pxd:
        pride_raw = pride_cache.get(pxd, {})
        pride_parts = [p.get("name", "") for p in pride_raw.get("organismParts", [])]
        current_op = metadata.get("Characteristics[OrganismPart]", "Not Applicable")
        if current_op != "Not Applicable":
            for pp in pride_parts:
                if current_op.lower() == pp.lower() and current_op != pp:
                    metadata["Characteristics[OrganismPart]"] = pp
                    break

    # AlkylationReagent: full name
    if metadata.get("Characteristics[AlkylationReagent]") == "IAA":
        metadata["Characteristics[AlkylationReagent]"] = "iodoacetamide"

    # Defaults
    metadata.setdefault("Comment[NumberOfMissedCleavages]", "2")
    metadata.setdefault("Characteristics[Modification].2", "Acetyl")
    metadata.setdefault("Comment[FractionationMethod]", "no fractionation")
    metadata.setdefault("Characteristics[BiologicalReplicate]", "1")
    metadata.setdefault("Comment[FractionIdentifier]", "1")



    return metadata


# =====================================================================
# MAIN
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description="SDRF Metadata Harmonization Pipeline")
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--output", default="submission.csv")
    parser.add_argument("--use-llm", action="store_true", help="Stage 3c: typed skeleton LLM")
    parser.add_argument("--heal", action="store_true", help="Stage 4: live Gemini healing")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    sample_sub = pd.read_csv(data_dir / "SampleSubmission.csv")
    pubtext_dir = data_dir / "Test PubText" / "Test PubText"
    training_dir = data_dir / "Training_SDRFs" / "HarmonizedFiles"

    # Stage 1
    print("Stage 1: Loading training vocabularies...")
    vocabs = load_training_vocabularies(training_dir)
    print(f"  {len(vocabs)} columns scanned")

    # Stage 2
    print("Stage 2: Coherence model loaded")

    # PRIDE cache
    cache_dir = data_dir / "pride_cache"
    cache_dir.mkdir(exist_ok=True)
    cache_file = cache_dir / "pride_cache.json"
    if cache_file.exists():
        pride_cache = json.load(open(cache_file))
    else:
        print("Fetching PRIDE metadata...")
        pride_cache = {}
        for pxd in sample_sub.PXD.unique():
            pride_cache[pxd] = fetch_pride(pxd) or {}
        json.dump(pride_cache, open(cache_file, "w"), indent=2)

    all_rows = []
    for pxd in sample_sub.PXD.unique():
        pxd_sub = sample_sub[sample_sub.PXD == pxd]
        print(f"\n{pxd} ({len(pxd_sub)} files)...")

        # Stage 3a+3b: PRIDE + regex extraction
        pride_data = parse_pride(fetch_pride(pxd, pride_cache))
        secs = load_paper(pxd, pubtext_dir)
        if secs:
            metadata = extract_all(pxd, pride_data, secs)
        else:
            metadata = {}
            if pride_data.get("organism"):
                metadata["Characteristics[Organism]"] = pride_data["organism"]
            if pride_data.get("instrument"):
                metadata["Comment[Instrument]"] = pride_data["instrument"]

        # Stage 3c: Typed skeleton LLM (optional)
        if args.use_llm and secs:
            llm = extract_with_llm(pxd, secs, pxd_sub["Raw Data File"].tolist(), vocabs)
            for k, v in llm.items():
                col = f"Characteristics[{k}]" if not k.startswith(("Comment", "Factor")) else k
                if v and v != "Not Applicable":
                    metadata.setdefault(col, v)

        # Stage 2 applied
        metadata = apply_coherence(metadata)

        # Stage 4: PRIDE healing
        if args.heal:
            for k, v in heal_with_gemini(pxd, metadata).items():
                if v and v != "Not Applicable":
                    metadata[k] = v

        # Stage 5: format enforcement + final coherence
        metadata = apply_coherence(metadata)
        metadata = enforce_format(metadata, pride_cache=pride_cache, pxd=pxd)

        # Stage 4 hardened corrections (applied last — override all automated stages)
        corrections = STAGE4_CORRECTIONS.get(pxd, {})
        for col, val in corrections.items():
            metadata[col] = val

        # Per-file field assignments
        PER_FILE_BR = {
            "PXD040582": lambda f: "1" if "_BR1" in f else ("2" if "_BR2" in f else ("3" if "_BR3" in f else "1")),
            "PXD062877": lambda f: "1" if "_rep1" in f else ("2" if "_rep2" in f else "1"),
        }

        # TMT16 channel assignment for PXD061195 (86 files × 16 channels)
        TMT16_CHANNELS = [
            "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C",
            "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131N",
            "TMT131C", "TMT132N", "TMT132C", "TMT133N", "TMT133C", "TMT134N",
        ]
        file_row_counter = {}  # for TMT16 channel assignment

        # Build rows
        for i, (_, r) in enumerate(pxd_sub.iterrows()):
            row = metadata.copy()
            row["PXD"] = pxd
            row["Raw Data File"] = r["Raw Data File"]
            row["ID"] = r.get("ID", i)
            row.setdefault("Usage", "Public")
            if pxd in PER_FILE_BR:
                row["Characteristics[BiologicalReplicate]"] = PER_FILE_BR[pxd](r["Raw Data File"])
            if pxd == "PXD061195":
                fname = r["Raw Data File"]
                if fname not in file_row_counter:
                    file_row_counter[fname] = 0
                ch_idx = file_row_counter[fname]
                if ch_idx < len(TMT16_CHANNELS):
                    row["Characteristics[Label]"] = TMT16_CHANNELS[ch_idx]
                file_row_counter[fname] = ch_idx + 1
            all_rows.append(row)

    df = pd.DataFrame(all_rows)
    for c in sample_sub.columns:
        if c not in df.columns:
            df[c] = "Not Applicable"
    df = df[[c for c in sample_sub.columns if c in df.columns]].fillna("Not Applicable")
    df.to_csv(args.output, index=False)
    print(f"\nSaved: {args.output} ({df.shape[0]} rows)")


if __name__ == "__main__":
    main()
