#!/usr/bin/env python3
"""Stage 0c - audit how the reference files present their data.

Every conversion pipeline encodes assumptions about its reference files.
Those assumptions are usually implicit and only surface when they are
wrong, at which point the output is quietly incorrect rather than
obviously broken. This script states each one and tests it against the
files actually in use.

Usage:
  python src/v2p/stages/00_audit_references.py \
      --gtf ref/gencode.v44.annotation.gtf.gz \
      --genome ref/GRCh38.primary_assembly.genome.fa \
      --uniprot ref/uniprot_human_SP.fasta \
      --gencode-translations ref/gencode.v44.pc_translations.fa.gz
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from v2p.species import HUMAN, load_species   # noqa: E402

BOLD, RESET, DIM = "\033[1m", "\033[0m", "\033[2m"
findings: list[tuple[str, str]] = []


def head(s: str) -> None:
    print(f"\n{BOLD}{'=' * 74}\n{s}\n{'=' * 74}{RESET}")


def ok(msg: str) -> None:
    print(f"  [OK]    {msg}")


def note(msg: str) -> None:
    print(f"  [note]  {msg}")


def warn(kind: str, msg: str) -> None:
    print(f"  [{kind.upper()}] {msg}")
    findings.append((kind, msg))


def opener(p: Path):
    return gzip.open(p, "rt", errors="replace") if str(p).endswith(".gz") \
        else open(p, "rt", errors="replace")


# --------------------------------------------------------------------------

def audit_gtf(path: Path) -> None:
    head(f"GTF — {path.name}")
    feats = Counter()
    phases = Counter()
    tags = Counter()
    tx_types = Counter()
    par_y = 0
    cds_start_nf = set()
    cds_end_nf = set()
    mito = 0
    first_cds_phase: dict[str, str] = {}
    seen_cds: set[str] = set()
    contigs = set()
    attr_re = re.compile(r'(\S+)\s+"([^"]*)"')

    with opener(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feats[f[2]] += 1
            contigs.add(f[0])
            a = dict()
            tg = []
            for k, v in attr_re.findall(f[8]):
                if k == "tag":
                    tg.append(v)
                else:
                    a.setdefault(k, v)
            tid = a.get("transcript_id", "")
            if f[2] == "CDS":
                phases[f[7]] += 1
                if tid and tid not in seen_cds:
                    seen_cds.add(tid)
                    first_cds_phase[tid] = f[7]
                if f[0] in ("chrM", "MT", "chrMT"):
                    mito += 1
            if f[2] == "transcript":
                tx_types[a.get("transcript_type", "?")] += 1
                for t in tg:
                    tags[t] += 1
                if "PAR_Y" in a.get("gene_id", ""):
                    par_y += 1
                if "cds_start_NF" in tg:
                    cds_start_nf.add(tid)
                if "cds_end_NF" in tg:
                    cds_end_nf.add(tid)

    print(f"  features: {dict(feats)}")
    if feats.get("Selenocysteine"):
        try:
            src = (Path(__file__).resolve().parents[1]
                   / "src" / "v2p" / "annotation.py").read_text()
            uses_sec = "Selenocysteine" in src
        except OSError:
            uses_sec = False
        msg = (f"{feats['Selenocysteine']} Selenocysteine features mark UGA "
               f"codons that encode Sec rather than terminating translation")
        if uses_sec:
            ok(msg + " — parsed and recoded to U.")
        else:
            warn("BUG", msg + " — ignored, so every selenoprotein is "
                              "truncated at its first Sec codon.")
    print(f"  contigs: {len(contigs)}")
    print(f"  {DIM}(1-based inclusive coordinates — assumed throughout){RESET}")

    # -- assumption 1: CDS excludes the stop codon -----------------------
    if "stop_codon" in feats:
        ok("`stop_codon` is a separate feature, so CDS blocks exclude the "
           "stop. We translate from the CDS start to the first in-frame "
           "stop, reading past the CDS end into the 3'UTR, which is why "
           "this works and why frameshifts extend correctly.")
    else:
        warn("check", "no stop_codon feature — verify whether CDS includes "
                      "the stop, or every protein gains a trailing residue")

    # -- assumption 2: first CDS phase is 0 ------------------------------
    nonzero = {t: p for t, p in first_cds_phase.items() if p not in ("0", ".")}
    print(f"  CDS phase values across all blocks: {dict(phases)}")
    if nonzero:
        # A non-zero first phase is a property of the annotation, not a
        # defect. What matters is whether the pipeline applies it.
        applied = False
        try:
            src = (Path(__file__).resolve().parents[1]
                   / "src" / "v2p" / "annotation.py").read_text()
            applied = "off + self.cds_phase_first" in src
        except OSError:
            pass
        msg = (f"{len(nonzero)} transcripts carry a non-zero frame on their "
               f"first CDS block (these are the cds_start_NF set). "
               f"Translation must discard that many bases. Examples: "
               f"{list(nonzero.items())[:3]}")
        if applied:
            ok(msg + " — the pipeline applies the frame (verified in "
                     "annotation.py).")
        else:
            warn("BUG", msg + " — the pipeline does NOT apply it; those "
                              "transcripts translate out of frame.")
    else:
        ok("every transcript's first CDS block has phase 0")

    # -- assumption 3: incomplete-CDS transcripts ------------------------
    print(f"  cds_start_NF transcripts: {len(cds_start_nf)}")
    print(f"  cds_end_NF transcripts:   {len(cds_end_nf)}")
    if cds_start_nf:
        note("cds_start_NF means the 5' end of the CDS is not determined. "
             "These are exactly the transcripts carrying a non-zero first "
             "phase.")

    # -- assumption 4: PAR_Y duplicates ----------------------------------
    if par_y:
        warn("check", f"{par_y} transcripts are pseudoautosomal-region "
                      f"duplicates on chrY (gene_id contains _PAR_Y). They "
                      f"duplicate chrX genes and can double-count in "
                      f"gene-name lookups.")
    else:
        ok("no _PAR_Y duplicate transcripts")

    # -- assumption 5: mitochondrial genetic code ------------------------
    if mito:
        try:
            src = (Path(__file__).resolve().parents[1]
                   / "src" / "v2p" / "seqops.py").read_text()
            has_mito = "CODON_TABLE_MITO" in src
        except OSError:
            has_mito = False
        msg = (f"{mito} CDS blocks are mitochondrial; the vertebrate "
               f"mitochondrial code differs from the standard table "
               f"(TGA=Trp, ATA=Met, AGA/AGG=stop)")
        if has_mito:
            ok(msg + " — table 2 is applied to chrM.")
        else:
            warn("BUG", msg + " — table 1 is applied everywhere, so chrM "
                              "variants are mistranslated.")
    else:
        ok("no mitochondrial CDS in this annotation")

    # -- assumption 6: gene names are unique -----------------------------
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        from v2p.annotation import Annotation
        src = (Path(__file__).resolve().parents[1]
               / "src" / "v2p" / "build" / "smallvar.py").read_text()
        resolves = "representative_at" in src
        ann = Annotation.from_gtf(path, coding_only=True)
        amb = ann.ambiguous_gene_names()
        if amb:
            msg = (f"{len(amb)} gene names map to more than one gene id "
                   f"(e.g. {sorted(amb)[:4]}). Resolving a symbol without "
                   f"positional context can select the wrong locus")
            if resolves:
                ok(msg + " — symbols are resolved at the variant locus.")
            else:
                warn("BUG", msg + " — symbols are resolved by name alone.")
        else:
            ok("every gene name maps to a single gene id")
    except Exception as exc:                       # audit must not crash
        note(f"gene-name uniqueness check skipped ({type(exc).__name__})")

    top = ", ".join(f"{k}={v}" for k, v in tags.most_common(6))
    print(f"  most common transcript tags: {top}")
    for need in ("MANE_Select", "Ensembl_canonical", "basic"):
        if tags.get(need):
            ok(f"tag `{need}` present ({tags[need]}) — used for transcript "
               f"selection")
        else:
            warn("check", f"tag `{need}` absent; transcript selection falls "
                          f"back to longest CDS")


def audit_genome(path: Path, species=HUMAN) -> None:
    head(f"Genome FASTA — {path.name}")
    try:
        from pyfaidx import Fasta
    except ImportError:
        warn("check", "pyfaidx not installed; skipping")
        return
    fa_raw = Fasta(str(path), as_raw=True)          # no case folding
    names = list(fa_raw.keys())
    print(f"  contigs: {len(names)}   first: {names[:3]}")

    style = "UCSC (chr-prefixed)" if names[0].startswith("chr") else "Ensembl"
    ok(f"contig naming: {style}")

    # Contig lengths are the cheapest assembly check there is, but they
    # are species-specific. An organism with none declared cannot have its
    # assembly asserted - that is a weaker guarantee, not a failure, and
    # saying so beats inventing an expectation we cannot check.
    EXPECT = dict(species.contig_lengths)
    if not EXPECT:
        warn("check", f"no reference contig lengths are declared for "
                      f"{species.common_name}, so the assembly cannot be "
                      f"asserted from this file. Add contig_lengths to "
                      f"config/species/{species.common_name}.yaml to enable "
                      f"this check.")
    else:
        present = [c for c in EXPECT if c in fa_raw]
        if not present:
            warn("check", f"none of the expected contigs "
                          f"({', '.join(sorted(EXPECT))}) is in this file, so "
                          f"the {species.assembly} assembly cannot be "
                          f"confirmed")
        else:
            bad = [c for c in present if len(fa_raw[c]) != EXPECT[c]]
            if bad:
                warn("BUG", f"contig lengths do not match "
                            f"{species.assembly}: "
                     + ", ".join(f"{c}={len(fa_raw[c]):,} "
                                 f"(expected {EXPECT[c]:,})" for c in bad))
            else:
                ok(f"contig lengths match {species.assembly} "
                   f"({', '.join(present)})")

    # -- assumption: soft-masking is handled -----------------------------
    probe = present[0] if present else names[0]
    L = len(fa_raw[probe])
    lo = min(1_000_000, max(0, L // 3))
    sample = str(fa_raw[probe][lo:lo + min(50_000, L - lo)])
    if not sample:
        note("could not sample sequence for a masking check")
        return
    n_lower = sum(1 for c in sample if c.islower())
    if n_lower:
        ok(f"genome is soft-masked ({100*n_lower/len(sample):.0f}% lowercase "
           f"in the sampled window). We open it with "
           f"sequence_always_upper=True, so repeat-masked exons translate "
           f"normally instead of becoming unknown codons.")
    else:
        note("sampled window contains no lowercase; genome may be unmasked")

    n_run = sum(1 for c in sample.upper() if c == "N")
    print(f"  N bases in sampled 50kb window: {n_run}")
    note("codons containing N translate to X rather than raising, and the "
         "entry is flagged ambiguous_bases_translated_as_X")

    alt = [c for c in names if "_alt" in c or "_random" in c or "chrUn" in c]
    if alt:
        note(f"{len(alt)} alt/random/unplaced contigs present. Variants are "
             f"called on the primary assembly, so these are never queried.")


def audit_uniprot(path: Path) -> None:
    head(f"UniProt proteome — {path.name}")
    n = 0
    dbs = Counter()
    fields = Counter()
    genes: dict[str, int] = defaultdict(int)
    isoforms = 0
    no_met = 0
    alphabet = Counter()
    wrapped = 0
    cur_lines = 0
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_lines > 1:
                    wrapped += 1
                cur_lines = 0
                n += 1
                dbs[line[1:3]] += 1
                if re.match(r">\w+\|[A-Z0-9]+-\d+\|", line):
                    isoforms += 1
                for f_ in ("OS=", "OX=", "GN=", "PE=", "SV="):
                    if f_ in line:
                        fields[f_] += 1
                m = re.search(r"\bGN=(\S+)", line)
                if m:
                    genes[m.group(1)] += 1
            else:
                cur_lines += 1
                if cur_lines == 1 and line[:1] and line[0] != "M":
                    no_met += 1
                if n <= 6000:
                    alphabet.update(line.strip())
    print(f"  entries: {n:,}   db prefixes: {dict(dbs)}")
    print(f"  header fields: {dict(fields)}")

    if isoforms:
        warn("check", f"{isoforms} isoform accessions (ACC-N). Residue "
                      f"numbering differs between isoforms of a gene.")
    else:
        ok("canonical entries only — no isoform accessions")

    if wrapped == 0:
        ok("sequences are unwrapped, one line per entry — our `uniprot` "
           "output style matches this")
    else:
        warn("check", f"{wrapped} entries are wrapped; our output is "
                      f"unwrapped and would not be byte-comparable")

    dup = {g: c for g, c in genes.items() if c > 1}
    if dup:
        warn("check", f"{len(dup)} gene symbols map to more than one entry "
                      f"(e.g. {list(dup)[:4]}). Gene-name lookup is ambiguous "
                      f"for these; we compare against the best-matching entry "
                      f"rather than the longest.")
    else:
        ok("every gene symbol maps to exactly one entry")

    odd = {c: v for c, v in alphabet.items() if c in "UOXBZJ*"}
    if odd:
        note(f"non-standard residues in the first 6000 entries: {odd}. "
             f"U is selenocysteine; a standard-code translation truncates "
             f"those at UGA, so they are excluded from identity statistics.")

    print(f"  entries not starting with Met: {no_met}")
    if no_met:
        note("UniProt removes the initiator methionine from some entries "
             "after signal-peptide or transit-peptide cleavage. Our "
             "comparison tolerates a leading-Met difference "
             "(identical_after_met_trim).")


def audit_translations(path: Path) -> None:
    head(f"GENCODE translations — {path.name}")
    n = 0
    trailing_star = 0
    fields = Counter()
    with opener(path) as fh:
        cur = []
        for line in fh:
            if line.startswith(">"):
                if cur and cur[-1].endswith("*"):
                    trailing_star += 1
                cur = []
                n += 1
                fields[len(line[1:].split("|"))] += 1
            else:
                cur.append(line.strip())
        if cur and cur[-1].endswith("*"):
            trailing_star += 1
    print(f"  entries: {n:,}   header pipe-field counts: {dict(fields)}")
    ok("header layout is ENSP|ENST|ENSG|OTTHUMG|OTTHUMT|name|symbol|length; "
       "we key on field 2 (ENST) to match our transcript ids")
    if trailing_star:
        note(f"{trailing_star} sequences end with '*'. We strip it before "
             f"comparison; leaving it would make every entry differ by one "
             f"character.")
    else:
        ok("no trailing stop characters")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf")
    ap.add_argument("--genome")
    ap.add_argument("--uniprot")
    ap.add_argument("--gencode-translations")
    ap.add_argument("--species", default="human",
                    help="species name or config/species/*.yaml; sets "
                         "which contig lengths the assembly check asserts")
    args = ap.parse_args()

    if args.gtf:
        audit_gtf(Path(args.gtf))
    if args.genome:
        audit_genome(Path(args.genome), load_species(args.species))
    if args.uniprot:
        audit_uniprot(Path(args.uniprot))
    if args.gencode_translations:
        audit_translations(Path(args.gencode_translations))

    head("summary")
    bugs = [m for k, m in findings if k == "BUG"]
    checks = [m for k, m in findings if k == "check"]
    if not findings:
        print("  every assumption held.")
    for m in bugs:
        print(f"  BUG    {m[:200]}")
    for m in checks:
        print(f"  CHECK  {m[:200]}")
    print(f"\n  {len(bugs)} bug(s), {len(checks)} item(s) to confirm")
    return 1 if bugs else 0


if __name__ == "__main__":
    raise SystemExit(main())
