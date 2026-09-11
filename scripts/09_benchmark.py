#!/usr/bin/env python3
"""Stage 9 - run the pipeline over the truth sets and score it.

Reports per-category recall and precision, and every disagreement with its
expected and observed value. A benchmark that only flatters is worthless,
so disagreements are listed in full rather than summarised away.

    python scripts/09_benchmark.py --ref ref/ --outdir benchmarks/results

Review status matters. ClinVar rows carry the submitter's assertion level,
and 73 of the 336 small-variant rows are "no assertion criteria provided" -
a submission with no stated evidence behind it. Scoring against those moves
the figure for reasons that have nothing to do with this tool, so they are
excluded by default and reported separately. `--min-review-status any`
scores everything.

Matching is on locus plus protein change, not on transcript id: the truth
set records RefSeq transcripts (NM_...) while the pipeline works in
GENCODE (ENST...), and mapping between them would introduce a second
source of error into the measurement. A truth row counts as recovered when
some output entry at that locus carries the expected protein change.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2p.provenance import RunLogger                       # noqa: E402

TRUTH = ROOT / "benchmarks" / "truth"

# ClinVar review status, weakest first.
REVIEW_RANK = {
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 2,
    "criteria provided, multiple submitters, no conflicts": 3,
    "reviewed by expert panel": 4,
    "practice guideline": 5,
}
LEVELS = {
    "any": -1,
    "asserted": 1,          # anything with stated criteria - the default
    "single": 2,
    "multiple": 3,
    "expert": 4,
    "guideline": 5,
}


def read_truth(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with open(path, encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def rank_of(row: dict) -> int:
    return REVIEW_RANK.get(
        (row.get("clinvar_review_status") or "").strip().lower(), 0)


def kv(header: str) -> dict[str, str]:
    out = {}
    for tok in header.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


_FS = re.compile(r"^([a-z])(\d+)(?:[a-z])?fs(?:\*\d+)?$")
# ClinVar `p.T946_L947ins*` - a stop inserted between two residues.
_INS_STOP = re.compile(r"^([a-z])(\d+)_([a-z])(\d+)ins\*$")
# v2p `p.L947*fs*0` or a plain `p.L947*` - that residue becomes a stop.
_STOP_AT = re.compile(r"^([a-z])(\d+)\*(?:fs\*\d+)?$")


def norm_change(p: str) -> str:
    """Normalise a protein change for comparison.

    Lower-cased, `p.` stripped, and `*`/`Ter`/`X` unified, because the two
    sources spell a stop three different ways and a spelling difference is
    not a disagreement worth reporting.
    """
    s = (p or "").strip()
    s = s[2:] if s.lower().startswith("p.") else s
    s = s.replace("Ter", "*").replace("ter", "*")
    if s.endswith("X"):
        s = s[:-1] + "*"
    return s.lower()


def equivalent(want: str, got: str) -> bool:
    """Whether two normalised protein changes describe the same event.

    Exact string equality is too strict to be honest here. ClinVar writes a
    frameshift as `p.M862fs` - first affected residue only - while v2p
    writes the fuller `p.M862Ifs*4`, naming the substituted residue and the
    distance to the new stop. Those are the same event described at
    different precision, and scoring the second as a miss would report a
    tool failure that did not happen.

    What is *not* forgiven: a different residue, a different position, or a
    different class of consequence.
    """
    if want == got:
        return True
    a, b = _FS.match(want), _FS.match(got)
    if a and b:
        return a.group(1) == b.group(1) and a.group(2) == b.group(2)

    # Synonymous, spelled two ways. ClinVar writes `p.G55=` naming the
    # unchanged residue; v2p writes `p.(=)` for the whole protein. Both
    # assert the same thing - no amino acid changed.
    def _syn(x: str) -> bool:
        return x == "(=)" or x.endswith("=")

    if _syn(want) and _syn(got):
        return True

    # A stop inserted between two residues terminates the protein at the
    # first of them, which is the same product as that residue becoming a
    # stop. ClinVar: `p.T946_L947ins*`. v2p: `p.L947*fs*0`. Same protein.
    ins = _INS_STOP.match(want) or _INS_STOP.match(got)
    other = got if _INS_STOP.match(want) else want
    if ins:
        end = _STOP_AT.match(other)
        if end and end.group(2) in (ins.group(2), ins.group(4)):
            return True
    return False


def write_inputs(folder: Path, small: list[dict], editing: list[dict]) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    if small:
        with open(folder / "truth_variants.vcf", "w", newline="\n") as fh:
            fh.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
            for c in sorted({r["chrom"] for r in small}):
                fh.write(f"##contig=<ID={c}>\n")
            fh.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">\n')
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for r in sorted(small, key=lambda x: (x["chrom"], int(x["pos"]))):
                fh.write(f"{r['chrom']}\t{r['pos']}\t.\t{r['ref']}\t"
                         f"{r['alt']}\t.\tPASS\tGENE={r['gene']}\n")
    if editing:
        cols = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
                "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene",
                "AAChange.refGene"]
        with open(folder / "truth_editing.txt", "w", newline="\n") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in editing:
                fh.write("\t".join([
                    r["chrom"], r["pos"], r["pos"], r["ref"], r["alt"],
                    "exonic", r["gene"], ".", "nonsynonymous SNV",
                    f"{r['gene']}:{r.get('transcript', '.')}:.:.:"
                    f"{r['expected_hgvs_p']}"]) + "\n")


def observed_changes(fasta: Path) -> dict[str, set[str]]:
    """locus -> the set of protein changes the pipeline produced there.

    Reads `tables/provenance.jsonl` when present, falling back to FASTA
    headers only when it is not. This is a correctness fix, not a
    preference: cross-class deduplication merges identical sequences and
    keeps a single header, so a locus whose protein duplicates another
    variant's vanishes from the headers. Scoring off headers alone
    under-counted recall - two HRAS variants both encoding p.F156L
    collapsed into one entry, and the benchmark called the second a miss
    when the protein had in fact been built.
    """
    prov = fasta.parent / "tables" / "provenance.jsonl"
    got: dict[str, set[str]] = defaultdict(set)
    if prov.is_file():
        with open(prov, encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                for v in json.loads(line).get("variants", []):
                    loc, pc = v.get("locus"), v.get("protein_change")
                    if not loc or not pc:
                        continue
                    m = re.match(r"^(chr[^:]+):(\d+)", loc)
                    if m:
                        got[f"{m.group(1)}:{m.group(2)}"].add(norm_change(pc))
        if got:
            return got

    with open(fasta, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            f = kv(line.rstrip("\n"))
            loc, pc = f.get("LOC"), f.get("PC")
            if not loc or not pc:
                continue
            m = re.match(r"^(chr[^:]+):(\d+)", loc)
            if m:
                got[f"{m.group(1)}:{m.group(2)}"].add(norm_change(pc))
    return got


def score(truth: list[dict], got: dict[str, set[str]], expect_key: str):
    """Return (per-class counters, list of disagreements)."""
    per: dict[str, Counter] = defaultdict(Counter)
    bad = []
    for r in truth:
        cls = r.get("consequence_class") or "editing"
        key = f"{r['chrom']}:{r['pos']}"
        want = norm_change(r.get(expect_key, ""))
        per[cls]["truth"] += 1
        seen = got.get(key)
        if not seen:
            per[cls]["no_entry"] += 1
            bad.append((r, want, ""))
        elif any(equivalent(want, g) for g in seen):
            per[cls]["hit"] += 1
        else:
            per[cls]["wrong"] += 1
            bad.append((r, want, ",".join(sorted(seen))[:60]))
    return per, bad


def table(per) -> list[str]:
    out = ["| category | truth | recovered | wrong change | no entry "
           "| recall | precision |",
           "|---|---:|---:|---:|---:|---:|---:|"]
    tot: Counter = Counter()
    for cls in sorted(per):
        c = per[cls]
        tot.update(c)
        attempted = c["hit"] + c["wrong"]
        rec = c["hit"] / c["truth"] if c["truth"] else 0.0
        pre = c["hit"] / attempted if attempted else 0.0
        out.append(f"| {cls} | {c['truth']} | {c['hit']} | {c['wrong']} | "
                   f"{c['no_entry']} | {rec:.1%} | {pre:.1%} |")
    att = tot["hit"] + tot["wrong"]
    out.append(f"| **all** | {tot['truth']} | {tot['hit']} | {tot['wrong']} | "
               f"{tot['no_entry']} | "
               f"{tot['hit'] / tot['truth'] if tot['truth'] else 0:.1%} | "
               f"{tot['hit'] / att if att else 0:.1%} |")
    return out


# v2p consequence -> the VEP term meaning the same thing. Milestone 7d:
# a second opinion, not a verdict. VEP and v2p read the same annotation
# but resolve transcripts independently, so a disagreement is a thing to
# look at rather than automatically a v2p bug.
VEP_EQUIV = {
    "missense": "missense_variant",
    "stop_gained": "stop_gained",
    "stop_lost": "stop_lost",
    "start_lost": "start_lost",
    "inframe_insertion": "inframe_insertion",
    "inframe_deletion": "inframe_deletion",
    "frameshift": "frameshift_variant",
    "synonymous": "synonymous_variant",
    "5_prime_UTR": "5_prime_UTR_variant",
    "3_prime_UTR": "3_prime_UTR_variant",
}


def vep_crosscheck(path: Path, fasta: Path):
    """Compare our consequence calls with VEP's, per locus.

    Returns (agree, disagree_rows, unmatched). A locus counts as agreeing
    when any v2p consequence at it maps to VEP's most severe term there.
    """
    ours: dict[str, set] = defaultdict(set)
    ours_tx: dict[str, dict[str, str]] = defaultdict(dict)
    with open(fasta, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            f = kv(line)
            loc, csq = f.get("LOC"), f.get("CSQ")
            if not loc or not csq:
                continue
            m = re.match(r"^(chr[^:]+):(\d+)", loc)
            if m:
                k = f"{m.group(1)}:{m.group(2)}"
                ours[k].add(csq)
                tx = (f.get("TX") or "").split(".")[0]
                if tx:
                    ours_tx[k][tx] = csq

    agree, disagree, unmatched = 0, [], 0
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = f"{row['chrom']}:{row['pos']}"
            theirs = (row.get("vep_consequence") or "").strip()
            mine = ours.get(key)
            if not mine or not theirs:
                unmatched += 1
                continue

            # Compare like with like. VEP's most severe consequence may sit
            # on a transcript we never translated, which manufactures a
            # disagreement out of a transcript-choice difference: MSH2
            # chr2:47471062 is missense on the MANE transcript - what
            # ClinVar and v2p both say - and stop_gained on another. When
            # the two tools share a transcript, score on that one.
            shared = None
            for pair in (row.get("vep_per_transcript") or "").split(";"):
                if ":" not in pair:
                    continue
                vtx, vcsq = pair.split(":", 1)
                vtx = vtx.split(".")[0]
                if vtx in ours_tx.get(key, {}):
                    if VEP_EQUIV.get(ours_tx[key][vtx]) == vcsq:
                        shared = True
                        break
                    shared = False if shared is None else shared
            if shared:
                agree += 1
                continue

            if any(VEP_EQUIV.get(c) == theirs for c in mine):
                agree += 1
            else:
                disagree.append((key, row.get("vep_gene", ""),
                                 ",".join(sorted(mine)), theirs))
    return agree, disagree, unmatched


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True,
                    help="directory holding genome, annotation, proteome")
    ap.add_argument("--outdir", default="benchmarks/results")
    ap.add_argument("--logdir", default="logs")
    ap.add_argument("--workdir", help="keep the intermediate run here")
    ap.add_argument("--vep-consequences", default="",
                    help="TSV from scripts/10_vep_annotate.py. Adds an "
                         "independent consequence cross-check ("
                         "milestone 7d). VEP is a second opinion, not "
                         "ground truth.")
    ap.add_argument("--min-review-status", default="asserted",
                    choices=sorted(LEVELS),
                    help="minimum ClinVar review status to score. Default "
                         "'asserted' excludes rows with no assertion "
                         "criteria; 'any' scores every row.")
    a = ap.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rl = RunLogger("09_benchmark", a.logdir)
    rl.add_params(min_review_status=a.min_review_status, ref=a.ref)

    small_all = read_truth(TRUTH / "cosmic_variants.tsv")
    editing = read_truth(TRUTH / "editing_sites.tsv")
    if not small_all and not editing:
        print(f"no truth files in {TRUTH}", file=sys.stderr)
        rl.close(status="error")
        return 2

    floor = LEVELS[a.min_review_status]
    small = [r for r in small_all if rank_of(r) >= floor]
    excluded = len(small_all) - len(small)
    rl.count("truth.small_total", len(small_all))
    rl.count("truth.small_scored", len(small))
    rl.count("truth.editing", len(editing))
    print(f"small variants: {len(small)} of {len(small_all)} scored "
          f"({excluded} below --min-review-status {a.min_review_status})")
    print(f"editing sites:  {len(editing)}")

    work = Path(a.workdir) if a.workdir else Path(tempfile.mkdtemp())
    inputs = work / "inputs"
    write_inputs(inputs, small_all, editing)   # every row is run once

    print("running the pipeline over the truth set ...")
    cmd = [sys.executable, str(ROOT / "scripts" / "v2p"), "run", str(inputs),
           "--ref", a.ref, "--outdir", str(work / "release"),
           "--name", "TRUTH", "--logdir", str(work / "logs"),
           "--decoys", "none", "--no-reference"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    fastas = sorted((work / "release").glob("*.target.fasta"))
    if not fastas:
        sys.stderr.write(p.stdout[-3000:] + "\n" + p.stderr[-3000:] + "\n")
        rl.close(status="error")
        print("the pipeline produced no FASTA; see the output above",
              file=sys.stderr)
        return 1

    got = observed_changes(fastas[0])
    per_s, bad_s = score(small, got, "expected_hgvs_p_1letter")
    per_e, bad_e = score(editing, got, "expected_hgvs_p")
    per_full, _ = score(small_all, got, "expected_hgvs_p_1letter")

    L = ["# v2p benchmark", "",
         f"Scored with `--min-review-status {a.min_review_status}`: "
         f"{len(small)} of {len(small_all)} small-variant rows "
         f"({excluded} excluded).", "",
         "Matching is on locus plus protein change. The truth set records "
         "RefSeq transcripts and the pipeline works in GENCODE, so a row "
         "counts as recovered when some entry at that locus carries the "
         "expected change.", "",
         "## Small variants (headline)", ""]
    L += table(per_s)
    L += ["", "## Small variants, all rows including unasserted (secondary)",
          "", f"All {len(small_all)} rows, including the {excluded} with no "
          "stated assertion criteria. Reported for completeness; the "
          "headline figure above is the defensible one.", ""]
    L += table(per_full)
    if editing:
        L += ["", "## A-to-I RNA editing", "",
              "No competing tool accepts an editing table, so this category "
              "has no comparator.", ""]
        L += table(per_e)

    if a.vep_consequences and Path(a.vep_consequences).is_file():
        ag, dis, un = vep_crosscheck(Path(a.vep_consequences), fastas[0])
        tot = ag + len(dis)
        rl.count("vep.agree", ag)
        rl.count("vep.disagree", len(dis))
        L += ["", "## Consequence cross-check against Ensembl VEP", "",
              f"{ag} of {tot} loci agree "
              f"({ag / tot if tot else 0:.1%}); {un} could not be compared "
              f"(no protein produced, or VEP returned no coding "
              f"consequence).", "",
              "VEP is an independent second opinion, not ground truth. Both "
              "read the same annotation but resolve transcripts "
              "independently, so a disagreement is a thing to look at.", ""]
        if dis:
            L += ["| locus | gene | v2p | VEP |", "|---|---|---|---|"]
            L += [f"| {k} | {g} | {m} | {t} |" for k, g, m, t in dis[:40]]
            if len(dis) > 40:
                L.append(f"| ... | | {len(dis) - 40} more | |")
        print(f"VEP cross-check: {ag}/{tot} agree "
              f"({ag / tot if tot else 0:.1%}), {len(dis)} disagree")

    L += ["", "## Disagreements", ""]
    if not (bad_s + bad_e):
        L.append("None.")
    else:
        L += ["| locus | gene | expected | observed |", "|---|---|---|---|"]
        for r, want, seen in (bad_s + bad_e):
            L.append(f"| {r['chrom']}:{r['pos']} | {r.get('gene', '')} | "
                     f"{want} | {seen or '(no protein produced)'} |")

    rep = out / "benchmark.md"
    rep.write_text("\n".join(L) + "\n", encoding="utf-8")
    rl.add_output("report", rep)

    print()
    print("\n".join(table(per_s)))
    hit = sum(per_s[c]["hit"] for c in per_s)
    n = sum(per_s[c]["truth"] for c in per_s)
    print(f"\nsmall-variant recall (headline): {hit / n if n else 0:.1%}")
    print(f"disagreements: {len(bad_s) + len(bad_e)} (all listed in {rep})")
    rl.count("small.hit", hit)
    rl.count("small.truth", n)
    rl.close()
    if not a.workdir:
        print(f"intermediate run kept at {work}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
