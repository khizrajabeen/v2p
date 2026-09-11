#!/usr/bin/env python3
"""Combinatorial proteoforms - co-occurring variants applied together.

Run:  python tests/test_combinatorial.py
Exit code 0 = all pass. No pytest dependency so it runs anywhere.

Three claims are load-bearing and each gets a positive and a negative
test:

1. Combining yields a protein carrying *both* changes, which neither
   single-variant entry contains.
2. The entry is emitted only when some tryptic peptide exists in the
   combined form and in neither single form - otherwise it is database
   bloat, and bloat costs FDR power.
3. Phase is honoured where the caller reports it. Two variants on
   opposite haplotypes are never combined, because no molecule carries
   both.
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

from make_fixture import build                              # noqa: E402
from v2p.annotation import Annotation, Genome               # noqa: E402
from v2p.build.combinatorial import (                       # noqa: E402
    COMBO_CLASS, CombinableVariant, _allele_changes,
    build_combinatorial_proteins, cooccurring_peptides, group_by_transcript,
    phase_groups,
)
from v2p.build.smallvar import build_small_variant_proteins  # noqa: E402
from v2p.fasta import write_fasta                           # noqa: E402
from v2p.parse.inputs import _info_af, _phase, parse_vcf    # noqa: E402
from v2p.seqops import table_for_contig, translate_orf      # noqa: E402

FIX = ROOT / "tests" / "fixtures"
PASS: list[str] = []
FAIL: list[str] = []

SWAP = {"A": "C", "C": "A", "G": "T", "T": "G"}

# Three clean tryptic fragments, each long enough to be a peptide.
REF_P = "MSSSSSSKLLLLLLLLKAAAAAAAAK"
NEAR_A = REF_P[:2] + "T" + REF_P[3:]          # residue 3, fragment 1
NEAR_B = REF_P[:5] + "T" + REF_P[6:]          # residue 6, fragment 1
FAR_B = REF_P[:19] + "G" + REF_P[20:]         # residue 20, fragment 3
COMBO_NEAR = NEAR_A[:5] + "T" + NEAR_A[6:]
COMBO_FAR = NEAR_A[:19] + "G" + NEAR_A[20:]


def check(name: str, cond: bool, detail: str = "") -> None:
    (PASS if cond else FAIL).append(name)
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name}{' | ' + detail if detail else ''}")


def gpos_plus(meta, gene: str, codon: int, base: int = 1) -> int:
    m = meta[gene] if gene in meta else next(
        v for k, v in meta.items() if k.split("@")[0] == gene)
    tx_off = m["utr5"] + 3 * (codon - 1) + (base - 1)
    seen = 0
    for s, e in m["exons"]:
        ln = e - s + 1
        if tx_off < seen + ln:
            return s + (tx_off - seen)
        seen += ln
    raise IndexError


def vcf_line(fmt: str, sample: str) -> list[str]:
    return ["chr1", "100", ".", "A", "G", ".", "PASS", ".", fmt, sample]


def test_peptide_necessity() -> None:
    """The filter that stops the database growing for nothing."""
    # Two variants in the same tryptic fragment: the fragment carrying
    # both exists in no single-variant form.
    both = cooccurring_peptides(COMBO_NEAR, [NEAR_A, NEAR_B], REF_P,
                                missed_cleavages=0)
    check("two variants in one peptide yield a peptide no single form has",
          bool(both), f"{sorted(both)}")

    # Two variants in different fragments, no missed cleavage: every
    # peptide of the combined form is already in one of the singles.
    apart = cooccurring_peptides(COMBO_FAR, [NEAR_A, FAR_B], REF_P,
                                 missed_cleavages=0)
    check("two variants in different peptides yield nothing new",
          apart == set(), f"{sorted(apart)}")

    # ... until missed cleavages let one peptide reach both.
    reach = cooccurring_peptides(COMBO_FAR, [NEAR_A, FAR_B], REF_P,
                                 missed_cleavages=2)
    check("missed cleavages can bring distant variants into one peptide",
          bool(reach), f"{len(reach)} peptide(s)")

    check("the reference's own peptides never count as new",
          cooccurring_peptides(REF_P, [], REF_P) == set())

    il = cooccurring_peptides(REF_P.replace("L", "I"), [], REF_P)
    check("an I/L-only difference is not new, because it is isobaric",
          il == set(), f"{sorted(il)}")


def test_allele_notation() -> None:
    """The entry describes the protein it contains, not its ingredients."""
    check("each changed residue is listed",
          _allele_changes("MAAA", "MVAG") == ["A2V", "A4G"],
          str(_allele_changes("MAAA", "MVAG")))
    check("a combination that changes no residue says so",
          _allele_changes("MAAA", "MAAA") == ["="])
    check("a truncation is reported as a stop, not a substitution",
          _allele_changes("MAAAA", "MA")[0].endswith("*"),
          str(_allele_changes("MAAAA", "MA")))
    long = _allele_changes("A" * 20, "C" * 20, limit=3)
    check("an implausibly long change list is capped",
          len(long) == 4 and long[-1] == "+17_more", str(long[-1]))


def test_allele_frequency() -> None:
    """A frequency filter must not act on silence."""
    check("AF is read from INFO", _info_af("DP=50;AF=0.25;MQ=60") == 0.25)
    check("VAF is accepted too, since somatic callers write it",
          _info_af("DP=50;VAF=0.98") == 0.98)
    check("a multi-allelic AF takes the first value",
          _info_af("AF=0.2,0.8") == 0.2)
    check("INFO without a frequency says nothing, not zero",
          _info_af("DP=50;MQ=60") is None)
    check("an empty INFO says nothing", _info_af("") is None)
    check("a malformed AF is not guessed at", _info_af("AF=high") is None)


def test_phase_grouping() -> None:
    """Phase is a fact when reported; its absence is a hypothesis."""
    def v(pos, ps="", haps=()):
        return CombinableVariant("chr1", pos, "A", "G", "SNV", "vcf", "",
                                 ps, frozenset(haps))

    same = phase_groups([v(10, "1", [0]), v(20, "1", [0])])
    check("two variants on the same haplotype form one phased group",
          len(same) == 1 and same[0][1] and len(same[0][2]) == 2,
          str([(g[0], len(g[2])) for g in same]))

    opp = phase_groups([v(10, "1", [0]), v(20, "1", [1])])
    check("two variants on opposite haplotypes are never combined",
          opp == [], str(opp))

    hom = phase_groups([v(10, "1", [0, 1]), v(20, "1", [1])])
    check("a homozygous variant joins both haplotypes",
          len(hom) == 1 and hom[0][0] == "1|1", str([g[0] for g in hom]))

    # Two variants each phased, but in different blocks: each block is
    # internally phased, their relative phase is not known. That is a
    # hypothesis, not a refusal.
    blocks = phase_groups([v(10, "1", [0]), v(20, "2", [0])])
    check("different phase sets are an unphased hypothesis, not a haplotype",
          len(blocks) == 1 and blocks[0][0] == "unphased"
          and not blocks[0][1], str([(g[0], g[1]) for g in blocks]))

    mixed = phase_groups([v(10, "1", [0]), v(20, "1", [0]), v(30)])
    labels = sorted(g[0] for g in mixed)
    check("a variant without phase makes an unphased hypothesis as well",
          labels == ["1|0", "unphased"], str(labels))

    plain = phase_groups([v(10), v(20)])
    check("with no phase at all there is one unphased group",
          len(plain) == 1 and not plain[0][1] and len(plain[0][2]) == 2)


def test_vcf_phase_parsing() -> None:
    """GT/PS reach the builder; a slashed GT is not phase."""
    check("a phased GT with a PS tag is phase",
          _phase(vcf_line("GT:PS", "0|1:12345"), 1) == ("12345", {1}),
          str(_phase(vcf_line("GT:PS", "0|1:12345"), 1)))
    check("an unphased GT is a genotype, not phase",
          _phase(vcf_line("GT", "0/1"), 1) == ("", frozenset()))
    check("a homozygous ALT sits on both haplotypes",
          _phase(vcf_line("GT:PS", "1|1:7"), 1) == ("7", {0, 1}))
    check("a phased GT without PS is one block",
          _phase(vcf_line("GT", "1|0"), 1) == ("*", {0}))
    check("the right ALT of a multi-allelic row is picked",
          _phase(vcf_line("GT", "0|2"), 2) == ("*", {1})
          and _phase(vcf_line("GT", "0|2"), 1) == ("", frozenset()))
    check("a sites-only VCF has no phase",
          _phase(["chr1", "1", ".", "A", "G", ".", "PASS", "."], 1)
          == ("", frozenset()))


VCF_HEAD = ("##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")


def write_vcf(path: Path, rows) -> Path:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(VCF_HEAD)
        for chrom, pos, ref, alt, fmt, sample in rows:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\t"
                     f"{fmt}\t{sample}\n")
    return path


def from_vcf(path: Path) -> list[CombinableVariant]:
    """Parse as the pipeline does, including the phase fields."""
    out = []
    for row in parse_vcf(path):
        p = row["payload"]
        out.append(CombinableVariant(
            p["chrom"], p["pos"], p["ref"], p["alt"], row["variant_class"],
            row["source"], row.get("confidence", ""),
            p.get("phase_set", ""), frozenset(p.get("haplotypes") or ())))
    return out


def test_phased_vcf_to_fasta(meta, ann, gen) -> None:
    """The three cases, from a real VCF through to a written FASTA.

    The opposite-haplotype case is asserted on the file's contents: no
    counter, no return value, the protein must simply not be there.
    """
    tmp = Path(tempfile.mkdtemp(prefix="v2p_phase_"))
    p3 = gpos_plus(meta, "GPLUS", 3, 1)
    p9 = gpos_plus(meta, "GPLUS", 9, 1)
    r3 = gen.fetch("chrT1", p3, p3)
    r9 = gen.fetch("chrT1", p9, p9)
    a3, a9 = SWAP[r3], SWAP[r9]

    def build(rows, name, **kw):
        vs = from_vcf(write_vcf(tmp / f"{name}.vcf", rows))
        recs = build_combinatorial_proteins(vs, ann, gen,
                                            transcript_mode="all", **kw)
        fa = tmp / f"{name}.fasta"
        write_fasta(recs, fa, style="uniprot")
        return recs, fa.read_text(encoding="utf-8")

    # 1. Same haplotype, same PS block -> one molecule carries both.
    same, same_fa = build(
        [("chrT1", p3, r3, a3, "GT:PS", "0|1:100"),
         ("chrT1", p9, r9, a9, "GT:PS", "0|1:100")], "same")
    check("a same-haplotype pair from a phased VCF combines", bool(same),
          f"{len(same)} record(s)")
    check("its phase is recorded in the FASTA header",
          "PHASE=phased" in same_fa,
          next((ln[-40:] for ln in same_fa.splitlines()
                if ln.startswith(">")), ""))
    combined_seq = same[0].sequence if same else ""

    # 2. Opposite haplotypes, same PS block -> no molecule carries both.
    opp, opp_fa = build(
        [("chrT1", p3, r3, a3, "GT:PS", "0|1:100"),
         ("chrT1", p9, r9, a9, "GT:PS", "1|0:100")], "opp")
    body = "".join(ln for ln in opp_fa.splitlines()
                   if not ln.startswith(">"))
    check("the opposite-haplotype protein is absent from the FASTA",
          combined_seq != "" and combined_seq not in body,
          f"fasta {len(opp_fa)} bytes, {len(opp)} record(s)")

    # 3. Different PS blocks -> relative phase unknown: a hypothesis.
    diff, diff_fa = build(
        [("chrT1", p3, r3, a3, "GT:PS", "0|1:100"),
         ("chrT1", p9, r9, a9, "GT:PS", "0|1:200")], "diff")
    check("variants in different PS blocks combine as a hypothesis",
          bool(diff) and "PHASE=unphased" in diff_fa,
          str(diff[0].notes) if diff else "none")

    diff_off, diff_off_fa = build(
        [("chrT1", p3, r3, a3, "GT:PS", "0|1:100"),
         ("chrT1", p9, r9, a9, "GT:PS", "0|1:200")], "diff_off",
        allow_unphased=False)
    check("--no-allow-unphased removes it from the FASTA",
          diff_off == [] and ">" not in diff_off_fa,
          f"{len(diff_off)} record(s)")

    # The AF filter, through a real VCF: the payload must carry the
    # frequency, and silence must survive it.
    rows = list(parse_vcf(write_vcf(tmp / "af.vcf", [
        ("chrT1", p3, r3, a3, "GT", "0/1")])))
    check("a variant with no INFO frequency is not assigned one",
          rows[0]["payload"]["af"] is None, str(rows[0]["payload"]["af"]))

    with open(tmp / "af2.vcf", "w", encoding="utf-8") as fh:
        fh.write(VCF_HEAD)
        fh.write(f"chrT1\t{p3}\t.\t{r3}\t{a3}\t.\tPASS\tDP=40;AF=0.12\t"
                 f"GT\t0/1\n")
    rows2 = list(parse_vcf(tmp / "af2.vcf"))
    check("a stated frequency reaches the payload",
          rows2[0]["payload"]["af"] == 0.12,
          str(rows2[0]["payload"]["af"]))

    # An unphased pair is still admitted by default: that is the common
    # case, and HCC1395's own VCF is sites-only.
    plain, _ = build([("chrT1", p3, r3, a3, "GT", "0/1"),
                      ("chrT1", p9, r9, a9, "GT", "0/1")], "plain")
    check("a slashed GT still combines by default, marked unphased",
          bool(plain) and plain[0].extra["phase"] == "unphased",
          str(plain[0].extra["phase"]) if plain else "none")

    shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    test_peptide_necessity()
    test_allele_notation()
    test_allele_frequency()
    test_phase_grouping()
    test_vcf_phase_parsing()

    meta = build()
    ann = Annotation.from_gtf(FIX / "mini.gtf")
    gen = Genome(FIX / "mini.fa")
    test_phased_vcf_to_fasta(meta, ann, gen)

    # Two substitutions on one transcript, several codons apart.
    p3 = gpos_plus(meta, "GPLUS", 3, 1)
    p9 = gpos_plus(meta, "GPLUS", 9, 1)
    r3 = gen.fetch("chrT1", p3, p3)
    r9 = gen.fetch("chrT1", p9, p9)
    a3, a9 = SWAP[r3], SWAP[r9]

    v3 = CombinableVariant("chrT1", p3, r3, a3, "SNV", "vcf")
    v9 = CombinableVariant("chrT1", p9, r9, a9, "SNV", "vcf")

    # What the single-variant path produces, for comparison.
    s3 = build_small_variant_proteins("chrT1", p3, r3, a3, ann, gen,
                                      variant_class="SNV", genes="GPLUS")
    s9 = build_small_variant_proteins("chrT1", p9, r9, a9, ann, gen,
                                      variant_class="SNV", genes="GPLUS")
    check("both single variants translate on their own",
          bool(s3) and bool(s9), f"{len(s3)}/{len(s9)}")

    # ------------------------------------------------------- grouping
    groups = group_by_transcript([v3, v9], ann, "all")
    check("two variants on one transcript are grouped",
          bool(groups) and all(len(val[1]) == 2 for val in groups.values()),
          f"{len(groups)} transcript(s)")
    check("a lone variant is not grouped",
          group_by_transcript([v3], ann, "all") == {})

    # ------------------------------------------------------- the point
    combo = build_combinatorial_proteins([v3, v9], ann, gen,
                                         transcript_mode="all")
    check("a combinatorial protein is produced", bool(combo),
          f"{len(combo)} record(s)")
    if not combo:
        print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
        return 1

    c = combo[0]
    single = {r.sequence for r in (s3 + s9)}
    check("the combined protein is not any single-variant protein",
          c.sequence not in single, f"len={len(c.sequence)}")

    only3 = s3[0].sequence
    only9 = s9[0].sequence
    diff3 = [i for i, (x, y) in enumerate(zip(only3, c.sequence)) if x != y]
    diff9 = [i for i, (x, y) in enumerate(zip(only9, c.sequence)) if x != y]
    check("the combined form differs from each single form at the other "
          "variant's residue", bool(diff3) and bool(diff9),
          f"vs-v3 {diff3[:3]}, vs-v9 {diff9[:3]}")

    check("its class is COMBO", c.variant_class == COMBO_CLASS,
          c.variant_class)
    check("both loci are recorded",
          c.extra["n_variants"] == 2 and len(c.extra["loci"]) == 2,
          str(c.extra.get("loci")))
    check("the peptides only the combination has are named",
          c.extra["n_cooccurring_peptides"] > 0
          and bool(c.extra["cooccurring_peptides"]),
          str(c.extra["cooccurring_peptides"][:1]))
    check("same-class variants are not flagged cross-evidence",
          c.extra["cross_evidence"] is False, str(c.extra["classes"]))

    # -------------------------------------------- phase, end to end
    check("without phase the entry is a hypothesis, in HGVS unphased form",
          "unphased" in c.notes and c.extra["phased"] is False
          and "(;)" in c.protein_change, c.protein_change)

    def phased(h3, h9):
        return build_combinatorial_proteins(
            [CombinableVariant("chrT1", p3, r3, a3, "SNV", "vcf", "", "77",
                               frozenset({h3})),
             CombinableVariant("chrT1", p9, r9, a9, "SNV", "vcf", "", "77",
                               frozenset({h9}))], ann, gen,
            transcript_mode="all")

    ph = phased(0, 0)
    check("phased input gives a phased entry in HGVS phased form",
          bool(ph) and ph[0].extra["phased"] is True
          and "phased" in ph[0].notes and "(;)" not in ph[0].protein_change,
          ph[0].protein_change if ph else "none")

    opp = phased(0, 1)
    check("variants on opposite haplotypes produce no protein", opp == [],
          f"{len(opp)} record(s)")

    # --------------------------------------------- two variants, one codon
    # Real HCC1395 calls put two variants in codon 225 of FKTN. Describing
    # them per variant claimed both D225N and D225E; the protein carries
    # a third residue that is neither.
    b1 = gpos_plus(meta, "GPLUS", 3, 1)
    b2 = gpos_plus(meta, "GPLUS", 3, 2)
    n1, n2 = gen.fetch("chrT1", b1, b1), gen.fetch("chrT1", b2, b2)
    codon = build_combinatorial_proteins(
        [CombinableVariant("chrT1", b1, n1, SWAP[n1], "SNV", "vcf"),
         CombinableVariant("chrT1", b2, n2, SWAP[n2], "SNV", "vcf")],
        ann, gen, transcript_mode="all")
    check("two variants in one codon are combined at all", bool(codon),
          f"{len(codon)} record(s)")
    if codon:
        m = (meta["GPLUS"] if "GPLUS" in meta else
             next(v for k, v in meta.items() if k.split("@")[0] == "GPLUS"))
        ref_prot = translate_orf(m["cds"], table_for_contig("chrT1")).protein
        pc = codon[0].protein_change
        aa3 = codon[0].sequence[2]
        check("one codon is described as one residue change, not two",
              pc.count(";") == 0 and pc == f"p.[{ref_prot[2]}3{aa3}]", pc)
        check("that residue is neither single-variant residue",
              aa3 != only3[2] and aa3 != ref_prot[2],
              f"combined {aa3}, single {only3[2]}, ref {ref_prot[2]}")

    # ------------------------------------ the case no other tool reaches
    # A somatic SNV and an ADAR edit on one transcript: a proteoform with
    # neither a purely genomic nor a purely transcriptomic basis.
    edit = CombinableVariant("chrT1", p9, r9, a9, "RNA_EDITING", "res_table")
    mixed = build_combinatorial_proteins([v3, edit], ann, gen,
                                         transcript_mode="all")
    check("a DNA variant and an RNA edit combine", bool(mixed))
    if mixed:
        m = mixed[0]
        check("the mix is flagged cross-evidence",
              m.extra["cross_evidence"] is True
              and m.extra["classes"] == ["RNA_EDITING", "SNV"],
              str(m.extra["classes"]))
        check("the evidence types are named in the notes",
              any(n.startswith("cross_evidence_") for n in m.notes),
              str(m.notes))

    # --------------------------------------------------------- refusals
    dup = CombinableVariant("chrT1", p3, r3, a3, "SNV", "vcf")
    over = build_combinatorial_proteins([v3, dup], ann, gen,
                                        transcript_mode="all")
    check("overlapping variants are refused, not mangled", over == [],
          f"{len(over)} record(s)")

    many = []
    for i in range(2, 12):
        gp = gpos_plus(meta, "GPLUS", i, 1)
        b = gen.fetch("chrT1", gp, gp)
        many.append(CombinableVariant("chrT1", gp, b, SWAP[b], "SNV", "vcf"))
    check("a transcript above the variant cap is skipped as an artefact",
          build_combinatorial_proteins(many, ann, gen, transcript_mode="all",
                                       max_variants=3) == [])

    check("an empty variant list yields nothing",
          build_combinatorial_proteins([], ann, gen) == [])

    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
