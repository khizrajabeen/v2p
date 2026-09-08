#!/usr/bin/env python3
"""Correctness tests against the synthetic reference.

Run:  python tests/test_pipeline.py
Exit code 0 = all pass. No pytest dependency so it runs anywhere.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

from make_fixture import build                             # noqa: E402
from v2p.annotation import Annotation, Genome              # noqa: E402
from v2p.build.fusion import build_fusion_proteins         # noqa: E402
from v2p.build.smallvar import build_small_variant_proteins  # noqa: E402
from v2p.build.splicing import (                           # noqa: E402
    event_junctions, transcript_junctions,
)
from v2p.fasta import HEADER_STYLES, write_fasta           # noqa: E402
from v2p.parse.inputs import (                             # noqa: E402
    parse_aachange, parse_suppa_event_id,
)
from v2p.seqops import (                                   # noqa: E402
    hgvs_p_to_change, revcomp, translate_orf,
)

FIX = ROOT / "tests" / "fixtures"
PASS: list[str] = []
FAIL: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    (PASS if cond else FAIL).append(f"{name}{' :: ' + detail if detail else ''}")
    print(f"[{'PASS' if cond else 'FAIL'}] {name}{' | ' + detail if detail else ''}")


def gpos_plus(meta, gene: str, codon: int, base: int = 1) -> int:
    """Genomic position of `base` of `codon` for a + strand fixture gene."""
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


def gpos_minus(meta, gene: str, codon: int, base: int = 1) -> int:
    m = meta[gene] if gene in meta else next(
        v for k, v in meta.items() if k.split("@")[0] == gene)
    tx_off = m["utr5"] + 3 * (codon - 1) + (base - 1)
    seen = 0
    for s, e in sorted(m["exons"], reverse=True):
        ln = e - s + 1
        if tx_off < seen + ln:
            return e - (tx_off - seen)
        seen += ln
    raise IndexError


def main() -> int:
    meta = build()
    ann = Annotation.from_gtf(FIX / "mini.gtf")
    gen = Genome(FIX / "mini.fa")

    # ---------------------------------------------------------------- unit
    check("revcomp", revcomp("ACGTN") == "NACGT")
    check("translate stop handling",
          translate_orf("ATGAAATAAGGG").protein == "MK"
          and translate_orf("ATGAAATAAGGG").stop_found)
    check("translate ambiguous -> X", translate_orf("ATGNNNAAA").protein == "MXK")
    check("hgvs 1-letter", hgvs_p_to_change("p.K509E") == ("K", 509, "E"))
    check("hgvs 3-letter", hgvs_p_to_change("p.Lys509Glu") == ("K", 509, "E"))
    aac = parse_aachange("POTEE:NM_001083538:exon14:c.A1525G:p.K509E")
    check("AAChange parse",
          aac and aac[0]["transcript"] == "NM_001083538"
          and aac[0]["protein"] == "p.K509E")
    ev = parse_suppa_event_id(
        "ENSG00000004455.18;RI:chr1:33007986:33010833-33013207:33013402:-")
    check("SUPPA RI id parse",
          ev and ev["event_type"] == "RI" and ev["strand"] == "-"
          and ev["coord_groups"] == [[33007986], [33010833, 33013207], [33013402]])
    ev2 = parse_suppa_event_id(
        "ENSG00000004534.15;AF:chr3:49940007:49940225-49962576:49940461:"
        "49940747-49962576:+")
    check("SUPPA AF id parse", ev2 and len(ev2["coord_groups"]) == 4)

    # -------------------------------------------------- plus-strand missense
    ref_plus = meta["GPLUS@chrT1"]["protein"]              # MASTKLYVWQDEGHRNFIPC
    p = gpos_plus(meta, "GPLUS", 5, 1)               # codon 5 = K (AAA)
    recs = build_small_variant_proteins("chrT1", p, "A", "C", ann, gen,
                                        variant_class="SNV", genes="GPLUS")
    ok = (len(recs) == 1 and recs[0].consequence == "missense"
          and recs[0].protein_change == "p.K5Q"
          and recs[0].sequence == ref_plus[:4] + "Q" + ref_plus[5:])
    check("plus-strand missense K5Q", ok,
          recs[0].protein_change if recs else "no record")

    # ------------------------------------------------- minus-strand missense
    ref_minus = meta["GMINUS@chrT2"]["protein"]            # MGKTVLYSEDAWQRNFHIP
    # codon 3 = K (AAA on transcript) -> genomic is TTT on the + strand
    pm = gpos_minus(meta, "GMINUS", 3, 1)
    obs = gen.fetch("chrT2", pm, pm)
    check("minus-strand ref base is complement", obs == "T", f"got {obs}")
    recs = build_small_variant_proteins("chrT2", pm, "T", "C", ann, gen,
                                        variant_class="RNA_EDITING",
                                        genes="GMINUS")
    ok = (len(recs) == 1 and recs[0].protein_change == "p.K3E"
          and recs[0].sequence == ref_minus[:2] + "E" + ref_minus[3:])
    check("minus-strand A-to-I edit K3E (genomic T>C)", ok,
          recs[0].protein_change if recs else "no record")

    # ----------------------------------------------------------- synonymous
    # codon 5 AAA -> AAG is still K
    p3 = gpos_plus(meta, "GPLUS", 5, 3)
    recs = build_small_variant_proteins("chrT1", p3, "A", "G", ann, gen,
                                        variant_class="SNV", genes="GPLUS")
    check("synonymous detected",
          recs and recs[0].consequence == "synonymous"
          and recs[0].sequence == ref_plus,
          recs[0].consequence if recs else "no record")

    # --------------------------------------------------------- stop gained
    # codon 7 = Y (TAT) -> TAA by changing base 3 T>A
    p7 = gpos_plus(meta, "GPLUS", 7, 3)
    recs = build_small_variant_proteins("chrT1", p7, "T", "A", ann, gen,
                                        variant_class="SNV", genes="GPLUS")
    check("stop gained truncates protein",
          recs and recs[0].consequence == "stop_gained"
          and recs[0].sequence == ref_plus[:6],
          f"{recs[0].consequence}/{recs[0].sequence}" if recs else "no record")

    # ------------------------- regression: stop_gained written as a deletion
    # A premature stop makes the alt protein end at `idx`, so the first
    # difference landed past its end and the protein change fell through to
    # the residue-deletion branch: a nonsense variant was reported as
    # `p.Y7del`, contradicting its own CSQ=stop_gained. HGVS spells it
    # `p.Tyr7Ter`, one-letter `p.Y7*`. Found by the benchmark, where it
    # showed as 0% recall across all 87 scored stop_gained rows.
    check("regression: stop_gained is written p.Y7*, not p.Y7del",
          recs and recs[0].protein_change == "p.Y7*",
          recs[0].protein_change if recs else "no record")

    # ---------------------------------------------------------- frameshift
    # delete 1 base at codon 6 -> everything downstream shifts
    p6 = gpos_plus(meta, "GPLUS", 6, 1)
    ref2 = gen.fetch("chrT1", p6 - 1, p6)            # VCF-style anchored del
    recs = build_small_variant_proteins("chrT1", p6 - 1, ref2, ref2[0], ann, gen,
                                        variant_class="INDEL", genes="GPLUS")
    fs = recs[0] if recs else None
    check("1bp deletion is frameshift",
          fs is not None and fs.consequence == "frameshift"
          and fs.sequence.startswith(ref_plus[:5])
          and fs.sequence != ref_plus,
          f"{fs.consequence} {fs.protein_change}" if fs else "no record")

    # ------------------------------------------------- in-frame 3bp deletion
    ref3 = gen.fetch("chrT1", p6 - 1, p6 + 2)        # anchor + 3 deleted bases
    recs = build_small_variant_proteins("chrT1", p6 - 1, ref3, ref3[0], ann, gen,
                                        variant_class="INDEL", genes="GPLUS")
    d = recs[0] if recs else None
    check("3bp deletion is in-frame and one residue shorter",
          d is not None and d.consequence == "inframe_deletion"
          and len(d.sequence) == len(ref_plus) - 1,
          f"{d.consequence} len={len(d.sequence)}" if d else "no record")

    # -------------------------------------------------- in-frame 3bp insertion
    anchor = gen.fetch("chrT1", p6 - 1, p6 - 1)
    recs = build_small_variant_proteins("chrT1", p6 - 1, anchor, anchor + "GGT",
                                        ann, gen, variant_class="INDEL",
                                        genes="GPLUS")
    ins = recs[0] if recs else None
    check("3bp insertion is in-frame and one residue longer",
          ins is not None and ins.consequence == "inframe_insertion"
          and len(ins.sequence) == len(ref_plus) + 1,
          f"{ins.consequence} len={len(ins.sequence)}" if ins else "no record")

    # ------------------------------------------------------ intronic is null
    recs = build_small_variant_proteins("chrT1", 250, "N", "N", ann, gen,
                                        variant_class="SNV", genes="GPLUS")
    check("intronic variant produces no protein", recs == [])

    # ----------------------------------------------------------- REF mismatch
    wrong = "A" if gen.fetch("chrT1", p, p) != "A" else "C"
    recs = build_small_variant_proteins("chrT1", p, wrong, "G", ann, gen,
                                        variant_class="SNV", genes="GPLUS")
    check("REF mismatch is flagged, not silent",
          recs and any(n.startswith("REF_MISMATCH") for n in recs[0].notes))

    # ------------------------------------------- positional lookup (no gene)
    # A sites-only VCF supplies no gene symbol; the same variant must give
    # the same protein when found by coordinate alone.
    recs_named = build_small_variant_proteins("chrT1", p, "A", "C", ann, gen,
                                              variant_class="SNV", genes="GPLUS")
    recs_pos = build_small_variant_proteins("chrT1", p, "A", "C", ann, gen,
                                            variant_class="SNV", genes="")
    check("positional lookup finds the transcript with no gene symbol",
          len(recs_pos) == 1 and recs_pos[0].sequence == recs_named[0].sequence
          and recs_pos[0].protein_change == "p.K5Q",
          f"{len(recs_pos)} record(s)")

    recs_pos_m = build_small_variant_proteins("chrT2", pm, "T", "C", ann, gen,
                                              variant_class="SNV", genes="")
    check("positional lookup works on the minus strand",
          len(recs_pos_m) == 1 and recs_pos_m[0].protein_change == "p.K3E")

    check("positional lookup in an intergenic gap returns nothing",
          build_small_variant_proteins("chrT1", 2900, "A", "C", ann, gen,
                                       genes="") == [])

    check("stale gene symbol falls back to coordinates",
          len(build_small_variant_proteins("chrT1", p, "A", "C", ann, gen,
                                           genes="NOT_A_REAL_GENE")) == 1)

    check("coding_only GTF parse keeps the coding transcripts",
          len(Annotation.from_gtf(FIX / "mini.gtf", coding_only=True).tx)
          == len(meta))

    # ---------------------------------------------------------------- fusion
    # GPLUS exonic bp at end of codon 10, GPART3 bp at start of its exon 2
    bp5 = gpos_plus(meta, "GPLUS", 10, 3)
    bp3 = meta["GPART3@chrT3"]["exons"][1][0]
    recs = build_fusion_proteins("GPLUS", "chrT1", bp5, "GPART3", "chrT3", bp3,
                                 ann, gen, try_both_orientations=False)
    f = recs[0] if recs else None
    check("fusion keeps 5' partner residues and junction index",
          f is not None and f.sequence.startswith(ref_plus[:10])
          and f.variant_pos_aa == 11
          and "junction_on_codon_boundary" in f.notes,
          f"{f.protein_change} {f.notes}" if f else "no record")
    check("fusion reports 3' partner frame status",
          f is not None and any(n in ("3p_native_frame", "3p_frameshifted")
                                for n in f.notes))

    # intronic 5' breakpoint should snap to the exon end
    recs = build_fusion_proteins("GPLUS", "chrT1", 200, "GPART3", "chrT3", bp3,
                                 ann, gen, try_both_orientations=False)
    check("intronic fusion breakpoint snaps to splice donor",
          recs and "bp5_intronic_snapped_to_exon_end" in recs[0].notes)

    # ---------------------------------------------------------- AS grammar
    t = ann.representative("GPLUS")
    jx = transcript_junctions(t)
    check("junction extraction", jx == {(152, 353)}, str(jx))
    se = event_junctions("SE", [[100, 200], [300, 400]])
    check("SE grammar", se == {"form1": [(100, 200), (300, 400)],
                               "form2": [(100, 400)]})
    check("bad arity rejected", event_junctions("SE", [[1, 2]]) is None)
    mx = event_junctions("MX", [[1, 2], [3, 4], [5, 6], [7, 8]])
    check("MX grammar", mx is not None and len(mx["form1"]) == 2)

    # SUPPA2 orders AF/AL groups differently per strand; both must parse
    afm = parse_suppa_event_id(
        "E;AF:chr7:26854890-26864363:26864590:26854890-26864842:26865113:-")
    afp = parse_suppa_event_id(
        "E;AF:chr19:29942261:29942664-29985223:29983233:29983316-29985223:+")
    alm = parse_suppa_event_id(
        "E;AL:chr3:149964904:149966264-149968358:149966268:149966586-149968358:-")
    jm = event_junctions("AF", afm["coord_groups"])
    jp = event_junctions("AF", afp["coord_groups"])
    jl = event_junctions("AL", alm["coord_groups"])
    check("AF minus-strand group order parses",
          jm == {"form1": [(26854890, 26864363)],
                 "form2": [(26854890, 26864842)]}, str(jm))
    check("AF plus-strand group order parses",
          jp == {"form1": [(29942664, 29985223)],
                 "form2": [(29983316, 29985223)]}, str(jp))
    check("AL minus-strand group order parses",
          jl == {"form1": [(149966264, 149968358)],
                 "form2": [(149966586, 149968358)]}, str(jl))
    check("AF with wrong pair count is rejected",
          event_junctions("AF", [[1], [2], [3], [4]]) is None)

    # ------------------------------------------------------- FASTA emission
    recs = build_small_variant_proteins("chrT1", p, "A", "C", ann, gen,
                                        variant_class="SNV", genes="GPLUS")
    for style in HEADER_STYLES:
        out = FIX / f"out.{style}.fasta"
        counts = write_fasta(list(recs), out, style=style)
        txt = out.read_text()
        body = [ln for ln in txt.splitlines() if not ln.startswith("#")]
        ok = (counts["_total"] == 1 and body[0].startswith(">")
              and body[1] == recs[0].sequence)
        vt = ("VT=SNV" in txt) or ("VariantType=SNV" in txt) or (".SNV." in txt)
        check(f"FASTA style '{style}' round-trips with variant type",
              ok and vt)

    # --------------------------------------- GTF frame / cds_start_NF
    # GENCODE sets a non-zero frame on the first CDS block of transcripts
    # whose 5' CDS end is undetermined. Ignoring it shifts the reading
    # frame and yields a wrong but plausible protein.
    t_plus = ann.representative("GPLUS")
    base_off = t_plus.cds_offset_in_tx()
    t_plus.cds_phase_first = 2
    shifted_off = t_plus.cds_offset_in_tx()
    t_plus.cds_phase_first = 0
    check("first-CDS frame is applied to the translation start",
          shifted_off == base_off + 2, f"{base_off} -> {shifted_off}")

    # --------------------------------- ambiguous gene symbol resolution
    amb = ann.ambiguous_gene_names()
    check("duplicate gene symbols are detected", "GDUP" in amb, str(amb))
    t5 = ann.representative_at("GDUP", "chrT5", 150)
    t6 = ann.representative_at("GDUP", "chrT6", 150)
    check("an ambiguous symbol resolves to the locus, not the first match",
          t5 is not None and t6 is not None and t5.tx_id != t6.tx_id
          and t5.chrom == "chrT5" and t6.chrom == "chrT6",
          f"{t5.tx_id if t5 else None} vs {t6.tx_id if t6 else None}")
    check("a locus with no copy of the gene returns nothing",
          ann.representative_at("GDUP", "chrT1", 150) is None)

    # a variant given the ambiguous symbol must still hit the right copy
    p5 = gpos_plus({"x": meta["GDUP@chrT5"]}, "x", 5, 1)
    r5 = build_small_variant_proteins("chrT5", p5, "A", "C", ann, gen,
                                      variant_class="SNV", genes="GDUP")
    check("ambiguous symbol + locus gives the correct protein",
          len(r5) == 1 and r5[0].protein_change == "p.K5Q"
          and r5[0].transcript.endswith("chrT5"),
          f"{r5[0].protein_change} on {r5[0].transcript}" if r5 else "none")

    # ------------------------------------- selenocysteine and chrM code
    from v2p.seqops import CODON_TABLE_MITO, table_for_contig
    t_sec = ann.representative("GSEC")
    check("GENCODE Selenocysteine features are parsed",
          t_sec is not None and len(t_sec.sec_sites) == 1,
          str(t_sec.sec_sites if t_sec else None))
    sec_prot, _ = None, None
    from v2p.build.splicing import _translate_transcript
    sec_prot, _n = _translate_transcript(t_sec, gen)
    check("a UGA annotated as Sec becomes U, not a truncation",
          sec_prot == meta["GSEC@chrT4"]["protein"],
          f"{sec_prot} vs {meta['GSEC@chrT4']['protein']}")
    check("without the Sec annotation the protein would truncate early",
          len(translate_orf(gen.blocks(t_sec.chrom, t_sec.exons, t_sec.strand)
                            [t_sec.cds_offset_in_tx():]).protein)
          < len(meta["GSEC@chrT4"]["protein"]))

    t_mito = ann.representative("GMITO")
    check("chrM selects the vertebrate mitochondrial code",
          table_for_contig(t_mito.chrom) is CODON_TABLE_MITO)
    mito_prot, _ = _translate_transcript(t_mito, gen)
    check("chrM translation uses table 2 (TGA=W)",
          mito_prot == meta["GMITO@chrM"]["protein"],
          f"{mito_prot} vs {meta['GMITO@chrM']['protein']}")

    # ------------------------------------------------ UTR vs synonymous
    # A UTR variant is exonic and leaves the protein unchanged, exactly
    # like a synonymous one. Conflating them inverts the nonsynonymous-
    # to-synonymous ratio, which is a standard sanity check on a somatic
    # call set.
    m_plus = meta["GPLUS@chrT1"]
    utr5_g = m_plus["exons"][0][0] + 5          # inside the 21 nt 5'UTR
    r5 = build_small_variant_proteins("chrT1", utr5_g,
                                      gen.fetch("chrT1", utr5_g, utr5_g),
                                      "G" if gen.fetch("chrT1", utr5_g, utr5_g) != "G" else "C",
                                      ann, gen, genes="GPLUS")
    check("a 5'UTR variant is labelled 5_prime_UTR, not synonymous",
          len(r5) == 1 and r5[0].consequence == "5_prime_UTR",
          r5[0].consequence if r5 else "no record")

    utr3_g = m_plus["exons"][1][1] - 5          # inside the 30 nt 3'UTR
    r3 = build_small_variant_proteins("chrT1", utr3_g,
                                      gen.fetch("chrT1", utr3_g, utr3_g),
                                      "G" if gen.fetch("chrT1", utr3_g, utr3_g) != "G" else "C",
                                      ann, gen, genes="GPLUS")
    check("a 3'UTR variant is labelled 3_prime_UTR, not synonymous",
          len(r3) == 1 and r3[0].consequence == "3_prime_UTR",
          r3[0].consequence if r3 else "no record")

    check("a true synonymous change is still called synonymous",
          build_small_variant_proteins("chrT1", p3, "A", "G", ann, gen,
                                       genes="GPLUS")[0].consequence
          == "synonymous")
    check("UTR and synonymous variants all leave the protein unchanged",
          r5[0].sequence == ref_plus and r3[0].sequence == ref_plus)

    # -------------------------------------- disposition classification
    # A variant that yields nothing leaves no trace in the FASTA, so the
    # reason has to be recorded at the point of decision.
    dsyn = build_small_variant_proteins("chrT1", p3, "A", "G", ann, gen,
                                        variant_class="SNV", genes="GPLUS")
    check("a synonymous variant still yields a record to classify",
          len(dsyn) == 1 and dsyn[0].consequence == "synonymous")
    check("an intronic position yields no record at all",
          build_small_variant_proteins("chrT1", 250, "A", "C", ann, gen,
                                       genes="GPLUS") == [])
    check("transcripts_at distinguishes intronic from intergenic",
          len(ann.transcripts_at("chrT1", 250)) >= 1
          and ann.transcripts_at("chrT1", 2900) == [])

    # ------------------------------------------------- peptides and NMD
    from v2p.peptides import (digest, decoy_pseudo_reverse, decoy_reverse,
                              novel_peptides, reference_peptide_space,
                              variant_spanning_peptides)
    from v2p.nmd import junction_offsets, predict_nmd

    d = digest("MAAAKBBBBRCCCCKPDDDDR", missed_cleavages=0, min_len=1)
    check("trypsin cleaves after K/R but not before P",
          "MAAAK" in d and "CCCCKPDDDDR" in d, str(sorted(d)))
    check("missed cleavages produce longer peptides",
          len(digest("MAAAKBBBBRCCCCR", missed_cleavages=2, min_len=1))
          > len(digest("MAAAKBBBBRCCCCR", missed_cleavages=0, min_len=1)))

    ref = ["MAAAKTTTTRQQQQK", "MSSSSKYYYYR"]
    space = reference_peptide_space(ref, min_len=4)
    check("a reference protein contributes no novel peptide",
          novel_peptides("MAAAKTTTTRQQQQK", space, min_len=4) == set())
    nv = novel_peptides("MAAAKTWTTRQQQQK", space, min_len=4)
    check("a substitution creates a novel peptide", len(nv) > 0, str(nv))
    check("I/L are equated (isobaric, indistinguishable by MS)",
          novel_peptides("MAAAKTTTTRQQQQK".replace("T", "I"),
                         reference_peptide_space(
                             ["MAAAKTTTTRQQQQK".replace("T", "L")], min_len=4),
                         min_len=4) == set())

    prot = "MAAAKTWTTRQQQQK"
    check("variant-spanning filter keeps only covering peptides",
          all(p_.find("W") >= 0 or prot.find(p_) <= 6
              for p_ in variant_spanning_peptides(prot, 7, nv)))

    check("pseudo-reverse decoy keeps tryptic C-termini",
          decoy_pseudo_reverse("MAAAKBBBBR").endswith("R")
          and decoy_pseudo_reverse("MAAAKBBBBR")[4] == "K",
          decoy_pseudo_reverse("MAAAKBBBBR"))
    check("pseudo-reverse preserves peptide count",
          len(digest(decoy_pseudo_reverse("MAAAKBBBBRCCCCR"), 0, 1))
          == len(digest("MAAAKBBBBRCCCCR", 0, 1)))
    check("plain reverse is a different decoy",
          decoy_reverse("MAAAK") == "KAAAM")

    check("junction offsets", junction_offsets([100, 200, 150]) == [100, 300])
    check("premature stop far upstream is NMD-likely",
          predict_nmd(100, [150, 200, 150])["nmd"] == "likely")
    check("stop in the last exon escapes NMD",
          predict_nmd(400, [150, 200, 150])["nmd"] == "escapes")
    check("stop within 50nt of the last junction escapes",
          predict_nmd(320, [150, 200, 150])["nmd"] == "escapes")
    check("single-exon transcripts always escape",
          predict_nmd(100, [500])["nmd"] == "escapes")
    check("an upstream indel shifts the junctions",
          predict_nmd(300, [150, 200, 150], variant_tx_offset=50,
                      length_delta=-60)["nmd"] == "escapes")

    # reference entries must look like SwissProt, not like variants
    from v2p.build.smallvar import ProteinRecord as PR
    ref_rec = PR(seq_id="REF_X", sequence="MASTKLYVWQ", variant_class="REFERENCE",
                 consequence="reference", gene="GPLUS", transcript="T1",
                 locus="chrT1:134A>C",
                 extra={"accession": "P00001", "description": "Test protein"})
    out = FIX / "out.refhdr.fasta"
    write_fasta([ref_rec], out, style="uniprot")
    h = out.read_text().splitlines()[0]
    check("reference entry uses sp| and omits the variant locus",
          h.startswith(">sp|P00001|GPLUS_HUMAN")
          and "chrT1:134" not in h and "variant protein" not in h, h[:90])

    # ------------------------------------------- input type detection
    # Files arrive named differently from every collaborator, so type
    # detection must run on content. A wrong guess here silently builds a
    # database from the wrong data.
    import gzip as _gz, tempfile
    from v2p.discover import (build_run_plan, detect_file, detect_folder,
                              FUSION as _F, GENOME as _G, RNA_EDITING as _R,
                              SMALL_VARIANTS as _S, SPLICING as _SP,
                              ANNOTATION as _A, PROTEOME as _P, INDEX as _I,
                              REGIONS as _RG, UNKNOWN as _U)

    tmp = Path(tempfile.mkdtemp())
    (tmp / "a.txt").write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t.\tPASS\t.\n")
    (tmp / "b.csv").write_text(
        "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\n"
        "chr1\t10\t.\tA\tG\tintronic\nchr2\t20\t.\tT\tC\tintronic\n")
    (tmp / "c").write_text(
        "AS,dPSI\nENSG1;SE:chr1:10-20:30-40:+,0.3\n"
        "ENSG2;AF:chr2:1:2-3:4:5-6:-,0.2\nENSG3;RI:chr3:1:2-3:4:+,0.1\n")
    (tmp / "d.csv").write_text(
        "tag,Gene1,Gene2,breakpoint1,breakpoint2\n"
        "t,AAA,BBB,chr1:100,chr5:200\nu,CCC,DDD,chr2:300,chr6:400\n")
    (tmp / "e.bed").write_text("chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t50\n")
    (tmp / "f.fa").write_text(">chr1\nACGTACGTACGTACGTACGT\n")
    (tmp / "g.fasta").write_text(
        ">sp|P1|X_HUMAN Protein OS=Homo sapiens OX=9606 GN=X\nMKWVTFISLL\n")
    with _gz.open(tmp / "h.gz", "wt") as fh:
        fh.write('chr1\ts\texon\t1\t9\t.\t+\t.\tgene_id "G"; transcript_id "T";\n')
    (tmp / "i.csi").write_bytes(b"CSI\x01" + b"\x00" * 16)
    (tmp / "j.csv").write_text("who,when\nkhizra,today\n")

    want = {"a.txt": _S, "b.csv": _R, "c": _SP, "d.csv": _F, "e.bed": _RG,
            "f.fa": _G, "g.fasta": _P, "h.gz": _A, "i.csi": _I, "j.csv": _U}
    got = {d.path.name: d.role for d in detect_folder(tmp)}
    for fn, role in want.items():
        check(f"detect {fn} -> {role}", got.get(fn) == role,
              f"got {got.get(fn)}")

    plan = build_run_plan(detect_folder(tmp))
    check("all four evidence types are found",
          set(plan["evidence_types"]) == {_S, _R, _F, _SP},
          str(plan["evidence_types"]))
    check("genome build is read from the VCF contig lengths",
          plan.get("genome_build") == "GRCh38", str(plan.get("genome_build")))
    check("a clean folder reports no problems", plan["problems"] == [],
          str(plan["problems"]))

    (tmp / "a2.txt").write_text((tmp / "a.txt").read_text())
    plan2 = build_run_plan(detect_folder(tmp))
    check("two equally plausible VCFs are reported, not silently picked",
          any("candidates of similar confidence" in p_
              for p_ in plan2["problems"]), str(plan2["problems"]))

    # companion index files share the shape of the data they index
    (tmp / "f.fa.fai").write_text("chr1\t20\t6\t60\t61\n")
    got2 = {d.path.name: d.role for d in detect_folder(tmp)}
    check("a .fai index is not mistaken for a BED", got2.get("f.fa.fai") == _I,
          str(got2.get("f.fa.fai")))
    check("a real BED is still a BED", got2.get("e.bed") == _RG)

    # the same file compressed is not a competing candidate
    with _gz.open(tmp / "f.fa.gz", "wt") as fh:
        fh.write((tmp / "f.fa").read_text())
    plan3 = build_run_plan(detect_folder(tmp))
    check("X and X.gz are one candidate, not an ambiguity",
          not any("reference genome" in p_ for p_ in plan3["problems"]),
          str(plan3["problems"]))
    check("the uncompressed copy is preferred for random access",
          plan3["selected"][_G].endswith("f.fa"),
          plan3["selected"].get(_G, ""))

    empty = Path(tempfile.mkdtemp())
    (empty / "notes.txt").write_text("nothing here\n")
    check("a folder with no evidence is refused",
          any("no variant evidence" in p_
              for p_ in build_run_plan(detect_folder(empty))["problems"]))

    (tmp / "a.txt").write_text(
        (tmp / "a.txt").read_text().replace("248956422", "249250621"))
    (tmp / "a2.txt").unlink()
    check("conflicting genome builds are detected",
          build_run_plan(detect_folder(tmp)).get("genome_build") == "GRCh37")

    # ------------------------------------- regression: run ignored --ref
    # `cmd_run` built the run plan from the input folder alone, so the
    # genome and annotation that --ref supplies were reported missing and
    # the documented invocation `v2p run <folder> --ref <dir>` refused to
    # run at all without --force. `cmd_detect` already scanned --ref;
    # `cmd_run` did not.
    import argparse as _ap

    # The CLI lives in the package so `pip install v2p` can expose it as a
    # console script; scripts/v2p is a shim over this same module.
    from v2p import cli as _cli

    rt = Path(tempfile.mkdtemp())
    (rt / "inputs").mkdir()
    (rt / "refs").mkdir()
    (rt / "inputs" / "calls.vcf").write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t.\tPASS\t.\n")
    (rt / "refs" / "genome.fa").write_text(">chr1\nACGTACGTACGTACGTACGT\n")
    (rt / "refs" / "proteome.fasta").write_text(
        ">sp|P1|X_HUMAN Protein OS=Homo sapiens OX=9606 GN=X\nMKWVTFISLL\n")
    with _gz.open(rt / "refs" / "genes.gtf.gz", "wt") as fh:
        fh.write('chr1\ts\texon\t1\t9\t.\t+\t.\t'
                 'gene_id "G"; transcript_id "T";\n')

    class _Started(Exception):
        """Stands in for the first pipeline stage, so the test stops exactly
        where the refusal used to happen."""

    def _stub(script, args, quiet=False, logdir=None):
        raise _Started(script)

    _cli._run = _stub

    def _ns(ref):
        return _ap.Namespace(
            folder=str(rt / "inputs"), ref=ref, outdir=str(rt / "out"),
            name="T", transcript_mode="all", header_style="uniprot",
            keep_unchanged=True, decoys="none", split_by_type=True,
            append_reference=True, no_recursive=False, force=False,
            min_agreement=0.90, logdir=str(rt / "logs"),
            skip_invariant=None)

    try:
        _cli.cmd_run(_ns(str(rt / "refs")))
        outcome = "refused without running anything"
    except _Started as e:
        outcome = f"started {e}"
    check("regression: run --ref no longer refuses the references it was "
          "given", outcome.startswith("started"), outcome)
    check("regression: run without --ref still refuses rather than guessing",
          _cli.cmd_run(_ns(None)) == 2)

    # ------------- regression: audit could not read its own package source
    # The audit checks what the pipeline does by reading our own modules.
    # The path was `parents[1] / "src" / "v2p" / ...`, correct while the
    # stage scripts sat in scripts/ and silently wrong once they moved into
    # src/v2p/stages/ - it resolved to src/v2p/src/v2p/..., the read raised
    # OSError, and a bare except turned "cannot check" into "the pipeline
    # does not do it". The audit then reported three long-fixed bugs as
    # live, which would have failed the reference-audit workflow forever.
    import importlib.util as _iu
    _spec = _iu.spec_from_file_location(
        "v2p_audit",
        ROOT / "src" / "v2p" / "stages" / "00_audit_references.py")
    _aud = _iu.module_from_spec(_spec)
    _spec.loader.exec_module(_aud)

    for _rel, _needle in (("annotation.py", "cds_phase_first"),
                          ("seqops.py", "CODON_TABLE_MITO"),
                          ("build/smallvar.py", "representative_at")):
        try:
            _got = _needle in _aud._pkg_source(_rel)
        except OSError as _e:
            _got = f"OSError: {_e}"
        check(f"regression: the audit can read its own {_rel}",
              _got is True, str(_got))

    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("FAILURES:")
        for f_ in FAIL:
            print("  -", f_)
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
