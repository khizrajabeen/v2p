#!/usr/bin/env python3
"""Milestone 1 acceptance test - the nine release invariants.

Run:  python tests/test_invariants.py
Exit code 0 = all pass. No pytest dependency, no network, no reference
files: the miniature release below is built from scratch in a temp dir.

The shape of the test is one positive and one negative case per check.
The positive case is the same clean release every time, asserted to produce
no violation of that check. The negative case injects exactly one fault and
asserts that check fires - and, just as importantly, that no *other* check
fires, because a check that reports everything reports nothing.
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2p.invariants import (                               # noqa: E402
    Layout, _intralocus_note, check_release, fasta_layout, format_report,
    header_fields, layout_matches, locate_release,
)
from v2p.provenance import sha256                          # noqa: E402

PASS: list[str] = []
FAIL: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    (PASS if cond else FAIL).append(name)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}"
          f"{' | ' + detail if detail else ''}")


# --------------------------------------------------------------------------
# a miniature release
# --------------------------------------------------------------------------

NAME = "MINI"

# Three variant entries, two wild-type counterparts, two appended reference
# proteins. Sequences are short and unwrapped, like the real output.
VARIANTS = [
    (f">vr|{NAME}_SNV_000001|AAA_HUMAN_SNV Alpha (p.K5E) OS=Homo sapiens "
     f"OX=9606 GN=AAA VT=SNV CSQ=missense TX=ENST00000000001.1 PC=p.K5E "
     f"POS=5", "MADQEFGHIKLMNPQRSTVWY"),
    (f">vr|{NAME}_INDEL_000002|BBB_HUMAN_INDEL Beta (p.G7fs) OS=Homo sapiens "
     f"OX=9606 GN=BBB VT=INDEL CSQ=frameshift TX=ENST00000000002.1 "
     f"PC=p.G7fs POS=7", "MCDEFGWWWWWWTT"),
    (f">vr|{NAME}_RNA_EDITING_000003|CCC_HUMAN_RNA_EDITING Gamma (p.Q10R) "
     f"OS=Homo sapiens OX=9606 GN=CCC VT=RNA_EDITING CSQ=missense "
     f"TX=ENST00000000003.1 PC=p.Q10R POS=10", "MSTUVWYACDRFGHIKL"),
]
WILDTYPES = [
    (f">sp|REF_AAA_ENST00000000001.1|AAA_HUMAN AAA reference protein "
     f"OS=Homo sapiens OX=9606 GN=AAA VT=REFERENCE CSQ=reference "
     f"TX=ENST00000000001.1", "MADQKFGHIKLMNPQRSTVWY"),
    (f">sp|REF_BBB_ENST00000000002.1|BBB_HUMAN BBB reference protein "
     f"OS=Homo sapiens OX=9606 GN=BBB VT=REFERENCE CSQ=reference "
     f"TX=ENST00000000002.1", "MCDEFGGHIKLMNPQ"),
]
APPENDED = [
    (">sp|P00001|DDD_HUMAN Delta OS=Homo sapiens OX=9606 GN=DDD",
     "MDDDEFGHIKLMNPQRST"),
    (">sp|P00002|EEE_HUMAN Epsilon OS=Homo sapiens OX=9606 GN=EEE",
     "MEEEFGHIKLMNPQRSTVW"),
]

# One manifest row per input variant, and one disposition row per manifest
# row. Some inputs produce no protein - that is a normal outcome, and they
# must still be accounted for.
MANIFEST_IDS = ["SNV|chr1:100A>G", "INDEL|chr1:200AT>A",
                "RNA_EDITING|chr2:300A>G", "SNV|chr3:400C>T",
                "AS_SE|chr4:500-600"]

VALIDATION_MD = """# UniProt validation report

Reference: `mini.fasta` - 2 entries, 2 genes.

## A. Annotation check - reference residue vs UniProt

2 transcript-level annotations across 2 unique sites.

## B1. Translation check vs GENCODE (authoritative)

2 transcripts compared against GENCODE's own translation of the same
transcript ids. Exact-agreement: **100.0%**.

| outcome | transcripts |
|---|---:|
| identical | 2 |

## B2. Comparison vs UniProt canonical (informational)

2 reference proteins. Same-protein rate: **100.0%**.
"""


def _write_fasta(path: Path, records, width: int = 0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for hdr, seq in records:
            fh.write(hdr + "\n")
            if width and len(seq) > width:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")


def _pseudo_reverse(seq: str) -> str:
    """Keep the last residue fixed and reverse the rest. Enough of a decoy
    for a count check; it need not match `peptides.py` exactly."""
    return seq[:-1][::-1] + seq[-1] if len(seq) > 1 else seq


def write_manifest(out: Path) -> None:
    """MANIFEST.txt in the exact grammar `06_package_release.py` writes."""
    man = [f"# MANIFEST - {NAME}", ""]
    for f in sorted(out.rglob("*")):
        if f.is_file() and f.name != "MANIFEST.txt":
            man.append(f"{sha256(f)}  {f.relative_to(out).as_posix()}  "
                       f"({f.stat().st_size:,} bytes)")
    (out / "MANIFEST.txt").write_text("\n".join(man) + "\n", encoding="utf-8")


def build_release(root: Path) -> tuple[Path, Path, Path, Path]:
    """A clean miniature release.

    Returns (outdir, manifest, disposition, proteome).
    """
    out = root / "release"
    work = out / "_work" / "tables"
    work.mkdir(parents=True, exist_ok=True)

    combined = VARIANTS + WILDTYPES + APPENDED
    _write_fasta(out / f"{NAME}.target.fasta", combined)
    _write_fasta(out / f"{NAME}.target_decoy.fasta",
                 combined + [(">DECOY_" + h[1:], _pseudo_reverse(s))
                             for h, s in combined])

    # The per-type split: one file per VT, plus the appended proteome.
    per: dict[str, list] = {}
    for hdr, seq in combined:
        vt = header_fields(hdr).get("VT")
        per.setdefault(vt if vt else "REFERENCE_PROTEOME", []).append(
            (hdr, seq))
    for cls, items in per.items():
        _write_fasta(out / "by_class" / f"{NAME}.{cls}.fasta", items)
    (out / "by_class" / "INDEX.md").write_text("# index\n", encoding="utf-8")

    # The entry table covers variant and wild-type entries only.
    rows = ["seq_id\tvariant_class\tgene"]
    for hdr, _s in VARIANTS + WILDTYPES:
        f = header_fields(hdr)
        rows.append(f"{hdr.split()[0][1:]}\t{f.get('VT', '')}\t"
                    f"{f.get('GN', '')}")
    (out / f"{NAME}.entries.tsv").write_text("\n".join(rows) + "\n",
                                             encoding="utf-8")

    manifest = work / "unified_variant_manifest.tsv"
    manifest.write_text(
        "variant_id\tvariant_class\tlocus\n"
        + "".join(f"{v}\tX\tchr1:1\n" for v in MANIFEST_IDS),
        encoding="utf-8")

    disposition = work / "disposition.tsv"
    disposition.write_text(
        "variant_id\tvariant_class\tn_proteins\toutcome\n"
        + "".join(f"{v}\tX\t1\tprotein_built\n" for v in MANIFEST_IDS),
        encoding="utf-8")

    (out / "qc").mkdir(exist_ok=True)
    (out / "qc" / "uniprot_validation.md").write_text(VALIDATION_MD,
                                                      encoding="utf-8")

    proteome = root / "mini_proteome.fasta"
    _write_fasta(proteome, APPENDED)

    write_manifest(out)
    return out, manifest, disposition, proteome


def paths_of(work: Path) -> tuple[Path, Path, Path, Path]:
    """The four `check_release` arguments for a copied release tree."""
    out = work / "release"
    tables = out / "_work" / "tables"
    return (out, tables / "unified_variant_manifest.tsv",
            tables / "disposition.tsv", work / "mini_proteome.fasta")


def run_checks(out: Path, manifest: Path, disposition: Path,
               proteome: Path) -> list:
    return check_release(out, manifest, disposition, proteome=proteome,
                         decoys_requested=True)


def fired(violations, severity: str | None = None) -> set[str]:
    """The check ids that reported something real.

    'not evaluated' notices are excluded: they say a check did not run,
    which is not the same as a check that found a fault.
    """
    return {v.check for v in violations
            if not v.message.startswith("not evaluated")
            and (severity is None or v.severity == severity)}


# --------------------------------------------------------------------------
# fault injectors - each breaks exactly one invariant
# --------------------------------------------------------------------------

def fault_i1(out: Path) -> None:
    """Put one sequence into two per-type files, so the split is no longer a
    partition. This is the shape of the double-counted wild-type bug."""
    src = out / "by_class" / f"{NAME}.SNV.fasta"
    dst = out / "by_class" / f"{NAME}.INDEL.fasta"
    with open(dst, "a", encoding="utf-8") as fh:
        fh.write(src.read_text(encoding="utf-8"))
    write_manifest(out)


def fault_i2(out: Path) -> None:
    """Drop one input variant's disposition row and duplicate another."""
    p = out / "_work" / "tables" / "disposition.tsv"
    lines = p.read_text(encoding="utf-8").splitlines()
    p.write_text("\n".join(lines[:-1] + [lines[1]]) + "\n", encoding="utf-8")
    write_manifest(out)


def fault_i3(out: Path) -> None:
    """Add a wild-type row to the entry table that no sequence backs, so the
    tabulated count no longer sums to the combined file."""
    with open(out / f"{NAME}.entries.tsv", "a", encoding="utf-8") as fh:
        fh.write("sp|REF_ZZZ_ENST00000000009.1|ZZZ_HUMAN\tREFERENCE\tZZZ\n")
    write_manifest(out)


def fault_i4(out: Path) -> None:
    """Give two entries the same accession. A search engine cannot tell them
    apart, so one of the two hits is unattributable.

    The collision is made by renaming rather than by adding an entry, so
    every count in the release is unchanged and this trips I4 alone.
    """
    # The clashing entry is also relabelled VT=SNV, so both copies of the id
    # sit in one per-type file. A duplicate id spread across two files would
    # be a genuine I1 overlap as well, and the point here is to isolate I4.
    dup_id = VARIANTS[0][0].split()[0]
    clash = (dup_id + " " + VARIANTS[2][0].split(" ", 1)[1].replace(
        "VT=RNA_EDITING", "VT=SNV"), VARIANTS[2][1])
    combined = [VARIANTS[0], VARIANTS[1], clash] + WILDTYPES + APPENDED
    _write_fasta(out / f"{NAME}.target.fasta", combined)
    _write_fasta(out / f"{NAME}.target_decoy.fasta",
                 combined + [(">DECOY_" + h[1:], _pseudo_reverse(s))
                             for h, s in combined])
    _write_fasta(out / "by_class" / f"{NAME}.SNV.fasta",
                 [VARIANTS[0], clash])
    (out / "by_class" / f"{NAME}.RNA_EDITING.fasta").unlink()
    write_manifest(out)


def fault_i5(out: Path) -> None:
    """Wrap the combined file at 10 while the supplied proteome template is
    unwrapped."""
    _write_fasta(out / f"{NAME}.target.fasta",
                 VARIANTS + WILDTYPES + APPENDED, width=10)
    write_manifest(out)


def fault_i6(out: Path) -> None:
    """Reference translations disagree with the annotation source."""
    (out / "qc" / "uniprot_validation.md").write_text(
        VALIDATION_MD.replace("Exact-agreement: **100.0%**.",
                              "Exact-agreement: **61.5%**."),
        encoding="utf-8")
    write_manifest(out)


def fault_i7(out: Path) -> None:
    """Label an entry `frameshift` whose protein is byte-identical to the
    reference of its own transcript. This is the UTR-as-frameshift bug."""
    hdr, _seq = VARIANTS[1]
    ref_seq = WILDTYPES[1][1]
    combined = ([VARIANTS[0], (hdr, ref_seq), VARIANTS[2]]
                + WILDTYPES + APPENDED)
    _write_fasta(out / f"{NAME}.target.fasta", combined)
    _write_fasta(out / "by_class" / f"{NAME}.INDEL.fasta", [(hdr, ref_seq)])
    _write_fasta(out / f"{NAME}.target_decoy.fasta",
                 combined + [(">DECOY_" + h[1:], _pseudo_reverse(s))
                             for h, s in combined])
    write_manifest(out)


def fault_i8(out: Path) -> None:
    """Change a file after the manifest was written - a stale manifest.

    `by_class/INDEX.md` is edited because no other check reads it, so the
    only thing that can notice is the checksum.
    """
    with open(out / "by_class" / "INDEX.md", "a", encoding="utf-8") as fh:
        fh.write("appended after the manifest was written\n")


def fault_i9(out: Path) -> None:
    """Drop one decoy, so the target and decoy spaces are different sizes and
    the FDR estimate is silently wrong."""
    p = out / f"{NAME}.target_decoy.fasta"
    text = p.read_text(encoding="utf-8").splitlines()
    cut = next(i for i, ln in enumerate(text) if ln.startswith(">DECOY_"))
    del text[cut:cut + 2]
    p.write_text("\n".join(text) + "\n", encoding="utf-8")
    write_manifest(out)


FAULTS = [("I1", fault_i1), ("I2", fault_i2), ("I3", fault_i3),
          ("I4", fault_i4), ("I5", fault_i5), ("I6", fault_i6),
          ("I7", fault_i7), ("I8", fault_i8), ("I9", fault_i9)]


# --------------------------------------------------------------------------

def main() -> int:
    root = Path(tempfile.mkdtemp(prefix="v2p_inv_"))
    try:
        # ---------------------------------------------------- unit: layout
        f = root / "u"
        _write_fasta(f / "unwrapped.fasta", APPENDED)
        _write_fasta(f / "wrapped.fasta", APPENDED, width=10)
        _write_fasta(f / "short.fasta", [(">sp|P1|X_HUMAN x", "MADQ")])
        (f / "ragged.fasta").write_text(
            ">sp|P1|X_HUMAN x\nMADQEFGHIK\nLMNP\nQRST\n", encoding="utf-8")

        check("layout: an unwrapped file has no wrap width",
              fasta_layout(f / "unwrapped.fasta").width is None)
        check("layout: a wrapped file reports its width",
              fasta_layout(f / "wrapped.fasta").width == 10)
        check("layout: uneven interior lines are ragged",
              fasta_layout(f / "ragged.fasta").ragged)
        check("layout: unwrapped output vs unwrapped template matches",
              layout_matches(fasta_layout(f / "unwrapped.fasta"),
                             fasta_layout(f / "unwrapped.fasta")) is None)
        check("layout: wrapped output vs unwrapped template does not match",
              layout_matches(fasta_layout(f / "wrapped.fasta"),
                             fasta_layout(f / "unwrapped.fasta")) is not None)
        check("layout: a file too short to reveal a width is not a mismatch",
              layout_matches(fasta_layout(f / "short.fasta"),
                             Layout(width=10, max_line=4)) is None,
              "every line fits inside the template width")
        check("layout: an unwrapped file longer than the template width is",
              layout_matches(fasta_layout(f / "unwrapped.fasta"),
                             Layout(width=10)) is not None)

        # ------------------------------------------------------ unit: locate
        clean = root / "clean"
        out, manifest, disposition, proteome = build_release(clean)
        rel = locate_release(out)
        check("locate: finds the combined FASTA",
              rel.combined is not None
              and rel.combined.name.endswith(".target.fasta"))
        check("locate: finds the decoy FASTA", rel.decoy is not None)
        check("locate: finds every per-type file", len(rel.per_type) == 5,
              f"{[p.name for p in rel.per_type]}")
        check("locate: finds the manifest and the validation report",
              rel.manifest_txt is not None and rel.validation_md is not None)

        two = root / "two" / "release"
        two.mkdir(parents=True)
        for n in ("A.target.fasta", "B.target.fasta"):
            (two / n).write_text(">x\nM\n", encoding="utf-8")
        check("locate: two releases in one folder is reported, not guessed",
              any("2 files match" in p for p in locate_release(two).problems))

        # ------------------------------------------- the nine positive cases
        base = run_checks(out, manifest, disposition, proteome)
        check("clean release: no error-severity violation",
              not fired(base, "error"), format_report(base))
        for cid, _fn in FAULTS:
            check(f"{cid} positive: the clean release does not trip {cid}",
                  cid not in fired(base))

        # ------------------------------------------- the nine negative cases
        for cid, injector in FAULTS:
            work = root / f"neg_{cid}"
            shutil.copytree(clean, work)
            args = paths_of(work)
            injector(args[0])
            vs = run_checks(*args)
            hit = fired(vs)
            check(f"{cid} negative: the injected fault trips {cid}",
                  cid in hit,
                  "; ".join(v.message for v in vs if v.check == cid)[:100])
            check(f"{cid} negative: it trips {cid} and nothing else",
                  hit == {cid}, f"fired={sorted(hit)}")

        # ------------------------------------------------- severity contract
        for cid, injector in FAULTS:
            work = root / f"sev_{cid}"
            shutil.copytree(clean, work)
            args = paths_of(work)
            injector(args[0])
            vs = run_checks(*args)
            if cid == "I7":
                check("I7 is a warning, so a mislabelled entry reports "
                      "without failing the run",
                      all(v.severity == "warning"
                          for v in vs if v.check == "I7")
                      and not fired(vs, "error"))
            else:
                check(f"{cid} is an error, so it fails the run",
                      cid in fired(vs, "error"))

        # ---------------------------- a check that cannot run must say so
        bare = root / "bare" / "release"
        bare.mkdir(parents=True)
        _write_fasta(bare / f"{NAME}.target.fasta", VARIANTS + WILDTYPES)
        vs = check_release(bare, None, None)
        skipped = {v.check for v in vs
                   if v.message.startswith("not evaluated")}
        check("a check with no evidence reports 'not evaluated', not silence",
              {"I1", "I2", "I5", "I6", "I8"} <= skipped,
              f"skipped={sorted(skipped)}")
        check("'not evaluated' is a warning, so it cannot fail a run",
              all(v.severity == "warning" for v in vs
                  if v.message.startswith("not evaluated")))

        # ------------------------- regression: PCDHGB1--PCDHGA10 fusion
        # HCC1395_FUSION_009036 tripped I7 on the real data: an
        # `in_frame_fusion` whose protein is identical to a reference
        # protein. Both breakpoints are 142 kb apart inside the clustered
        # protocadherin gamma locus, where every family member splices onto
        # the same constant exons, so the call - split-read only,
        # confidence 57.41 - is a caller artifact, not a translation fault.
        # I7 must say so; the classification logic must not change.
        pcdh_loc = "chr5:141352669::chr5:141494807"
        check("regression: intra-cluster breakpoints are named as the "
              "likely cause",
              (_intralocus_note(pcdh_loc) or "").startswith(
                  "both breakpoints on chr5, 142 kb apart"),
              str(_intralocus_note(pcdh_loc)))
        check("regression: a genuinely interchromosomal fusion gets no note",
              _intralocus_note("chr9:130854064::chr22:23290413") is None)
        check("regression: two distant breakpoints on one chromosome get no "
              "note",
              _intralocus_note("chr5:1000000::chr5:90000000") is None)
        check("regression: a small-variant LOC is not read as a fusion",
              _intralocus_note("chr1:1049980G>C") is None)

        pc = root / "pcdh"
        shutil.copytree(clean, pc)
        args = paths_of(pc)
        fusion_hdr = (
            f">vr|{NAME}_FUSION_009036|PCDHGB1--PCDHGA10_HUMAN_FUSION "
            f"Protocadherin gamma-B1 OS=Homo sapiens OX=9606 "
            f"GN=PCDHGB1--PCDHGA10 VT=FUSION CSQ=in_frame_fusion "
            f"TX=ENST00000523390.2--ENST00000398610.3 LOC={pcdh_loc} "
            f"POS=804 ALSO=REFERENCE")
        combined = (VARIANTS + [(fusion_hdr, APPENDED[0][1])]
                    + WILDTYPES + APPENDED)
        _write_fasta(args[0] / f"{NAME}.target.fasta", combined)
        _write_fasta(args[0] / "by_class" / f"{NAME}.FUSION.fasta",
                     [(fusion_hdr, APPENDED[0][1])])
        vs = [v for v in check_release(args[0], proteome=args[3])
              if v.check == "I7"]
        check("regression: I7 fires on the PCDHGB1--PCDHGA10 entry",
              len(vs) == 1 and vs[0].severity == "warning")
        check("regression: I7 names the fusion caller, not the translation",
              "points at the fusion caller" in vs[0].message,
              vs[0].message[-90:])
        check("regression: the offending id carries its breakpoint distance",
              any("142 kb apart" in e for e in vs[0].examples),
              str(vs[0].examples))

        # ------------------------------------- --skip-invariant escape hatch
        for cid, injector in FAULTS:
            work = root / f"skip_{cid}"
            shutil.copytree(clean, work)
            args = paths_of(work)
            injector(args[0])
            vs = check_release(args[0], args[1], args[2], proteome=args[3],
                               decoys_requested=True, skip_checks=[cid])
            check(f"skip: --skip-invariant {cid} stops {cid} failing the run",
                  not fired(vs, "error"),
                  "; ".join(v.message for v in vs
                            if v.severity == "error")[:90])
            check(f"skip: {cid} is still reported as turned off, not silent",
                  any(v.check == cid and v.message.startswith("skipped")
                      for v in vs))

        check("skip: a lower-case id is accepted",
              any(v.check == "I5" and v.message.startswith("skipped")
                  for v in check_release(out, skip_checks=["i5"])))
        check("skip: the summary line says checks were turned off",
              "1 check(s) turned off"
              in format_report(check_release(out, skip_checks=["I5"])))
        try:
            check_release(out, skip_checks=["I99"])
            typo = "accepted"
        except ValueError as e:
            typo = str(e)
        check("skip: a typo'd id raises rather than quietly disabling "
              "nothing", typo.startswith("unknown invariant id"), typo[:70])

        # I6 must not go quiet if the report format changes under it.
        blind = root / "blind"
        shutil.copytree(clean, blind)
        args = paths_of(blind)
        (args[0] / "qc" / "uniprot_validation.md").write_text(
            VALIDATION_MD.replace("Exact-agreement: **100.0%**.",
                                  "Agreement was fine, honestly."),
            encoding="utf-8")
        write_manifest(args[0])
        check("I6 fails loudly when it can no longer read the agreement rate",
              "I6" in fired(run_checks(*args), "error"))

    finally:
        shutil.rmtree(root, ignore_errors=True)

    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    for name in FAIL:
        print(f"  FAILED: {name}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
