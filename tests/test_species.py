#!/usr/bin/env python3
"""Milestone 4 acceptance test - species and annotation independence.

Run:  python tests/test_species.py
Exit code 0 = all pass. No pytest dependency so it runs anywhere.

M4 is a refactor, so the load-bearing assertion is the negative one: with
no species given, every header is byte-for-byte what it was when the
species was hard-coded. The mouse cases prove the parameter is real and
not decorative.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2p.build.smallvar import ProteinRecord                # noqa: E402
from v2p.fasta import HEADER_STYLES                         # noqa: E402
from v2p.species import (                                   # noqa: E402
    HUMAN, MOUSE, available, from_mapping, load_species,
)

PASS: list[str] = []
FAIL: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    (PASS if cond else FAIL).append(name)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}{' | ' + detail if detail else ''}")


def raises(fn, want: str = "") -> tuple[bool, str]:
    try:
        fn()
    except ValueError as e:
        return (want.lower() in str(e).lower()), str(e)
    except Exception as e:                       # noqa: BLE001
        return False, f"raised {type(e).__name__}: {e}"
    return False, "did not raise"


def record(**kw) -> ProteinRecord:
    base = dict(seq_id="X1", variant_class="SNV", consequence="missense",
                gene="tp53", transcript="ENST1", protein_change="p.R1G",
                locus="chr1:100A>G", sequence="MKV", source="vcf",
                confidence="PASS")
    base.update(kw)
    return ProteinRecord(**base)


def main() -> int:
    # ------------------------------------------------------- the defaults
    check("the default species is human",
          load_species(None) == HUMAN and load_species("") == HUMAN)
    check("human and mouse are both available",
          {"human", "mouse"} <= set(available()), str(available()))
    check("mouse carries its own taxon and suffix",
          MOUSE.taxon_id == 10090 and MOUSE.entry_suffix == "MOUSE"
          and MOUSE.assembly == "GRCm39",
          f"{MOUSE.taxon_id}/{MOUSE.entry_suffix}/{MOUSE.assembly}")
    check("entry_name appends the suffix",
          HUMAN.entry_name("TP53") == "TP53_HUMAN"
          and MOUSE.entry_name("Trp53") == "Trp53_MOUSE")
    check("os_ox renders the UniProt fragment",
          HUMAN.os_ox() == "OS=Homo sapiens OX=9606", HUMAN.os_ox())

    # ------------------------------------ the shipped YAML matches the code
    # If config/species/human.yaml drifted from v2p.species.HUMAN, a run
    # would behave differently depending on how the species was named.
    for name, want in (("human", HUMAN), ("mouse", MOUSE)):
        f = ROOT / "config" / "species" / f"{name}.yaml"
        if not f.is_file():
            check(f"config/species/{name}.yaml exists", False, "missing")
            continue
        got = load_species(f)
        check(f"config/species/{name}.yaml matches the built-in {name}",
              got == want,
              f"{got.taxon_id}/{got.entry_suffix} vs "
              f"{want.taxon_id}/{want.entry_suffix}")

    check("a species can be loaded by bare name",
          load_species("mouse").taxon_id == 10090)

    # ------------------------------------------------------- the refactor
    # The negative case that makes M4 a refactor rather than a change.
    r = record()
    for style in sorted(HEADER_STYLES):
        fmt = HEADER_STYLES[style]
        check(f"{style}: omitting species equals passing human explicitly",
              fmt(r, "P") == fmt(r, "P", HUMAN))

    uni = HEADER_STYLES["uniprot"]
    check("human variant header still says Homo sapiens OX=9606",
          "OS=Homo sapiens OX=9606" in uni(r, "P")
          and "TP53_HUMAN_SNV" in uni(r, "P"),
          uni(r, "P")[:90])

    # ------------------------------------------------------- mouse output
    m = uni(r, "P", MOUSE)
    check("mouse variant header carries OX=10090",
          "OS=Mus musculus OX=10090" in m, m[:90])
    check("mouse variant entry name uses _MOUSE",
          "TP53_MOUSE_SNV" in m, m[:60])
    check("no human marker survives in mouse output",
          "_HUMAN" not in m and "9606" not in m and "Homo" not in m)

    ref = record(variant_class="REFERENCE", consequence="reference")
    mr = uni(ref, "P", MOUSE)
    check("mouse reference header uses _MOUSE and OX=10090",
          "TP53_MOUSE" in mr and "OX=10090" in mr, mr[:90])

    peff = HEADER_STYLES["peff"]
    check("peff carries the mouse taxonomy keys",
          "\\TaxName=Mus musculus" in peff(r, "P", MOUSE)
          and "\\NcbiTaxId=10090" in peff(r, "P", MOUSE))

    # --------------------------------------------- loading from a mapping
    fly = from_mapping({"scientific_name": "Drosophila melanogaster",
                        "taxon_id": 7227, "entry_suffix": "drome",
                        "assembly": "BDGP6"})
    check("a third species needs no code change",
          fly.entry_name("Act5C") == "Act5C_DROME"
          and fly.os_ox() == "OS=Drosophila melanogaster OX=7227",
          fly.os_ox())
    # The audit must degrade, not fail, when it has nothing to assert.
    check("a species with no contig lengths is legal",
          fly.contig_lengths == {} and bool(HUMAN.contig_lengths))

    # ---------------------------------------------------------- rejection
    ok, msg = raises(lambda: from_mapping({"taxon-id": 7227}), "taxon-id")
    check("a misspelled key is refused, not ignored", ok, msg[:70])
    ok, msg = raises(lambda: from_mapping({"taxon_id": "not a number"}),
                     "whole number")
    check("a non-numeric taxon id is refused", ok, msg[:70])
    ok, msg = raises(lambda: from_mapping({"contig_lengths": "chr1"}),
                     "mapping")
    check("contig_lengths must be a mapping", ok, msg[:70])
    ok, msg = raises(lambda: load_species("nosuchbeast"), "unknown species")
    check("an unknown species name is refused", ok, msg[:70])
    ok, msg = raises(lambda: load_species("/no/such/file.yaml"), "no species")
    check("a missing species file is refused", ok, msg[:70])

    tmp = Path(tempfile.mkdtemp())
    y = tmp / "zfish.yaml"
    y.write_text("scientific_name: Danio rerio\ntaxon_id: 7955\n"
                 "entry_suffix: DANRE\nassembly: GRCz11\n"
                 "contig_lengths:\n  chr1: 59578282\n", encoding="utf-8")
    z = load_species(y)
    check("a species loads from an arbitrary YAML path",
          z.taxon_id == 7955 and z.entry_name("actb") == "actb_DANRE"
          and z.contig_lengths == {"chr1": 59578282},
          f"{z.taxon_id}/{z.contig_lengths}")

    check("mitochondrial contigs and code are per species",
          HUMAN.is_mito("chrM") and HUMAN.mito_code == 2
          and not HUMAN.is_mito("chr1"))

    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
