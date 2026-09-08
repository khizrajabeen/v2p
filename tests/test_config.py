#!/usr/bin/env python3
"""Milestone 2 acceptance test - configuration and reproducibility.

Run:  python tests/test_config.py
Exit code 0 = all pass. No pytest dependency so it runs anywhere.

The reproducibility half of M2 - two runs from one config producing a
byte-identical FASTA - needs the 18 GB reference, so it lives in
`make reproducibility` rather than here. What *is* here is the property
that makes it possible: a config file and a command line resolve to one
set of values, deterministically, with no key silently ignored.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2p.config import (                                   # noqa: E402
    SCHEMA, ConfigError, config_from_namespace, default_config,
    load_config, loads_config,
)

PASS: list[str] = []
FAIL: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    (PASS if cond else FAIL).append(name)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}{' | ' + detail if detail else ''}")


def raises(fn, want: str = "") -> tuple[bool, str]:
    """True if `fn` raises ConfigError, with `want` in the message."""
    try:
        fn()
    except ConfigError as e:
        return (want.lower() in str(e).lower()), str(e)
    except Exception as e:                       # noqa: BLE001
        return False, f"raised {type(e).__name__} not ConfigError: {e}"
    return False, "did not raise"


def main() -> int:
    # ------------------------------------------------------------ defaults
    d = default_config()
    check("every schema section has defaults",
          set(d) == set(SCHEMA), f"{sorted(d)}")
    check("a default is the documented one",
          d["translate"]["transcript_mode"] == "all"
          and d["output"]["header_style"] == "uniprot",
          f"{d['translate']['transcript_mode']}/{d['output']['header_style']}")
    check("an empty config is legal and yields pure defaults",
          loads_config("").to_dict() == d)
    check("a comment-only config is legal",
          loads_config("# nothing here\n").to_dict() == d)

    # ------------------------------------------------- accepted overrides
    cfg = loads_config(
        "translate:\n"
        "  transcript_mode: representative\n"
        "  keep_unchanged: false\n"
        "output:\n"
        "  name: mine\n"
        "  decoys: none\n"
        "validate:\n"
        "  min_translation_agreement: 0.75\n")
    check("a string key is read",
          cfg.get("translate", "transcript_mode") == "representative")
    check("a bool key is read as a bool",
          cfg.get("translate", "keep_unchanged") is False)
    check("a float key is read as a float",
          cfg.get("validate", "min_translation_agreement") == 0.75)
    check("an untouched key keeps its default",
          cfg.get("output", "header_style") == "uniprot")

    # `provided` is what makes precedence decidable: a key the file did not
    # mention must not silently beat a command-line option.
    check("only the keys the file set are marked provided",
          cfg.overrides() == {"transcript_mode": "representative",
                              "keep_unchanged": False, "name": "mine",
                              "decoys": "none", "min_agreement": 0.75},
          str(sorted(cfg.overrides())))
    check("a key set to its own default is still marked provided",
          "transcript_mode" in loads_config(
              "translate:\n  transcript_mode: all\n").overrides())

    # --------------------------------------------------------- rejections
    # Every one of these is a config that would otherwise run and quietly
    # not mean what it says.
    ok, msg = raises(lambda: loads_config("nosuchsection:\n  a: 1\n"),
                     "nosuchsection")
    check("an unknown section is refused", ok, msg[:70])

    ok, msg = raises(lambda: loads_config("translate:\n  nosuchkey: 1\n"),
                     "nosuchkey")
    check("an unknown key is refused", ok, msg[:70])

    ok, msg = raises(
        lambda: loads_config("translate:\n  transcript_mode: representive\n"),
        "representative")
    check("a near-miss value suggests the right one", ok, msg[:80])

    ok, msg = raises(
        lambda: loads_config("translate:\n  transcript_mode: sideways\n"))
    check("a value outside the allowed set is refused", ok, msg[:70])

    ok, msg = raises(
        lambda: loads_config("translate:\n  keep_unchanged: maybe\n"))
    check("a non-bool for a bool key is refused", ok, msg[:70])

    ok, msg = raises(
        lambda: loads_config(
            "validate:\n  min_translation_agreement: quite high\n"))
    check("a non-number for a float key is refused", ok, msg[:70])

    ok, msg = raises(lambda: loads_config("translate: 3\n"))
    check("a section that is not a mapping is refused", ok, msg[:70])

    ok, msg = raises(lambda: load_config(ROOT / "no" / "such" / "file.yaml"))
    check("a missing config file is refused", ok, msg[:70])

    # ------------------------------------------------------- round trips
    # --write-config then --config must be a no-op. If it is not, a run
    # record cannot be turned back into the run.
    ns = argparse.Namespace(
        folder="data", ref="ref", outdir="out", name="RT",
        transcript_mode="representative", header_style="peff",
        keep_unchanged=False, decoys="shuffle", split_by_type=False,
        append_reference=False, logdir="lg", min_agreement=0.8,
        skip_invariant=["I5"], force=True, genetic_code="auto",
        species="mouse", include_noncanonical=True, nc_min_aa=45,
        nc_any_start=True,
        small_variants=None, rna_editing=None, fusion_calls=None,
        splicing=None, genome=None, annotation=None, proteome=None,
        translations=None)
    a = config_from_namespace(ns)
    b = loads_config(a.to_yaml())
    diffs = "; ".join(
        f"{s}.{k}: {a.to_dict()[s][k]!r} != {b.to_dict()[s][k]!r}"
        for s in a.to_dict() for k in a.to_dict()[s]
        if a.to_dict()[s][k] != b.to_dict()[s][k])
    check("a command line round-trips through YAML unchanged",
          a.to_dict() == b.to_dict(), diffs[:110])
    check("the emitted YAML is byte-stable across two emissions",
          a.to_yaml() == config_from_namespace(ns).to_yaml())
    check("a list key survives the round trip",
          b.get("validate", "skip_invariant") == ["I5"],
          repr(b.get("validate", "skip_invariant")))
    check("the species round-trips",
          b.get("translate", "species") == "mouse",
          str(b.get("translate", "species")))
    check("the non-canonical flags round-trip, including the int",
          b.get("translate", "include_noncanonical") is True
          and b.get("translate", "noncanonical_min_aa") == 45
          and b.get("translate", "noncanonical_any_start") is True,
          str(b.to_dict()["translate"]))
    ok, msg = raises(
        lambda: loads_config(
            "translate:\n  noncanonical_min_aa: thirty\n"),
        "whole number")
    check("a non-integer ORF length is refused", ok, msg[:70])
    check("a false bool survives the round trip and is not read as absent",
          b.get("translate", "keep_unchanged") is False
          and b.get("output", "split_by_type") is False)

    tmp = Path(tempfile.mkdtemp())
    p = tmp / "run.yaml"
    p.write_text(a.to_yaml(), encoding="utf-8")
    check("a written config loads back from disk",
          load_config(p).to_dict() == a.to_dict())
    check("the loaded config records where it came from",
          str(load_config(p).source or "").endswith("run.yaml"),
          str(load_config(p).source))

    # ------------------------------------------------ provenance coverage
    # BUILD_SPEC M2: "every config key appears in the provenance JSON".
    # to_dict() is what the CLI hands to RunLogger.add_params, so the test
    # is that it covers the schema exactly.
    flat = {f"{s}.{k}" for s, vals in a.to_dict().items() for k in vals}
    want = {f"{s}.{k.name}" for s, keys in SCHEMA.items() for k in keys}
    check("every config key reaches the provenance record",
          flat == want, f"missing {sorted(want - flat)}")

    # ------------------------------------------------- the shipped default
    # config/params.yaml used to be a dead file with a false header comment
    # and defaults that contradicted the code. Whatever ships must load.
    shipped = ROOT / "config" / "params.yaml"
    if not shipped.exists():
        check("no orphaned config/params.yaml is shipped", True, "absent")
    else:
        try:
            got = load_config(shipped)
            check("the shipped config/params.yaml actually loads", True,
                  f"{len(got.overrides())} keys set")
            check("the shipped config agrees with the code's defaults",
                  got.get("translate", "transcript_mode")
                  == default_config()["translate"]["transcript_mode"],
                  f"file={got.get('translate', 'transcript_mode')} "
                  f"code={default_config()['translate']['transcript_mode']}")
        except ConfigError as e:
            check("the shipped config/params.yaml actually loads", False,
                  str(e)[:90])

    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
