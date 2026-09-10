"""v2p command line - variant sequences to protein sequences.

Three verbs instead of nine scripts:

    v2p detect <folder>                     what is in here
    v2p run    <folder> --ref <dir>         the whole conversion
    v2p audit  --ref <dir>                  check the reference interface

`detect` never writes anything, so it is safe to run on someone else's
data to see what the pipeline would make of it.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from v2p.config import (
    ConfigError, config_from_namespace, load_config,
)
from v2p.discover import (
    ANNOTATION, FUSION, GENOME, PROTEOME, REGIONS, RNA_EDITING,
    ROLE_LABEL, SMALL_VARIANTS, SPLICING, TRANSLATIONS,
    Detection, build_run_plan, detect_folder,
)
from v2p.invariants import (
    CHECK_TITLES, check_release, format_report,
)
from v2p.provenance import RunLogger

# The pipeline stages are separate executable scripts, run as subprocesses
# so each one records its own provenance. Finding them differs between a
# source checkout and an installed package, and guessing wrong produces a
# confusing "No such file" from deep inside a subprocess, so resolve it
# once and say plainly what was searched.
_SENTINEL = "01_parse_inputs.py"


def _stage_candidates() -> list[Path]:
    """Where the stage scripts might be, best source first.

    `importlib.resources` is the supported way to find a file that ships
    inside a package, and it is what makes an installed wheel work. The
    repo-layout fallback is for running straight from a checkout where the
    package has not been installed at all.
    """
    cands: list[Path] = []
    try:
        from importlib.resources import files
        cands.append(Path(str(files("v2p.stages"))))
    except (ImportError, ModuleNotFoundError, TypeError):
        # A zipimported or otherwise non-filesystem package has no usable
        # path. Fall through to the layout probes rather than failing here.
        pass
    here = Path(__file__).resolve()
    cands.append(here.parent / "stages")            # src/v2p/stages
    cands.append(here.parents[2] / "scripts")       # pre-1.0 checkout layout
    return cands


def stage_dir() -> Path:
    """Directory holding the numbered pipeline stage scripts."""
    cands = _stage_candidates()
    for cand in cands:
        if (cand / _SENTINEL).is_file():
            return cand
    raise FileNotFoundError(
        "cannot find the v2p pipeline stage scripts. Searched:\n  "
        + "\n  ".join(str(c) for c in cands)
        + f"\n\nExpected to find {_SENTINEL} in one of these. If you are "
          "running from a source checkout, run from the repository root; if "
          "from an installed package, the install is incomplete - reinstall "
          "with `pip install v2p`.")


def _run(script: str, args: list[str], quiet: bool = False,
         logdir: str | None = None) -> None:
    # Every stage takes --logdir and defaults it to logs/. Passing the run's
    # own value through keeps one run's provenance in one place; without it,
    # `v2p run --logdir X` scattered the stage JSONs into logs/ and left
    # only the invariant record in X.
    if logdir is not None:
        args = [*args, "--logdir", logdir]
    cmd = [sys.executable, str(stage_dir() / script), *args]
    if quiet:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    else:
        subprocess.run(cmd, check=True)


# --------------------------------------------------------------------------

def cmd_detect(args: argparse.Namespace) -> int:
    dets = detect_folder(args.folder, recursive=not args.no_recursive)
    if getattr(args, "ref", None):
        # Reference files normally live outside the input folder. Without
        # scanning them, detect reports a missing genome that is not
        # actually missing, and a warning nobody can act on is a warning
        # people learn to ignore.
        dets += detect_folder(args.ref, recursive=False)
    plan = build_run_plan(dets)

    if args.json:
        print(json.dumps({
            "detections": [{"path": str(d.path), "role": d.role,
                            "confidence": d.confidence,
                            "records": d.n_records,
                            "evidence": d.evidence, "detail": d.detail}
                           for d in dets],
            "plan": plan,
        }, indent=2))
        return 1 if plan["problems"] else 0

    print(f"\n{'file':<38} {'identified as':<40} conf  records")
    print("-" * 96)
    for d in sorted(dets, key=lambda x: (x.role, x.path.name)):
        n = f"{d.n_records:,}" if d.n_records else ""
        print(f"{d.path.name[:38]:<38} {ROLE_LABEL[d.role][:40]:<40} "
              f"{d.confidence:>4.0%}  {n:>9}")
        if args.verbose:
            for e in d.evidence:
                print(f"{'':<38}   - {e}")

    print("\nrun plan")
    print("-" * 96)
    for role, path in plan["selected"].items():
        print(f"  {ROLE_LABEL[role]:<40} {Path(path).name}")
        for alt in plan["alternatives"].get(role, []):
            print(f"  {'':<40}   (also: {Path(alt).name})")
    if plan.get("genome_build"):
        print(f"  {'genome build':<40} {plan['genome_build']}")

    if plan["problems"]:
        print("\nproblems")
        print("-" * 96)
        for p in plan["problems"]:
            print(f"  ! {p}")
        return 1
    print("\nno problems. `v2p run` would convert "
          + ", ".join(ROLE_LABEL[e] for e in plan["evidence_types"]) + ".")
    return 0


# Config keys that name one file each, and the role each one fills. A
# path given here is an assertion by the user, so it replaces detection
# for that role rather than competing with it.
_EXPLICIT_PATHS = (
    ("genome", GENOME), ("annotation", ANNOTATION),
    ("proteome", PROTEOME), ("translations", TRANSLATIONS),
    ("small_variants", SMALL_VARIANTS), ("rna_editing", RNA_EDITING),
    ("fusion_calls", FUSION), ("splicing", SPLICING),
)


def _apply_explicit_paths(args: argparse.Namespace,
                          dets: list) -> list | None:
    """Replace detections for any role the user named a file for.

    Returns the amended detection list, or None if a named file is not
    there - which is a hard stop, because carrying on would build a
    release from whatever detection happened to find instead.

    Substituting into the detection list rather than into the finished
    plan means `build_run_plan` still decides everything: an explicitly
    named file is the only candidate for its role, so the ambiguity
    warning and the "no reference genome found" problem both disappear
    on their own instead of having to be filtered out afterwards.
    """
    named: dict[str, Path] = {}
    for attr, role in _EXPLICIT_PATHS:
        val = getattr(args, attr, None)
        if val:
            named[role] = Path(val)
    if not named:
        return dets

    for role, p in named.items():
        if not p.exists():
            print(f"cannot run - {ROLE_LABEL[role]}: no such file: {p}",
                  file=sys.stderr)
            return None

    out = [d for d in dets if d.role not in named]
    for role, p in sorted(named.items()):
        # Reuse the real detection when the named file was also found by
        # scanning, so its record count and genome build survive.
        prior = next((d for d in dets
                      if d.role == role and d.path == p), None)
        out.append(prior or Detection(
            path=p, role=role, confidence=1.0,
            evidence=["named explicitly, so detection was not consulted"]))
    return out


def cmd_run(args: argparse.Namespace) -> int:
    dets: list = []
    if getattr(args, "folder", None):
        dets = detect_folder(args.folder, recursive=not args.no_recursive)
    if getattr(args, "ref", None):
        # Same reason as in `detect`: the reference files live outside the
        # input folder, so a plan built from the input folder alone reports
        # the genome and annotation as missing - the very files --ref
        # supplies. Without this the documented invocation
        # `v2p run <folder> --ref <dir>` refuses to run without --force.
        dets += detect_folder(args.ref, recursive=False)
    dets = _apply_explicit_paths(args, dets)
    if dets is None:
        return 2
    if not dets:
        print("cannot run - no input folder and no explicit input paths; "
              "give a folder, or name the files in a --config file",
              file=sys.stderr)
        return 2
    plan = build_run_plan(dets)
    sel = plan["selected"]

    ref = Path(args.ref) if args.ref else None
    if ref:
        for role, pats in ((GENOME, ("*.fa", "*.fasta")),
                           (ANNOTATION, ("*.gtf", "*.gtf.gz", "*.gff*")),
                           (PROTEOME, ("*uniprot*", "*sprot*", "*.fasta")),
                           (TRANSLATIONS, ("*translations*",))):
            if role in sel:
                continue
            for pat in pats:
                hits = sorted(ref.glob(pat))
                if hits:
                    d = detect_folder(hits[0].parent)
                    for x in d:
                        if x.path == hits[0] and x.role == role:
                            sel[role] = str(hits[0])
                            break
                if role in sel:
                    break

    missing = [r for r in (GENOME, ANNOTATION, PROTEOME) if r not in sel]
    if missing:
        print("cannot run - missing: "
              + ", ".join(ROLE_LABEL[m] for m in missing), file=sys.stderr)
        print("supply them with --ref <dir>", file=sys.stderr)
        return 2
    if not plan["evidence_types"]:
        print("cannot run - no variant evidence in " + str(args.folder),
              file=sys.stderr)
        return 2
    if plan["problems"] and not args.force:
        for p in plan["problems"]:
            print(f"! {p}", file=sys.stderr)
        print("refusing to run; re-run with --force to override",
              file=sys.stderr)
        return 2

    out = Path(args.outdir)
    work = out / "_work"
    work.mkdir(parents=True, exist_ok=True)
    print(f"converting {', '.join(ROLE_LABEL[e] for e in plan['evidence_types'])}")

    # -- stage 1 ---------------------------------------------------------
    a = ["--outdir", str(work)]
    for role, flag in ((SMALL_VARIANTS, "--vcf"), (RNA_EDITING, "--res"),
                       (FUSION, "--fusion"), (SPLICING, "--as-lr")):
        if role in sel:
            a += [flag, sel[role]]
    print("  [1/5] parsing inputs")
    _run("01_parse_inputs.py", a, quiet=True, logdir=args.logdir)

    manifest = work / "tables" / "unified_variant_manifest.tsv"
    print("  [2/5] input QC")
    _run("03_qc_report.py", ["--manifest", str(manifest),
                             "--outdir", str(work)], quiet=True,
         logdir=args.logdir)

    print("  [3/5] translating")
    b = ["--manifest", str(manifest),
         "--genome", sel[GENOME], "--gtf", sel[ANNOTATION],
         "--uniprot", sel[PROTEOME],
         "--header-style", args.header_style,
         "--transcript-mode", args.transcript_mode,
         "--include-reference",
         "--species", args.species,
         *(["--include-noncanonical"] if args.include_noncanonical else []),
         "--nc-min-aa", str(args.nc_min_aa),
         *(["--combine-variants"] if args.combine_variants else []),
         *([] if args.allow_unphased else ["--no-allow-unphased"]),
         "--combine-max", str(args.combine_max),
         "--emit-disposition", str(work / "tables" / "disposition.tsv"),
         "--outdir", str(work)]
    if args.keep_unchanged:
        b.append("--keep-synonymous")
    _run("02_build_protein_fasta.py", b, quiet=True, logdir=args.logdir)

    built = (work / "fasta" /
             f"HCC1395_variant_proteins.{args.header_style}."
             f"{args.transcript_mode}.fasta")
    if not built.exists():
        cands = sorted((work / "fasta").glob("*.fasta"))
        if not cands:
            print("translation produced no FASTA", file=sys.stderr)
            return 1
        built = cands[0]

    print("  [4/5] validating")
    v = ["--uniprot", sel[PROTEOME],
         "--protein-fasta", str(built), "--outdir", str(work),
         "--min-identity-rate", str(args.min_agreement)]
    rec = work / "tables" / "res_recoding_sites.tsv"
    if rec.exists():
        v += ["--recoding", str(rec)]
    if TRANSLATIONS in sel:
        v += ["--gencode-translations", sel[TRANSLATIONS]]
    try:
        _run("04_validate_uniprot.py", v, quiet=True, logdir=args.logdir)
    except subprocess.CalledProcessError:
        print("  ! validation below threshold - see "
              f"{work}/qc/uniprot_validation.md", file=sys.stderr)
        if not args.force:
            return 1
    _run("08_recovery_report.py",
         ["--manifest", str(manifest),
          "--disposition", str(work / "tables" / "disposition.tsv"),
          "--outdir", str(work)], quiet=True, logdir=args.logdir)

    print("  [5/5] packaging")
    pk = ["--fasta", str(built), "--uniprot", sel[PROTEOME],
          "--name", args.name, "--outdir", str(out),
          "--extra-file", str(work / "qc")]
    if args.decoys != "none":
        pk += ["--decoy", args.decoys]
    if args.split_by_type:
        pk.append("--split-by-class")
    if args.append_reference:
        pk.append("--append-reference")
    _run("06_package_release.py", pk, quiet=True, logdir=args.logdir)

    # -- release invariants ---------------------------------------------
    # Run after packaging and before declaring success. Every check reads;
    # none writes, so MANIFEST.txt stays valid for the checksum check. The
    # provenance JSON lands in --logdir, outside the release, for the same
    # reason.
    print("  checking release invariants")
    rl = RunLogger("v2p_run_invariants", args.logdir)
    skips = sorted(set(args.skip_invariant or []))
    rl.add_params(outdir=str(out), min_agreement=args.min_agreement,
                  decoys=args.decoys, skip_invariant=skips)
    # The whole resolved configuration, every key, defaulted ones
    # included - so the provenance record answers "what settings produced
    # this?" without the reader having to know what the defaults were at
    # the time. This is the same object --write-config emits, so a run
    # record can be turned back into a runnable config file.
    rl.add_params(config=config_from_namespace(args).to_dict())
    for cid in skips:
        # Recorded per check as well as in params, so grepping a log for a
        # disabled check finds it.
        rl.log.warning("invariant %s turned off by --skip-invariant", cid)
        rl.count(f"skipped.{cid}")

    try:
        violations = check_release(
            out,
            work / "tables" / "unified_variant_manifest.tsv",
            work / "tables" / "disposition.tsv",
            proteome=Path(sel[PROTEOME]),
            min_agreement=args.min_agreement,
            decoys_requested=args.decoys != "none",
            skip_checks=skips)
    except ValueError as e:
        rl.log.error("%s", e)
        rl.close(status="error")
        print(f"\n{e}", file=sys.stderr)
        return 2

    print()
    print(format_report(violations))
    for v in violations:
        rl.count(f"{v.severity}.{v.check}")
    rl.add_params(violations=[{"check": v.check, "severity": v.severity,
                               "message": v.message, "examples": v.examples}
                              for v in violations])

    errors = [v for v in violations if v.severity == "error"]
    rl.close(status="error" if errors else "ok")
    if errors:
        # --force overrides detection problems and the validation
        # threshold. It deliberately does not override these: a release
        # that contradicts itself should not be shipped on a flag. The
        # per-check escape hatch is --skip-invariant, which is recorded.
        print(f"\n{len(errors)} invariant violation(s) - this release is not "
              f"internally consistent", file=sys.stderr)
        return 1

    print(f"\ndone -> {out}/")
    fa = sorted(out.glob("*.target.fasta"))
    if fa:
        n = sum(1 for line in open(fa[0]) if line.startswith(">"))
        print(f"  {fa[0].name}: {n:,} sequences")
    return 0


def cmd_audit(args: argparse.Namespace) -> int:
    ref = Path(args.ref)
    a = []
    for role, flag in ((ANNOTATION, "--gtf"), (GENOME, "--genome"),
                       (PROTEOME, "--uniprot"),
                       (TRANSLATIONS, "--gencode-translations")):
        for d in detect_folder(ref, recursive=False):
            if d.role == role:
                a += [flag, str(d.path)]
                break
    if not a:
        print(f"no reference files recognised in {ref}", file=sys.stderr)
        return 2
    try:
        _run("00_audit_references.py", a)
    except subprocess.CalledProcessError as e:
        return e.returncode
    return 0


def _add_run_arguments(r: argparse.ArgumentParser,
                       suppress: bool = False) -> None:
    """Declare `run`'s arguments.

    With ``suppress=True`` every default becomes ``argparse.SUPPRESS``, so
    parsing the same argv a second time yields a namespace holding only
    what the user actually typed. That is what makes the precedence rule
    - defaults, then the config file, then the command line - decidable:
    without it there is no way to tell `--decoys pseudo_reverse` from the
    default of the same name, and a config file would be unable to change
    any option that has a default.
    """
    def A(*names, **kw):
        if suppress:
            kw["default"] = argparse.SUPPRESS
        r.add_argument(*names, **kw)

    A("folder", nargs="?",
      help="folder of variant calls. May be omitted when a --config file "
           "names the inputs.")
    A("--ref", help="directory holding genome, annotation, proteome")
    A("--outdir", default="v2p_output")
    A("--name", default="variant_proteome")
    A("--transcript-mode", default="all",
      choices=["all", "representative"])
    A("--combine-variants", action="store_true",
      help="also emit one protein per haplotype carrying two or more "
           "co-occurring variants. A tryptic peptide spanning two variants "
           "exists only in the combined form, so a single-variant database "
           "cannot identify it at any FDR; combinations yielding no such "
           "peptide are dropped rather than inflating the database. Phase "
           "is honoured where the caller reports GT/PS, and variants on "
           "opposite haplotypes are never combined; without phase the entry "
           "is a hypothesis, marked unphased. OFF by default.")
    A("--combine-max", type=int, default=8,
      help="most variants combined on one transcript (default 8). Above "
           "this a transcript is far likelier to be an alignment artefact "
           "than a real proteoform.")
    A("--no-allow-unphased", dest="allow_unphased", action="store_false",
      help="refuse to combine variants the caller did not phase. Default "
           "is to allow them, since sites-only and unphased VCFs are the "
           "common case, but every entry records PHASE=phased or "
           "PHASE=unphased either way, so hypotheses can be filtered from "
           "the database without rebuilding it.")
    A("--include-noncanonical", action="store_true",
      help="also three-frame translate non-coding transcripts (lncRNA, "
           "pseudogene) into NC_* entries. OFF by default: it multiplies "
           "database size several-fold, and an inflated search space costs "
           "sensitivity at a fixed FDR.")
    A("--nc-min-aa", type=int, default=30,
      help="minimum non-canonical ORF length in residues (default 30)")
    A("--species", default="human",
      help="species name (human, mouse, ...) or a path to a "
           "config/species/*.yaml. Sets the entry-name suffix and the "
           "OS=/OX= fields in every header. Default human.")
    A("--header-style", default="uniprot",
      choices=["uniprot", "peff", "pvac", "descriptive"])
    A("--keep-unchanged", action="store_true", default=True,
      help="keep synonymous and UTR variants (default)")
    A("--drop-unchanged", dest="keep_unchanged", action="store_false")
    A("--decoys", default="pseudo_reverse",
      choices=["none", "reverse", "pseudo_reverse", "shuffle"])
    A("--split-by-type", action="store_true", default=True)
    A("--no-split", dest="split_by_type", action="store_false")
    A("--append-reference", action="store_true", default=True)
    A("--no-reference", dest="append_reference", action="store_false")
    A("--no-recursive", action="store_true")
    A("--min-agreement", type=float, default=0.90,
      help="minimum agreement between our reference "
           "translations and the annotation source's own "
           "(default 0.90). Governs both the validation gate "
           "and invariant I6, so the two cannot drift apart.")
    A("--logdir", default="logs",
      help="where stage logs and provenance JSONs are written "
           "(default logs/). Keep this outside --outdir, or the "
           "release manifest goes stale as soon as it is "
           "written.")
    A("--skip-invariant", action="append", metavar="ID",
      choices=sorted(CHECK_TITLES),
      help="turn off one release invariant; repeatable. For a "
           "check that is wrong about your data - I5 fires when "
           "the supplied proteome is wrapped and you want "
           "unwrapped output. Every skip is recorded in the "
           "provenance JSON and printed in the report, so a "
           "release built with a check off cannot be mistaken "
           "for one that passed it.")
    A("--force", action="store_true",
      help="proceed despite detection problems or a failed "
           "validation gate. Does not override the release "
           "invariants - use --skip-invariant for that.")

    # -- reproducibility -------------------------------------------------
    A("--config", metavar="FILE",
      help="run configuration in YAML, per docs/TOOL_DESIGN.md. Anything "
           "also given on the command line wins over the file.")
    A("--write-config", metavar="FILE",
      help="write the configuration this command line implies to FILE and "
           "exit without running. Round-trips: feeding the result back "
           "through --config reproduces the same settings.")

    # Naming a file settles what detection would otherwise have to guess,
    # and lets a config file be a complete description of a run.
    A("--genome", help="reference genome FASTA, overriding detection")
    A("--annotation", help="GTF/GFF annotation, overriding detection")
    A("--proteome", help="reference proteome FASTA, overriding detection")
    A("--translations",
      help="annotation source's own protein translations, for invariant I6")
    A("--small-variants", help="VCF of small variants, overriding detection")
    A("--rna-editing", help="RNA editing table, overriding detection")
    A("--fusion-calls", help="gene fusion calls, overriding detection")
    A("--splicing", help="alternative splicing events, overriding detection")


def _build_parser(suppress: bool = False) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="v2p", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("detect", help="identify the contents of a folder")
    d.add_argument("folder")
    d.add_argument("--json", action="store_true")
    d.add_argument("--verbose", "-v", action="store_true",
                   help="show the evidence for each decision")
    d.add_argument("--ref", help="also scan this directory for reference files")
    d.add_argument("--no-recursive", action="store_true")
    if not suppress:
        d.set_defaults(func=cmd_detect)

    r = sub.add_parser("run", help="convert a folder of variant calls")
    _add_run_arguments(r, suppress=suppress)
    if not suppress:
        r.set_defaults(func=cmd_run)

    a = sub.add_parser("audit", help="check how the reference files are read")
    a.add_argument("--ref", required=True)
    if not suppress:
        a.set_defaults(func=cmd_audit)
    return ap


def _merge_config(ns: argparse.Namespace,
                  argv: list[str] | None) -> argparse.Namespace:
    """Fold a --config file into an already-parsed `run` namespace.

    Precedence is defaults, then the file, then the command line. A key
    the file does not mention leaves the default alone, which is why
    RunConfig records which keys were actually present.
    """
    if not getattr(ns, "config", None):
        return ns
    typed = set(vars(_build_parser(suppress=True).parse_args(argv)))
    cfg = load_config(ns.config)
    for dest, value in cfg.overrides().items():
        if dest in typed:
            continue
        setattr(ns, dest, value)
    return ns


def main(argv: list[str] | None = None) -> int:
    ns = _build_parser().parse_args(argv)

    if getattr(ns, "cmd", None) == "run":
        try:
            ns = _merge_config(ns, argv)
        except ConfigError as e:
            print(f"config error: {e}", file=sys.stderr)
            return 2
        dest = getattr(ns, "write_config", None)
        if dest:
            p = Path(dest)
            if p.parent != Path(""):
                p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(config_from_namespace(ns).to_yaml(), encoding="utf-8")
            print(f"wrote {p}")
            return 0

    return ns.func(ns)


if __name__ == "__main__":
    raise SystemExit(main())
