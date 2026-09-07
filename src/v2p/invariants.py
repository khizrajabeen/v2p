"""Mechanical checks on a finished release, run before declaring success.

Every bug that reached a shipped file during development produced
plausible-looking *wrong output*, not an error: double-counted wild-types,
per-type files that did not sum to the combined file, UTR variants labelled
frameshift, duplicated accessions, a stale manifest. None was caught by
looking at the sequences, because the sequences looked fine.

These nine checks are the ones that would have caught them. They are
deliberately mechanical - they compare a release against itself and against
the inputs it claims to describe, and they know no biology. A check that
needs judgement belongs in `validate.py`, not here.

The checks never write to the release. `MANIFEST.txt` is written as the last
step of packaging, so anything written afterwards would make I8 fail against
a release that is in fact intact.

A check that cannot be evaluated - because the artefact it needs was not
produced by this run - reports a `warning` saying so, rather than passing
silently. Silence is where wrong answers hide.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from .provenance import sha256

__all__ = [
    "Violation", "check_release", "format_report",
    "SEQUENCE_CHANGING", "CHECK_TITLES",
]

# Consequences that change the amino-acid sequence by definition. An entry
# carrying one of these whose protein equals its own reference is either
# mislabelled or should not have been emitted. Splicing form labels are
# deliberately absent: an alternative form that is itself an annotated
# isoform legitimately equals that isoform's reference protein.
SEQUENCE_CHANGING = frozenset({
    "missense", "stop_gained", "stop_lost", "start_lost",
    "inframe_insertion", "inframe_deletion", "frameshift",
    "in_frame_fusion", "frameshift_fusion",
})

CHECK_TITLES = {
    "I1": "per-type files partition the combined file",
    "I2": "every input variant appears exactly once in the disposition",
    "I3": "variant + wild-type + reference counts sum to the combined file",
    "I4": "every sequence id is unique",
    "I5": "sequence line layout matches the supplied proteome template",
    "I6": "reference translations agree with the annotation source",
    "I7": "no sequence-changing entry equals its own reference",
    "I8": "MANIFEST.txt checksums match the files on disk",
    "I9": "decoy count equals target count",
}

ERROR = "error"
WARNING = "warning"

_MAX_EXAMPLES = 5
_NOT_EVALUATED = "not evaluated"
_SKIPPED = "skipped"


@dataclass
class Violation:
    check: str
    severity: str        # "error" | "warning"
    message: str
    examples: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        s = f"[{self.severity}] {self.check}: {self.message}"
        if self.examples:
            s += "\n      e.g. " + ", ".join(self.examples)
        return s


def _v(check: str, severity: str, message: str, examples=()) -> Violation:
    return Violation(check=check, severity=severity, message=message,
                     examples=[str(e) for e in list(examples)[:_MAX_EXAMPLES]])


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def read_fasta(path: Path):
    """Yield (seq_id, header, sequence, line_lengths) for each entry.

    `seq_id` is the first whitespace-delimited token after the '>', which is
    the same id `06_package_release.py` writes into the entry table.
    """
    name = None
    lines: list[str] = []

    def _emit():
        tok = name.split()
        return (tok[0][1:] if tok else "", name, "".join(lines),
                [len(x) for x in lines])

    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            raw = raw.rstrip("\n")
            if raw.startswith("#"):
                continue
            if raw.startswith(">"):
                if name is not None:
                    yield _emit()
                name, lines = raw, []
            elif name is not None:
                lines.append(raw.strip())
    if name is not None:
        yield _emit()


def header_fields(header: str) -> dict[str, str]:
    """The KEY=VALUE tokens of a header, as a dict."""
    out: dict[str, str] = {}
    for tok in header.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def _fasta_ids(path: Path) -> list[str]:
    return [i for i, _h, _s, _l in read_fasta(path)]


def _tsv_column(path: Path, column: str) -> list[str]:
    """Values of `column`, falling back to the first column if it is absent."""
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        try:
            head = next(r)
        except StopIteration:
            return []
        idx = head.index(column) if column in head else 0
        return [row[idx] for row in r if row]


# --------------------------------------------------------------------------
# release layout
# --------------------------------------------------------------------------

@dataclass
class Release:
    outdir: Path
    combined: Path | None = None
    decoy: Path | None = None
    per_type: list[Path] = field(default_factory=list)
    entries: Path | None = None
    manifest_txt: Path | None = None
    validation_md: Path | None = None
    problems: list[str] = field(default_factory=list)


def locate_release(outdir: Path) -> Release:
    """Find the pieces of a release without assuming the dataset name.

    Ambiguity is a problem to report, not a coin to flip: two files matching
    `*.target.fasta` means the caller pointed at two releases in one folder.
    """
    outdir = Path(outdir)
    rel = Release(outdir=outdir)
    if not outdir.is_dir():
        rel.problems.append(f"{outdir} is not a directory")
        return rel

    combined = sorted(p for p in outdir.glob("*.target.fasta") if p.is_file())
    if len(combined) == 1:
        rel.combined = combined[0]
    elif not combined:
        rel.problems.append(f"no *.target.fasta in {outdir}")
    else:
        rel.problems.append(
            f"{len(combined)} files match *.target.fasta in {outdir}: "
            + ", ".join(p.name for p in combined))

    decoy = sorted(p for p in outdir.glob("*.target_decoy.fasta")
                   if p.is_file())
    if len(decoy) == 1:
        rel.decoy = decoy[0]
    elif decoy:
        rel.problems.append(
            f"{len(decoy)} files match *.target_decoy.fasta in {outdir}")

    bc = outdir / "by_class"
    if bc.is_dir():
        rel.per_type = sorted(p for p in bc.glob("*.fasta") if p.is_file())

    entries = sorted(p for p in outdir.glob("*.entries.tsv") if p.is_file())
    if len(entries) == 1:
        rel.entries = entries[0]

    man = outdir / "MANIFEST.txt"
    if man.is_file():
        rel.manifest_txt = man

    val = outdir / "qc" / "uniprot_validation.md"
    if val.is_file():
        rel.validation_md = val
    return rel


# --------------------------------------------------------------------------
# line layout
# --------------------------------------------------------------------------

@dataclass
class Layout:
    """How a FASTA lays its sequences out across lines."""
    width: int | None = None      # interior line length; None if never wrapped
    max_line: int = 0
    ragged: bool = False
    detail: str = ""

    def describe(self) -> str:
        if self.ragged:
            return f"ragged ({self.detail})"
        if self.width is None:
            return "unwrapped"
        return f"wrapped at {self.width}"


def fasta_layout(path: Path) -> Layout:
    """Classify a FASTA as unwrapped, wrapped at a fixed width, or ragged.

    A file whose sequences all fit on one line is 'unwrapped'. It is
    compatible with any wrap width no shorter than its longest line, because
    nothing in it could reveal the width - see `layout_matches`.
    """
    lay = Layout()
    widths: set[int] = set()
    for sid, _h, _s, lens in read_fasta(path):
        if not lens:
            continue
        lay.max_line = max(lay.max_line, max(lens))
        if len(lens) == 1:
            continue
        interior = set(lens[:-1])
        if len(interior) > 1:
            lay.ragged = True
            lay.detail = f"{sid} has interior lines of {sorted(interior)}"
            return lay
        w = interior.pop()
        if lens[-1] > w:
            lay.ragged = True
            lay.detail = f"{sid} has a final line longer than {w}"
            return lay
        widths.add(w)
        if len(widths) > 1:
            lay.ragged = True
            lay.detail = f"wrap widths {sorted(widths)} in one file"
            return lay
    lay.width = widths.pop() if widths else None
    return lay


def layout_matches(out: Layout, tmpl: Layout) -> str | None:
    """None if `out` is consistent with `tmpl`, else why it is not."""
    if out.ragged:
        return f"output layout is {out.describe()}"
    if tmpl.ragged:
        return None                      # nothing coherent to match against
    if out.width == tmpl.width:
        return None
    if tmpl.width is None:               # template unwrapped, output wrapped
        return f"template is unwrapped but the output is {out.describe()}"
    if out.width is None:                # output unwrapped, template wrapped
        if out.max_line > tmpl.width:
            return (f"template is {tmpl.describe()} but the output is "
                    f"unwrapped (longest line {out.max_line})")
        return None                      # indeterminate: every line would fit
    return f"template is {tmpl.describe()} but the output is {out.describe()}"


# --------------------------------------------------------------------------
# the checks
# --------------------------------------------------------------------------

def check_i1_partition(rel: Release) -> list[Violation]:
    """I1 - the per-type files partition the combined file."""
    if rel.combined is None:
        return [_v("I1", WARNING, f"{_NOT_EVALUATED}: no combined FASTA found")]
    if not rel.per_type:
        return [_v("I1", WARNING,
                   f"{_NOT_EVALUATED}: no per-type files (by_class/*.fasta); "
                   f"the release was built without splitting")]

    combined = set(_fasta_ids(rel.combined))
    per = {p.name: set(_fasta_ids(p)) for p in rel.per_type}
    union: set[str] = set().union(*per.values()) if per else set()
    out: list[Violation] = []

    missing = combined - union
    if missing:
        out.append(_v("I1", ERROR,
                      f"{len(missing)} sequence(s) in {rel.combined.name} "
                      f"appear in no per-type file", sorted(missing)))
    extra = union - combined
    if extra:
        out.append(_v("I1", ERROR,
                      f"{len(extra)} sequence(s) in the per-type files are "
                      f"absent from {rel.combined.name}", sorted(extra)))

    names = sorted(per)
    overlaps = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            both = per[a] & per[b]
            if both:
                overlaps.append(f"{a} & {b}: {len(both)} shared "
                                f"(e.g. {sorted(both)[0]})")
    if overlaps:
        out.append(_v("I1", ERROR,
                      f"{len(overlaps)} per-type file pair(s) overlap; the "
                      f"split must be a partition", overlaps))
    return out


def check_i2_disposition(manifest: Path | None,
                         disposition: Path | None) -> list[Violation]:
    """I2 - every input variant appears exactly once in the disposition."""
    if manifest is None or not Path(manifest).is_file():
        return [_v("I2", WARNING,
                   f"{_NOT_EVALUATED}: no unified variant manifest supplied")]
    if disposition is None or not Path(disposition).is_file():
        return [_v("I2", WARNING,
                   f"{_NOT_EVALUATED}: no disposition table supplied")]

    want = _tsv_column(Path(manifest), "variant_id")
    got = _tsv_column(Path(disposition), "variant_id")
    out: list[Violation] = []

    seen: dict[str, int] = {}
    for vid in got:
        seen[vid] = seen.get(vid, 0) + 1
    dupes = sorted(k for k, n in seen.items() if n > 1)
    if dupes:
        out.append(_v("I2", ERROR,
                      f"{len(dupes)} variant(s) appear more than once in "
                      f"{Path(disposition).name}", dupes))

    want_set, got_set = set(want), set(seen)
    absent = sorted(want_set - got_set)
    if absent:
        out.append(_v("I2", ERROR,
                      f"{len(absent)} input variant(s) have no disposition; "
                      f"every input must have a recorded outcome", absent))
    unknown = sorted(got_set - want_set)
    if unknown:
        out.append(_v("I2", ERROR,
                      f"{len(unknown)} disposition row(s) name a variant that "
                      f"is not in the input manifest", unknown))
    return out


def _classify(header: str, seq_id: str) -> str:
    """Which of the three buckets a combined-file entry belongs to."""
    vt = header_fields(header).get("VT")
    if vt == "REFERENCE":
        return "wildtype"
    if vt:
        return "variant"
    if seq_id.startswith(("sp|", "tr|")):
        # Appended reference-proteome entries keep their pristine UniProt
        # header, so they carry no VT. That is the only legitimate way for
        # an entry in the combined file to have none.
        return "reference_proteome"
    return "unclassified"


def check_i3_counts(rel: Release) -> list[Violation]:
    """I3 - variant + wild-type + reference counts sum to the combined file."""
    if rel.combined is None:
        return [_v("I3", WARNING, f"{_NOT_EVALUATED}: no combined FASTA found")]

    buckets: dict[str, list[str]] = {"variant": [], "wildtype": [],
                                     "reference_proteome": [],
                                     "unclassified": []}
    for sid, hdr, _s, _l in read_fasta(rel.combined):
        buckets[_classify(hdr, sid)].append(sid)
    total = sum(len(v) for v in buckets.values())
    out: list[Violation] = []

    if buckets["unclassified"]:
        out.append(_v("I3", ERROR,
                      f"{len(buckets['unclassified'])} entries in "
                      f"{rel.combined.name} carry neither a VT= field nor an "
                      f"sp|/tr| reference id, so they count towards no "
                      f"category", buckets["unclassified"]))

    n_var = len(buckets["variant"])
    n_wt = len(buckets["wildtype"])
    n_ref = len(buckets["reference_proteome"])

    # The entry table covers variant and wild-type entries; the appended
    # reference proteome is not tabulated. Double-counted wild-types show up
    # here as an entry table longer than the two buckets it describes.
    if rel.entries is not None:
        n_rows = len(_tsv_column(rel.entries, "seq_id"))
        if n_rows != n_var + n_wt:
            out.append(_v("I3", ERROR,
                          f"{rel.entries.name} has {n_rows} rows but the "
                          f"combined file holds {n_var} variant + {n_wt} "
                          f"wild-type = {n_var + n_wt} tabulated entries"))
        elif n_rows + n_ref != total:
            out.append(_v("I3", ERROR,
                          f"{n_rows} tabulated + {n_ref} appended reference "
                          f"proteins = {n_rows + n_ref}, but "
                          f"{rel.combined.name} holds {total} sequences"))

    # The per-type split must agree about how many of each there are.
    by_name = {p.name: p for p in rel.per_type}
    for suffix, expected, label in ((".REFERENCE.fasta", n_wt, "wild-type"),
                                    (".REFERENCE_PROTEOME.fasta", n_ref,
                                     "appended reference")):
        hit = [p for n, p in by_name.items() if n.endswith(suffix)]
        if len(hit) == 1:
            got = len(_fasta_ids(hit[0]))
            if got != expected:
                out.append(_v("I3", ERROR,
                              f"{hit[0].name} holds {got} sequences but the "
                              f"combined file holds {expected} {label} "
                              f"entries"))
    return out


def check_i4_unique_ids(rel: Release) -> list[Violation]:
    """I4 - every sequence id is unique, in every FASTA of the release."""
    targets = [p for p in (rel.combined, rel.decoy) if p is not None]
    targets += rel.per_type
    if not targets:
        return [_v("I4", WARNING, f"{_NOT_EVALUATED}: no FASTA files found")]

    out: list[Violation] = []
    for path in targets:
        seen: dict[str, int] = {}
        for sid in _fasta_ids(path):
            seen[sid] = seen.get(sid, 0) + 1
        dupes = sorted(k for k, n in seen.items() if n > 1)
        if dupes:
            out.append(_v("I4", ERROR,
                          f"{path.name} has {len(dupes)} duplicated sequence "
                          f"id(s); a search engine cannot tell them apart",
                          dupes))
    return out


def check_i5_layout(rel: Release, proteome: Path | None) -> list[Violation]:
    """I5 - sequence line layout matches the supplied proteome template."""
    if rel.combined is None:
        return [_v("I5", WARNING, f"{_NOT_EVALUATED}: no combined FASTA found")]

    tmpl_path = None
    if proteome is not None and Path(proteome).is_file():
        tmpl_path = Path(proteome)
    else:
        hit = [p for p in rel.per_type
               if p.name.endswith(".REFERENCE_PROTEOME.fasta")]
        if len(hit) == 1:
            tmpl_path = hit[0]
    if tmpl_path is None:
        return [_v("I5", WARNING,
                   f"{_NOT_EVALUATED}: no reference proteome to use as the "
                   f"layout template")]

    tmpl = fasta_layout(tmpl_path)
    out: list[Violation] = []
    checked = [rel.combined] + ([rel.decoy] if rel.decoy else []) + rel.per_type
    for path in checked:
        if path == tmpl_path:
            continue
        why = layout_matches(fasta_layout(path), tmpl)
        if why:
            out.append(_v("I5", ERROR,
                          f"{path.name} does not match the layout of "
                          f"{tmpl_path.name}: {why}"))
    return out


_B1 = re.compile(r"^##\s+B1\..*$", re.MULTILINE)
_AGREE = re.compile(r"Exact-agreement:\s*\*\*([0-9.]+)%\*\*")


def check_i6_agreement(rel: Release, min_agreement: float) -> list[Violation]:
    """I6 - reference translations agree with the annotation source's own."""
    if rel.validation_md is None:
        return [_v("I6", WARNING,
                   f"{_NOT_EVALUATED}: qc/uniprot_validation.md is not in the "
                   f"release, so the translation check did not run")]

    text = rel.validation_md.read_text(encoding="utf-8")
    m = _B1.search(text)
    if not m:
        return [_v("I6", WARNING,
                   f"{_NOT_EVALUATED}: the validation report has no GENCODE "
                   f"translation section; the annotation source's own "
                   f"translations were not supplied to this run")]

    section = text[m.end():]
    nxt = re.search(r"^##\s", section, re.MULTILINE)
    if nxt:
        section = section[:nxt.start()]

    hit = _AGREE.search(section)
    if not hit:
        # The section exists, so the check was attempted. Failing to read its
        # result must not be indistinguishable from passing it.
        return [_v("I6", ERROR,
                   f"the translation check ran but its agreement rate could "
                   f"not be read from {rel.validation_md.name}; the report "
                   f"format has changed and this check is now blind")]

    rate = float(hit.group(1)) / 100.0
    if rate < min_agreement:
        return [_v("I6", ERROR,
                   f"reference translations agree with the annotation source "
                   f"at {rate:.1%}, below the {min_agreement:.0%} threshold")]
    return []


# Two fusion breakpoints this close together on one chromosome are far more
# likely to be a read-through or an intra-cluster rearrangement than a true
# fusion of two distant genes. Clustered loci - the protocadherin gamma
# family, the immunoglobulin and HLA regions - splice many variable exons
# onto a shared constant region, so a caller can report a "fusion" whose
# protein is simply an existing family member. 1 Mb comfortably spans the
# largest of those clusters.
_INTRALOCUS_WINDOW = 1_000_000


def _fusion_breakpoints(loc: str) -> tuple[str, int, str, int] | None:
    """Parse a fusion `LOC=<chrom>:<pos>::<chrom>:<pos>` into its two ends.

    Returns None for any other LOC shape, including the small-variant form
    `chr1:1049980G>C`, which has no second breakpoint.
    """
    if "::" not in loc:
        return None
    left, _, right = loc.partition("::")
    ends = []
    for half in (left, right):
        chrom, _, pos = half.rpartition(":")
        if not chrom or not pos.isdigit():
            return None
        ends.append((chrom, int(pos)))
    return ends[0][0], ends[0][1], ends[1][0], ends[1][1]


def _intralocus_note(loc: str) -> str | None:
    """A note naming the likely cause when both breakpoints share a locus."""
    ends = _fusion_breakpoints(loc)
    if ends is None:
        return None
    c1, p1, c2, p2 = ends
    if c1 != c2:
        return None
    gap = abs(p2 - p1)
    if gap > _INTRALOCUS_WINDOW:
        return None
    return (f"both breakpoints on {c1}, {gap / 1000:.0f} kb apart - likely a "
            f"read-through or intra-cluster event rather than a fusion of two "
            f"distant genes")


def check_i7_unchanged(rel: Release) -> list[Violation]:
    """I7 - no sequence-changing entry has the sequence of its own reference."""
    if rel.combined is None:
        return [_v("I7", WARNING, f"{_NOT_EVALUATED}: no combined FASTA found")]

    ref_by_tx: dict[str, str] = {}
    candidates: list[tuple[str, dict[str, str], str]] = []
    for sid, hdr, seq, _l in read_fasta(rel.combined):
        f = header_fields(hdr)
        if f.get("VT") == "REFERENCE":
            tx = f.get("TX")
            if tx:
                ref_by_tx[tx] = seq
        elif f.get("CSQ") in SEQUENCE_CHANGING:
            candidates.append((sid, f, seq))

    bad: list[str] = []
    n_intralocus = 0
    for sid, f, seq in candidates:
        # Cross-class deduplication merges a variant sequence with an
        # identical reference and keeps the variant label, recording the
        # other in ALSO=. Either signal means the protein did not change.
        also = f.get("ALSO", "").split(",") if f.get("ALSO") else []
        tx = f.get("TX", "")
        if "REFERENCE" not in also and not (tx and ref_by_tx.get(tx) == seq):
            continue
        note = _intralocus_note(f.get("LOC", ""))
        if note:
            n_intralocus += 1
            bad.append(f"{sid} (CSQ={f.get('CSQ', '')}; {note})")
        else:
            bad.append(f"{sid} (CSQ={f.get('CSQ', '')})")

    if not bad:
        return []

    msg = (f"{len(bad)} entry(ies) labelled with a sequence-changing "
           f"consequence have a protein identical to a reference protein; "
           f"either the label or the entry is wrong")
    if n_intralocus:
        # Naming the likely cause is the difference between a warning
        # someone acts on and a warning someone learns to ignore.
        msg += (f". {n_intralocus} of them have both fusion breakpoints "
                f"within {_INTRALOCUS_WINDOW // 1000} kb on one chromosome, "
                f"which points at the fusion caller rather than at this "
                f"pipeline - check the call's supporting evidence before "
                f"treating it as a translation fault")
    return [_v("I7", WARNING, msg, bad)]


_MAN_LINE = re.compile(r"^([0-9a-fA-F]{64})\s\s(.+?)\s\s\([\d,]+ bytes\)\s*$")


def check_i8_manifest(rel: Release) -> list[Violation]:
    """I8 - MANIFEST.txt checksums match the files on disk."""
    if rel.manifest_txt is None:
        return [_v("I8", WARNING,
                   f"{_NOT_EVALUATED}: the release has no MANIFEST.txt")]

    listed: dict[str, str] = {}
    unreadable = []
    for line in rel.manifest_txt.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        m = _MAN_LINE.match(line)
        if m:
            listed[m.group(2).replace("\\", "/")] = m.group(1).lower()
        else:
            unreadable.append(line[:80])

    out: list[Violation] = []
    if unreadable:
        out.append(_v("I8", ERROR,
                      f"{len(unreadable)} line(s) of {rel.manifest_txt.name} "
                      f"are not '<sha256>  <path>  (<n> bytes)'", unreadable))

    missing, changed = [], []
    for rel_path, want in sorted(listed.items()):
        f = rel.outdir / rel_path
        if not f.is_file():
            missing.append(rel_path)
        elif sha256(f).lower() != want:
            changed.append(rel_path)
    if missing:
        out.append(_v("I8", ERROR,
                      f"{len(missing)} file(s) in the manifest are not on "
                      f"disk", missing))
    if changed:
        out.append(_v("I8", ERROR,
                      f"{len(changed)} file(s) differ from their manifest "
                      f"checksum", changed))

    on_disk = {str(p.relative_to(rel.outdir)).replace("\\", "/")
               for p in rel.outdir.rglob("*")
               if p.is_file() and p.name != rel.manifest_txt.name}
    unlisted = sorted(on_disk - set(listed))
    if unlisted:
        out.append(_v("I8", ERROR,
                      f"{len(unlisted)} file(s) in the release are absent "
                      f"from the manifest; it is stale", unlisted))
    return out


def check_i9_decoys(rel: Release,
                    decoys_requested: bool | None) -> list[Violation]:
    """I9 - decoy count equals target count when decoys were requested."""
    if decoys_requested is None:
        decoys_requested = rel.decoy is not None
    if not decoys_requested:
        return []
    if rel.decoy is None:
        return [_v("I9", ERROR,
                   "decoys were requested but no *.target_decoy.fasta was "
                   "written")]
    if rel.combined is None:
        return [_v("I9", WARNING,
                   f"{_NOT_EVALUATED}: no combined FASTA to compare against")]

    n_target = len(_fasta_ids(rel.combined))
    ids = _fasta_ids(rel.decoy)
    n_decoy = sum(1 for i in ids if i.startswith("DECOY_"))
    n_other = len(ids) - n_decoy
    out: list[Violation] = []
    if n_decoy != n_target:
        out.append(_v("I9", ERROR,
                      f"{rel.decoy.name} holds {n_decoy} decoys for "
                      f"{n_target} targets; FDR estimation assumes one decoy "
                      f"per target"))
    if n_other != n_target:
        out.append(_v("I9", ERROR,
                      f"{rel.decoy.name} holds {n_other} target entries but "
                      f"{rel.combined.name} holds {n_target}"))
    return out


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------

def check_release(outdir: Path, manifest: Path | None = None,
                  disposition: Path | None = None, *,
                  proteome: Path | None = None,
                  min_agreement: float = 0.90,
                  decoys_requested: bool | None = None,
                  skip_checks: Iterable[str] = ()) -> list[Violation]:
    """Run every invariant against a finished release.

    `outdir` is the release directory. `manifest` is the unified variant
    manifest the run started from and `disposition` the per-variant outcome
    table; both are needed for I2 and either may be None.

    The keyword-only arguments refine checks that would otherwise have to
    guess: `proteome` is the template for I5, `min_agreement` the threshold
    for I6, and `decoys_requested` says whether decoys were asked for, so a
    missing decoy file is an error rather than an untested case.

    `skip_checks` names checks not to run. A skipped check still returns a
    warning saying it was skipped, so a release built with one turned off
    cannot be mistaken for one that passed it. Unknown ids raise, because a
    typo that quietly disables nothing - or quietly disables the wrong
    thing - is the failure this whole module exists to prevent.

    Returns every violation found. An empty list means the release is
    internally consistent; it does not mean the biology is right.
    """
    skip = {s.strip().upper() for s in skip_checks if s and s.strip()}
    unknown = sorted(skip - set(CHECK_TITLES))
    if unknown:
        raise ValueError(
            f"unknown invariant id(s) {', '.join(unknown)}; "
            f"choose from {', '.join(sorted(CHECK_TITLES))}")

    rel = locate_release(Path(outdir))
    out: list[Violation] = [_v("release", ERROR, p) for p in rel.problems]
    runners = {
        "I1": lambda: check_i1_partition(rel),
        "I2": lambda: check_i2_disposition(manifest, disposition),
        "I3": lambda: check_i3_counts(rel),
        "I4": lambda: check_i4_unique_ids(rel),
        "I5": lambda: check_i5_layout(rel, proteome),
        "I6": lambda: check_i6_agreement(rel, min_agreement),
        "I7": lambda: check_i7_unchanged(rel),
        "I8": lambda: check_i8_manifest(rel),
        "I9": lambda: check_i9_decoys(rel, decoys_requested),
    }
    for cid in sorted(runners):
        if cid in skip:
            out.append(_v(cid, WARNING,
                          f"{_SKIPPED}: turned off for this run with "
                          f"--skip-invariant {cid}"))
        else:
            out += runners[cid]()
    return out


def format_report(violations: list[Violation]) -> str:
    """A human-readable summary, listing every check so silence is visible."""
    by_check: dict[str, list[Violation]] = {}
    for v in violations:
        by_check.setdefault(v.check, []).append(v)

    lines = ["release invariants", "-" * 72]
    for cid in sorted(CHECK_TITLES):
        hits = by_check.get(cid, [])
        if not hits:
            state = "ok"
        elif any(h.severity == ERROR for h in hits):
            state = "FAIL"
        elif all(h.message.startswith(_SKIPPED) for h in hits):
            state = "OFF"
        elif all(h.message.startswith(_NOT_EVALUATED) for h in hits):
            state = "skip"
        else:
            state = "warn"
        lines.append(f"  {cid}  {state:<4}  {CHECK_TITLES[cid]}")
        for h in hits:
            lines.append(f"        {h.severity}: {h.message}")
            for e in h.examples:
                lines.append(f"          - {e}")
    for v in by_check.get("release", []):
        lines.append(f"  release  FAIL  {v.message}")

    n_err = sum(1 for v in violations if v.severity == ERROR)
    n_warn = sum(1 for v in violations
                 if v.severity == WARNING
                 and not v.message.startswith((_NOT_EVALUATED, _SKIPPED)))
    n_off = sum(1 for v in violations if v.message.startswith(_SKIPPED))
    lines.append("-" * 72)
    tail = f"  {n_err} error(s), {n_warn} warning(s)"
    if n_off:
        # A run with checks turned off must say so on the summary line, not
        # only in the body where it can be scrolled past.
        tail += f", {n_off} check(s) turned off"
    lines.append(tail)
    return "\n".join(lines)
