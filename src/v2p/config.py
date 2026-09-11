"""Run configuration: the YAML schema documented in docs/USAGE.md.

One config file plus the input checksums must fully determine a run, so
this module is the single place that knows

  * what the legal keys are,
  * what type each one holds and what values it may take,
  * what its default is,
  * which ``v2p run`` command-line option it corresponds to.

Everything else - the CLI, ``--write-config``, and the provenance record -
is generated from :data:`SCHEMA`. A key added there is accepted by the
loader, emitted by ``--write-config`` and recorded in the provenance JSON
without any further wiring, which is the only way those three stay in
agreement.

Paths are used exactly as written, interpreted relative to the working
directory the command runs in. They are deliberately *not* resolved
against the config file's own directory: silently rewriting a path is the
kind of helpfulness that makes a run hard to explain afterwards.

YAML parsing prefers PyYAML. When PyYAML is absent a small parser covers
the documented subset - two levels of mapping, scalars, inline lists - so
that ``v2p run --config`` works in a standard-library-only environment,
which the rest of the parsing and QC stages already promise.
"""

from __future__ import annotations

import difflib
import re
from dataclasses import dataclass, field
from typing import Any

__all__ = [
    "ConfigError", "Key", "SCHEMA", "RunConfig",
    "default_config", "load_config", "loads_config", "config_from_namespace",
]


class ConfigError(ValueError):
    """A config file that cannot be trusted to describe a run.

    Raised for unknown sections and keys, wrong types, values outside the
    allowed set, and malformed YAML. Never raised for a *missing* key: an
    omitted key takes its documented default.
    """


@dataclass(frozen=True)
class Key:
    """One configuration key.

    ``dest`` is the ``argparse`` destination on ``v2p run`` that this key
    sets, which is what lets the config file and the command line be
    merged by one rule instead of by hand per option.
    """

    name: str
    dest: str
    type: str                 # str | path | bool | int | float | list
    default: Any = None
    choices: tuple[str, ...] | None = None
    comment: str = ""


# The order here is the order --write-config emits, so it is also the
# order a reader meets the options in. Keep it stable.
SCHEMA: dict[str, tuple[Key, ...]] = {
    "inputs": (
        Key("folder", "folder", "path", None,
            comment="directory of variant calls; contents detected by type"),
        Key("small_variants", "small_variants", "path", None,
            comment="explicit path, overriding detection"),
        Key("rna_editing", "rna_editing", "path", None),
        Key("fusion", "fusion_calls", "path", None),
        Key("splicing", "splicing", "path", None),
    ),
    "reference": (
        Key("dir", "ref", "path", None,
            comment="directory holding genome, annotation and proteome"),
        Key("genome", "genome", "path", None,
            comment="explicit path, overriding what --ref finds"),
        Key("annotation", "annotation", "path", None),
        Key("proteome", "proteome", "path", None),
        Key("translations", "translations", "path", None,
            comment="annotation source's own protein translations"),
    ),
    "translate": (
        Key("transcript_mode", "transcript_mode", "str", "all",
            ("all", "representative")),
        Key("keep_unchanged", "keep_unchanged", "bool", True,
            comment="keep synonymous and UTR variants"),
        # Only 'auto' exists: table 2 on the mitochondrial contig, table 1
        # elsewhere. Listing a value the pipeline cannot honour would put
        # a setting in the provenance record that did not happen, so the
        # loader rejects anything else rather than accepting and ignoring
        # it.
        Key("genetic_code", "genetic_code", "str", "auto", ("auto",),
            comment="auto = table 2 on the mitochondrial contig, 1 elsewhere"),
        Key("species", "species", "str", "human",
            comment="species name, or a path to a config/species/*.yaml"),
        # Off by default and recorded either way. A run that turned
        # three-frame translation on must say so in its own provenance,
        # or the database size cannot be explained later.
        Key("include_noncanonical", "include_noncanonical", "bool", False,
            comment="three-frame translate non-coding transcripts"),
        Key("noncanonical_min_aa", "nc_min_aa", "int", 30,
            comment="minimum non-canonical ORF length in residues"),
        Key("noncanonical_any_start", "nc_any_start", "bool", False,
            comment="keep ORFs that do not begin at ATG"),
        # Combining changes what is in the database, so a config that
        # omitted these would reproduce a different release from the run
        # that wrote it.
        Key("combine_variants", "combine_variants", "bool", False,
            comment="emit proteins carrying all co-occurring variants"),
        Key("combine_max", "combine_max", "int", 8,
            comment="most variants combined on one transcript"),
        Key("allow_unphased", "allow_unphased", "bool", True,
            comment="combine variants the caller did not phase"),
        Key("min_af", "min_af", "float", None,
            comment="drop variants below this INFO allele frequency"),
        Key("max_combinatorial_fraction", "max_combinatorial_fraction",
            "float", None,
            comment="fail if combinatorial entries exceed this share"),
    ),
    "output": (
        Key("dir", "outdir", "path", "v2p_output"),
        Key("name", "name", "str", "variant_proteome"),
        Key("header_style", "header_style", "str", "uniprot",
            ("uniprot", "peff", "pvac", "descriptive")),
        Key("decoys", "decoys", "str", "pseudo_reverse",
            ("none", "reverse", "pseudo_reverse", "shuffle")),
        Key("split_by_type", "split_by_type", "bool", True),
        Key("append_reference", "append_reference", "bool", True),
        Key("logdir", "logdir", "path", "logs",
            comment="keep outside output.dir or the release manifest goes stale"),
    ),
    "validate": (
        Key("min_translation_agreement", "min_agreement", "float", 0.90),
        Key("skip_invariant", "skip_invariant", "list", None,
            comment="release invariants to turn off, e.g. [I5]"),
        Key("force", "force", "bool", False,
            comment="proceed despite detection problems or a failed gate"),
    ),
}

_BY_SECTION_KEY: dict[tuple[str, str], Key] = {
    (sec, k.name): k for sec, keys in SCHEMA.items() for k in keys
}
_BY_DEST: dict[str, tuple[str, Key]] = {
    k.dest: (sec, k) for sec, keys in SCHEMA.items() for k in keys
}

# `inputs: auto` is the documented way of saying "detect everything".
AUTO = "auto"


def default_config() -> dict[str, dict[str, Any]]:
    """The full config with every key at its documented default."""
    return {sec: {k.name: k.default for k in keys}
            for sec, keys in SCHEMA.items()}


# ---------------------------------------------------------------- typing

def _did_you_mean(name: str, options) -> str:
    near = difflib.get_close_matches(name, list(options), n=1, cutoff=0.6)
    return f"; did you mean '{near[0]}'?" if near else ""


def _coerce(key: Key, raw: Any, where: str) -> Any:
    """Validate and convert one value, or raise ConfigError naming it."""
    if raw is None:
        return None

    if key.type == "bool":
        if isinstance(raw, bool):
            return raw
        if isinstance(raw, str) and raw.strip().lower() in (
                "true", "false", "yes", "no", "on", "off"):
            return raw.strip().lower() in ("true", "yes", "on")
        raise ConfigError(
            f"{where}: expected true or false, got {raw!r}")

    if key.type == "int":
        if isinstance(raw, bool) or isinstance(raw, (list, dict)):
            raise ConfigError(f"{where}: expected a whole number, got {raw!r}")
        try:
            return int(str(raw).strip())
        except (TypeError, ValueError):
            raise ConfigError(
                f"{where}: expected a whole number, got {raw!r}") from None

    if key.type == "float":
        if isinstance(raw, bool) or isinstance(raw, (list, dict)):
            raise ConfigError(f"{where}: expected a number, got {raw!r}")
        try:
            val = float(raw)
        except (TypeError, ValueError):
            raise ConfigError(
                f"{where}: expected a number, got {raw!r}") from None
        if key.name == "min_translation_agreement" and not 0.0 <= val <= 1.0:
            raise ConfigError(
                f"{where}: must be a fraction between 0 and 1, got {val!r}")
        return val

    if key.type == "list":
        if isinstance(raw, (list, tuple)):
            items = [str(x).strip() for x in raw]
        elif isinstance(raw, str):
            items = [x.strip() for x in raw.split(",")]
        else:
            raise ConfigError(f"{where}: expected a list, got {raw!r}")
        items = [x for x in items if x]
        return items or None

    # str / path
    if isinstance(raw, (bool, list, dict, tuple)):
        raise ConfigError(f"{where}: expected a string, got {raw!r}")
    val = str(raw)
    if key.choices and val not in key.choices:
        raise ConfigError(
            f"{where}: '{val}' is not one of {', '.join(key.choices)}"
            + _did_you_mean(val, key.choices))
    return val


# ---------------------------------------------------------------- parsing

_COMMENT = re.compile(r"(?:^|\s)#.*$")


def _mini_scalar(text: str) -> Any:
    t = text.strip()
    if t == "" or t in ("null", "~", "None"):
        return None
    if len(t) >= 2 and t[0] == t[-1] and t[0] in "'\"":
        return t[1:-1]
    if t.startswith("[") and t.endswith("]"):
        inner = t[1:-1].strip()
        if not inner:
            return []
        return [_mini_scalar(x) for x in inner.split(",")]
    low = t.lower()
    if low in ("true", "yes", "on"):
        return True
    if low in ("false", "no", "off"):
        return False
    try:
        return int(t)
    except ValueError:
        pass
    try:
        return float(t)
    except ValueError:
        pass
    return t


def _mini_yaml(text: str, where: str) -> dict:
    """Parse the documented subset without PyYAML.

    Two levels only: top-level ``section:`` headers with indented
    ``key: value`` pairs underneath, plus top-level ``key: value`` for the
    ``inputs: auto`` short form. Anything else is rejected rather than
    guessed at, because a config file that parses to something other than
    what it says is worse than one that will not parse at all.
    """
    doc: dict[str, Any] = {}
    section: str | None = None
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = _COMMENT.sub("", raw).rstrip()
        if not line.strip():
            continue
        if line.lstrip().startswith("- "):
            raise ConfigError(
                f"{where}: line {lineno}: block lists are not supported "
                f"without PyYAML; write it inline as [a, b]")
        indent = len(line) - len(line.lstrip())
        body = line.strip()
        if ":" not in body:
            raise ConfigError(
                f"{where}: line {lineno}: expected 'key: value', got {body!r}")
        name, _, value = body.partition(":")
        name, value = name.strip(), value.strip()
        if indent == 0:
            if value == "":
                section = name
                doc[name] = {}
            else:
                doc[name] = _mini_scalar(value)
                section = None
        else:
            if section is None:
                raise ConfigError(
                    f"{where}: line {lineno}: indented key '{name}' is not "
                    f"under any section")
            doc[section][name] = _mini_scalar(value)
    return doc


def _parse(text: str, where: str) -> dict:
    try:
        import yaml                                   # noqa: PLC0415
    except ImportError:
        return _mini_yaml(text, where)
    try:
        doc = yaml.safe_load(text)
    except Exception as exc:                          # yaml.YAMLError et al
        raise ConfigError(f"{where}: not valid YAML: {exc}") from exc
    if doc is None:
        return {}
    if not isinstance(doc, dict):
        raise ConfigError(
            f"{where}: expected a mapping of sections, got {type(doc).__name__}")
    return doc


# ---------------------------------------------------------------- the type

@dataclass
class RunConfig:
    """A validated configuration, plus which keys the file actually set.

    ``provided`` matters for precedence: a key the file did not mention
    must not override a command-line option, and the only way to tell the
    difference between "absent" and "set to the default" is to record it.
    """

    values: dict[str, dict[str, Any]] = field(default_factory=default_config)
    provided: set[str] = field(default_factory=set)
    source: str | None = None

    def get(self, section: str, name: str) -> Any:
        return self.values[section][name]

    def to_dict(self) -> dict[str, dict[str, Any]]:
        """Every key with its resolved value - what goes in provenance."""
        return {sec: dict(vals) for sec, vals in self.values.items()}

    def overrides(self) -> dict[str, Any]:
        """``argparse`` dest -> value, for the keys the file actually set."""
        out: dict[str, Any] = {}
        for sec, keys in SCHEMA.items():
            for k in keys:
                if f"{sec}.{k.name}" in self.provided:
                    out[k.dest] = self.values[sec][k.name]
        return out

    def to_yaml(self) -> str:
        """Emit the config as YAML. Deterministic: same values, same bytes."""
        lines = [
            "# v2p run configuration.",
            "# Generated by `v2p run --write-config`. Every key below is",
            "# recorded in the run's provenance JSON.",
            "# Paths are relative to the directory the command runs in.",
        ]
        for sec, keys in SCHEMA.items():
            lines.append("")
            lines.append(f"{sec}:")
            for k in keys:
                line = f"  {k.name}: {_render(self.values[sec][k.name])}"
                line = line.rstrip()
                if k.comment:
                    line += f"    # {k.comment}"
                lines.append(line)
        return "\n".join(lines) + "\n"


def _render(val: Any) -> str:
    if val is None:
        return ""
    if isinstance(val, bool):
        return "true" if val else "false"
    if isinstance(val, float):
        return repr(val)
    if isinstance(val, (list, tuple)):
        return "[" + ", ".join(str(x) for x in val) + "]"
    text = str(val)
    if text == "" or text != text.strip() or text[0] in "'\"[{#":
        return f"'{text}'"
    return text


# ---------------------------------------------------------------- loading

def loads_config(text: str, where: str = "config") -> RunConfig:
    """Parse and validate config text. See :func:`load_config`."""
    doc = _parse(text, where)
    cfg = RunConfig(values=default_config(), source=where)

    for sec_name, body in doc.items():
        if sec_name not in SCHEMA:
            raise ConfigError(
                f"{where}: unknown section '{sec_name}'; valid sections are "
                + ", ".join(SCHEMA) + _did_you_mean(sec_name, SCHEMA))
        # `inputs: auto` is the documented short form for "detect it all".
        if isinstance(body, str) and body.strip().lower() == AUTO:
            continue
        if body is None:
            continue
        if not isinstance(body, dict):
            raise ConfigError(
                f"{where}: section '{sec_name}' must be a mapping of keys, "
                f"got {type(body).__name__}")
        for key_name, raw in body.items():
            key = _BY_SECTION_KEY.get((sec_name, str(key_name)))
            if key is None:
                valid = [k.name for k in SCHEMA[sec_name]]
                raise ConfigError(
                    f"{where}: unknown key '{sec_name}.{key_name}'; "
                    f"'{sec_name}' accepts " + ", ".join(valid)
                    + _did_you_mean(str(key_name), valid))
            cfg.values[sec_name][key.name] = _coerce(
                key, raw, f"{where}: {sec_name}.{key.name}")
            cfg.provided.add(f"{sec_name}.{key.name}")
    return cfg


def load_config(path) -> RunConfig:
    """Read and validate a config file.

    Raises :class:`ConfigError` with the file name and the offending key
    for anything wrong with it, and :class:`FileNotFoundError` if it is
    not there.
    """
    from pathlib import Path                          # noqa: PLC0415
    p = Path(path)
    try:
        text = p.read_text(encoding="utf-8")
    except FileNotFoundError:
        raise ConfigError(f"no such config file: {p}") from None
    return loads_config(text, where=str(p))


def config_from_namespace(ns) -> RunConfig:
    """The config an already-parsed ``v2p run`` command line implies.

    Used both by ``--write-config`` and by the provenance record, so what
    a run reports and what ``--write-config`` emits cannot drift apart.
    """
    cfg = RunConfig(values=default_config(), source="command line")
    for dest, (sec, key) in _BY_DEST.items():
        if not hasattr(ns, dest):
            continue
        val = getattr(ns, dest)
        if val is None:
            continue
        if key.type == "list":
            val = sorted({str(x) for x in val}) or None
        elif key.type == "path" or key.type == "str":
            val = str(val)
        cfg.values[sec][key.name] = val
        cfg.provided.add(f"{sec}.{key.name}")
    return cfg
