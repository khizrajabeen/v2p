"""Species-specific facts, so the pipeline is not human-only.

What actually varies between species is small and dull: the name and taxon
id that go in a header, the suffix in a UniProt entry name (`TP53_HUMAN`
against `Trp53_MOUSE`), which contig carries the mitochondrial genome, and
the contig lengths the reference audit asserts against. None of it belongs
hard-coded in a header formatter.

The default is human, and every function here defaults to it, so code that
does not care about species keeps working and produces byte-identical
output. That is deliberate: M4 is a refactor, and a refactor that changes
output has failed.

Contig lengths are optional. A species with none cannot have its assembly
asserted, which is a weaker guarantee rather than an error - the audit
says so instead of inventing an expectation it cannot check.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "Species", "HUMAN", "MOUSE", "RAT", "ZEBRAFISH", "BUILTIN",
    "SUPPORTED_MITO_CODES",
    "from_mapping", "load_species", "species_dir", "available",
]


# The genetic codes seqops.py actually implements. Declaring any other in
# a species file is refused rather than accepted and ignored: an
# invertebrate (table 5) or yeast (table 3) mitochondrion would otherwise
# be translated with the vertebrate table and produce plausible, wrong
# protein - the exact failure mode this project keeps meeting.
SUPPORTED_MITO_CODES = {
    1: "standard",
    2: "vertebrate mitochondrial",
}


@dataclass(frozen=True)
class Species:
    """Everything the pipeline needs to know about an organism."""

    scientific_name: str = "Homo sapiens"
    taxon_id: int = 9606
    # Suffix in a UniProt entry name: the HUMAN in `TP53_HUMAN`.
    entry_suffix: str = "HUMAN"
    common_name: str = "human"
    assembly: str = "GRCh38"
    # Contigs translated with the mitochondrial code. Naming differs
    # between references, so several spellings are accepted.
    mito_contigs: tuple[str, ...] = ("chrM", "MT", "M", "chrMT")
    # NCBI translation table for those contigs. 2 is the vertebrate
    # mitochondrial code; invertebrates and yeast differ.
    mito_code: int = 2
    # Optional. Empty means the audit cannot assert the assembly, which it
    # reports rather than treating as a pass.
    contig_lengths: dict[str, int] = field(default_factory=dict)

    def is_mito(self, contig: str) -> bool:
        return contig in self.mito_contigs

    def entry_name(self, gene: str) -> str:
        """`TP53` -> `TP53_HUMAN`. Gene case is left as the source had it."""
        return f"{gene}_{self.entry_suffix}"

    def os_ox(self) -> str:
        """The `OS=... OX=...` fragment every UniProt-style header carries."""
        return f"OS={self.scientific_name} OX={self.taxon_id}"


HUMAN = Species(
    scientific_name="Homo sapiens", taxon_id=9606, entry_suffix="HUMAN",
    common_name="human", assembly="GRCh38",
    contig_lengths={"chr1": 248956422, "chr2": 242193529,
                    "chr17": 83257441},
)

MOUSE = Species(
    scientific_name="Mus musculus", taxon_id=10090, entry_suffix="MOUSE",
    common_name="mouse", assembly="GRCm39",
    contig_lengths={"chr1": 195154279, "chr2": 181755017,
                    "chr11": 121973369},
)

RAT = Species(
    scientific_name="Rattus norvegicus", taxon_id=10116, entry_suffix="RAT",
    common_name="rat", assembly="mRatBN7.2",
)

ZEBRAFISH = Species(
    scientific_name="Danio rerio", taxon_id=7955, entry_suffix="DANRE",
    common_name="zebrafish", assembly="GRCz11",
)

BUILTIN = {"human": HUMAN, "mouse": MOUSE, "rat": RAT,
           "zebrafish": ZEBRAFISH}


def species_dir() -> Path:
    """`config/species/` in the source tree, if it is there."""
    return Path(__file__).resolve().parents[2] / "config" / "species"


def available() -> list[str]:
    """Names `load_species` will accept, built-in and on disk."""
    names = set(BUILTIN)
    d = species_dir()
    if d.is_dir():
        names.update(p.stem for p in d.glob("*.yaml"))
    return sorted(names)


def _as_int(val, key: str) -> int:
    try:
        return int(str(val).strip())
    except (TypeError, ValueError):
        raise ValueError(f"species: {key} must be a whole number, "
                         f"got {val!r}") from None


def from_mapping(data: dict, where: str = "species") -> Species:
    """Build a Species from a parsed YAML mapping, rejecting unknown keys.

    Unknown keys are refused rather than ignored: a misspelled `taxon-id`
    that silently left the taxon at 9606 would put the wrong organism in
    every header of the run.
    """
    known = set(Species.__dataclass_fields__)
    unknown = sorted(set(data) - known)
    if unknown:
        raise ValueError(
            f"{where}: unknown key(s) {', '.join(unknown)}; "
            f"choose from {', '.join(sorted(known))}")

    lengths: dict[str, int] = {}
    raw = data.get("contig_lengths") or {}
    if not isinstance(raw, dict):
        raise ValueError(f"{where}: contig_lengths must be a mapping")
    for k, v in raw.items():
        lengths[str(k)] = _as_int(v, f"contig_lengths.{k}")

    mito = data.get("mito_contigs") or list(Species.mito_contigs)
    if isinstance(mito, str):
        mito = [x.strip() for x in mito.split(",") if x.strip()]

    code = _as_int(data.get("mito_code", Species.mito_code), "mito_code")
    if code not in SUPPORTED_MITO_CODES:
        # Accepting a code seqops cannot apply would put a setting in the
        # provenance record that never happened - the same failure as a
        # config file whose keys do nothing. Refuse instead, and name what
        # is implemented so the message is actionable.
        raise ValueError(
            f"{where}: mito_code {code} is not implemented. "
            f"seqops.py provides "
            f"{', '.join(f'{k} ({v})' for k, v in SUPPORTED_MITO_CODES.items())}. "
            f"Add the table to seqops.py before declaring it here.")

    return Species(
        scientific_name=str(data.get("scientific_name",
                                     Species.scientific_name)),
        taxon_id=_as_int(data.get("taxon_id", Species.taxon_id), "taxon_id"),
        entry_suffix=str(data.get("entry_suffix",
                                  Species.entry_suffix)).upper(),
        common_name=str(data.get("common_name", Species.common_name)),
        assembly=str(data.get("assembly", Species.assembly)),
        mito_contigs=tuple(str(x) for x in mito),
        mito_code=code,
        contig_lengths=lengths,
    )


def load_species(name) -> Species:
    """Resolve a species by name, or by path to a YAML file.

    `None` and `"human"` both give the human default, so an unspecified
    species behaves exactly as before this module existed.
    """
    if name is None:
        return HUMAN
    if isinstance(name, Species):        # tolerated for convenience
        return name

    text = str(name).strip()
    if not text:
        return HUMAN

    from .config import _mini_yaml       # the same tiny parser as run config

    p = Path(text)
    if p.suffix in (".yaml", ".yml") or p.is_file():
        if not p.is_file():
            raise ValueError(f"no species file at {p}")
        return from_mapping(_mini_yaml(p.read_text(encoding="utf-8"),
                                       str(p)), str(p))

    key = text.lower()
    on_disk = species_dir() / f"{key}.yaml"
    if on_disk.is_file():
        return from_mapping(
            _mini_yaml(on_disk.read_text(encoding="utf-8"), str(on_disk)),
            str(on_disk))
    if key in BUILTIN:
        return BUILTIN[key]
    raise ValueError(f"unknown species {text!r}; "
                     f"choose from {', '.join(available())} "
                     f"or give a path to a YAML file")
