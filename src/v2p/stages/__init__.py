"""The numbered pipeline stages.

Each stage is a standalone executable script rather than a function, so that
every stage records its own provenance JSON — inputs with checksums,
parameters, counters, environment — and so a stage can be re-run on its own
without re-running the pipeline. `v2p.cli` invokes them as subprocesses.

They live inside the package so that a `pip install v2p` can find them.
Their filenames start with digits, so they are not importable as modules and
were never meant to be: `v2p.cli.stage_dir()` locates this directory and
runs them by path.
"""

from __future__ import annotations

__all__ = ["STAGE_FILES"]

# The stages the CLI runs, in pipeline order. Listed explicitly so a missing
# or misnamed file is a clear error at packaging time rather than a
# subprocess failure halfway through a run.
STAGE_FILES = (
    "00_audit_references.py",
    "00_inspect_inputs.py",
    "01_parse_inputs.py",
    "02_build_protein_fasta.py",
    "03_qc_report.py",
    "04_validate_uniprot.py",
    "05_compare_tools.py",
    "06_package_release.py",
    "08_recovery_report.py",
)
