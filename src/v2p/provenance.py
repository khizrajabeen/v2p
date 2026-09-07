"""Run provenance: logging, input checksums, environment capture.

Every stage of the pipeline opens a RunLogger. It writes a human-readable
log to logs/<stage>.<run_id>.log and a machine-readable JSON sidecar so a
third party can reproduce the run byte-for-byte.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


def sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    """SHA-256 of a file, streamed so large references don't blow up RAM."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def _pip_freeze() -> list[str]:
    try:
        out = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True, text=True, timeout=120, check=False,
        )
        return sorted(x for x in out.stdout.splitlines() if x.strip())
    except Exception:  # pragma: no cover - environment dependent
        return []


def _git_commit(repo_root: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=20, check=False,
        )
        return out.stdout.strip() or None
    except Exception:  # pragma: no cover
        return None


@dataclass
class RunLogger:
    """Structured logger + provenance recorder for one pipeline stage."""

    stage: str
    log_dir: Path
    run_id: str = field(default_factory=lambda: time.strftime("%Y%m%dT%H%M%S"))
    _inputs: dict[str, dict[str, Any]] = field(default_factory=dict, init=False)
    _outputs: dict[str, dict[str, Any]] = field(default_factory=dict, init=False)
    _params: dict[str, Any] = field(default_factory=dict, init=False)
    _counters: dict[str, int] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        self.log_dir = Path(self.log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_path = self.log_dir / f"{self.stage}.{self.run_id}.log"
        self.json_path = self.log_dir / f"{self.stage}.{self.run_id}.provenance.json"
        self._t0 = time.time()

        self.log = logging.getLogger(f"v2p.{self.stage}.{self.run_id}")
        self.log.setLevel(logging.DEBUG)
        self.log.handlers.clear()
        self.log.propagate = False
        fmt = logging.Formatter(
            "%(asctime)s\t%(levelname)-7s\t%(message)s", "%Y-%m-%dT%H:%M:%S"
        )
        fh = logging.FileHandler(self.log_path, encoding="utf-8")
        fh.setFormatter(fmt)
        fh.setLevel(logging.DEBUG)
        sh = logging.StreamHandler(sys.stderr)
        sh.setFormatter(fmt)
        sh.setLevel(logging.INFO)
        self.log.addHandler(fh)
        self.log.addHandler(sh)

        self.log.info("stage=%s run_id=%s", self.stage, self.run_id)
        self.log.info("python=%s", sys.version.replace("\n", " "))
        self.log.info("platform=%s", platform.platform())

    # -- recording -------------------------------------------------------
    def add_input(self, role: str, path: str | Path, checksum: bool = True) -> Path:
        p = Path(path)
        rec: dict[str, Any] = {"path": str(p.resolve()), "exists": p.exists()}
        if p.exists():
            rec["bytes"] = p.stat().st_size
            if checksum and p.stat().st_size < 4 * (1 << 30):
                rec["sha256"] = sha256(p)
        self._inputs[role] = rec
        self.log.info(
            "input[%s]=%s bytes=%s sha256=%s",
            role, rec["path"], rec.get("bytes"), rec.get("sha256", "skipped")[:16],
        )
        return p

    def add_output(self, role: str, path: str | Path) -> Path:
        p = Path(path)
        rec: dict[str, Any] = {"path": str(p.resolve()), "exists": p.exists()}
        if p.exists():
            rec["bytes"] = p.stat().st_size
            rec["sha256"] = sha256(p)
        self._outputs[role] = rec
        self.log.info("output[%s]=%s bytes=%s", role, rec["path"], rec.get("bytes"))
        return p

    def add_params(self, **kw: Any) -> None:
        self._params.update(kw)
        for k, v in kw.items():
            self.log.info("param %s=%r", k, v)

    def count(self, key: str, n: int = 1) -> None:
        self._counters[key] = self._counters.get(key, 0) + n

    def counters(self) -> dict[str, int]:
        return dict(self._counters)

    # -- finalisation ----------------------------------------------------
    def close(self, status: str = "ok") -> Path:
        for k in sorted(self._counters):
            self.log.info("count %s=%d", k, self._counters[k])
        elapsed = round(time.time() - self._t0, 3)
        doc = {
            "stage": self.stage,
            "run_id": self.run_id,
            "status": status,
            "elapsed_seconds": elapsed,
            "started_utc": time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", time.gmtime(self._t0)
            ),
            "command_line": sys.argv,
            "cwd": os.getcwd(),
            "python": sys.version,
            "platform": platform.platform(),
            "git_commit": _git_commit(Path(__file__).resolve().parents[2]),
            "params": self._params,
            "inputs": self._inputs,
            "outputs": self._outputs,
            "counters": self._counters,
            "pip_freeze": _pip_freeze(),
            "log_file": str(self.log_path),
        }
        self.json_path.write_text(json.dumps(doc, indent=2), encoding="utf-8")
        self.log.info("status=%s elapsed=%.3fs provenance=%s",
                      status, elapsed, self.json_path)
        for h in list(self.log.handlers):
            h.flush()
            h.close()
            self.log.removeHandler(h)
        return self.json_path
