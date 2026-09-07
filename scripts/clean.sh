#!/usr/bin/env bash
# Remove intermediate and superseded files, keeping inputs, code and the
# release. Run from the project root.
#
#   bash scripts/clean.sh --dry-run     show what would go
#   bash scripts/clean.sh               delete it
#
# Never touches: data/, ref/, src/, scripts/, tests/, config/, docs/,
# release/, or the current results/ and logs/ trees.

set -euo pipefail
DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

say() { if [ "$DRY" = 1 ]; then echo "  would remove  $1"; else rm -rf "$1" && echo "  removed  $1"; fi; }

echo "== superseded patch scripts (already applied) =="
for f in apply_v11.py apply_v12.py apply_v13.py apply_v14.py; do
  [ -e "$f" ] && say "$f"
done

echo "== pre-patch backups =="
find src scripts tests -name "*.v1?.bak" 2>/dev/null | while read -r f; do say "$f"; done

echo "== python caches =="
find . -name "__pycache__" -type d 2>/dev/null | while read -r f; do say "$f"; done
find . -name "*.pyc" 2>/dev/null | while read -r f; do say "$f"; done

echo "== scratch and test fixtures (regenerated on demand) =="
for f in results/_dryrun results/_itest tests/fixtures/mini.fa tests/fixtures/mini.gtf \
         tests/fixtures/mini.fa.fai tests/fixtures/out.*.fasta patch_v1.1; do
  [ -e "$f" ] && say "$f"
done

# Logs are provenance: every stage records inputs, parameters, environment
# and checksums, so they are the reproducibility record. Only the runs that
# predate the current code are dead weight.
echo "== logs older than the newest stage-02 run =="
NEWEST=$(ls -t logs/02_build_protein_fasta.*.log 2>/dev/null | head -1 || true)
if [ -n "$NEWEST" ]; then
  STAMP=$(basename "$NEWEST" | sed 's/.*\.\([0-9T]*\)\.log/\1/')
  echo "  keeping runs at or after $STAMP"
  for f in logs/*.log logs/*.provenance.json; do
    [ -e "$f" ] || continue
    s=$(basename "$f" | sed 's/.*\.\([0-9]\{8\}T[0-9]\{6\}\)\..*/\1/')
    [ "$s" \< "$STAMP" ] && say "$f"
  done
fi

echo
echo "Kept: data/ ref/ src/ scripts/ tests/ config/ docs/ release/ results/ logs/(current)"
[ "$DRY" = 1 ] && echo "Dry run — nothing was deleted."
du -sh . 2>/dev/null | sed 's/^/Project size: /'
