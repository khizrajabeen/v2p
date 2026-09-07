"""HCC1395 variant-to-protein conversion toolkit."""

# Single source of truth for the version.
#
# pyproject.toml reads this line -- and only this line -- via
# [tool.setuptools.dynamic] version = {attr = "v2p.__version__"}. setuptools
# resolves that from the module's AST without importing the package, which is
# why this assignment must stay a plain literal: no f-string, no
# importlib.metadata lookup, no computation. Anything cleverer here and the
# build backend falls back to importing v2p, which drags in the whole
# dependency chain at build time.
#
# Do not read the version from importlib.metadata either. The package is
# routinely used straight from a source checkout (tests/ and scripts/v2p both
# put src/ on sys.path rather than installing), where no distribution
# metadata exists.
__version__ = "1.0.0"

# Parsed form, for comparisons that should not be string comparisons.
__version_info__ = tuple(int(part) for part in __version__.split("."))

__all__ = ["__version__", "__version_info__"]
