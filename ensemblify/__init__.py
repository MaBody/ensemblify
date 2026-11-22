from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

ROOT = Path(__file__).resolve().parent

try:
    __version__ = version("ensemblify")
except PackageNotFoundError:  # running from a source tree without installing
    __version__ = "unknown"
