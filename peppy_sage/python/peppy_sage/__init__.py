# peppy_sage/__init__.py

# Import the compiled PyO3 extension module
# (the same name as your PyO3 crate, built by maturin)
from ._peppy_sage import (
    PyPeptide,
    PyIndexedDatabase,
    PyProcessedSpectrum,
    PyScorer,
    PyTolerance,
    PyPrecursor,
    PyKind,
)

# Define the _rust alias here for internal use by other modules.
# This makes it accessible via 'from . import _rust' in other package files.
_rust = type('RustBindings', (object,), {
    'PyPeptide': PyPeptide,
    'PyIndexedDatabase': PyIndexedDatabase,
    'PyProcessedSpectrum': PyProcessedSpectrum,
    'PyScorer': PyScorer,
    'PyTolerance': PyTolerance,
    'PyPrecursor': PyPrecursor,
    'PyKind': PyKind,
})

# Import the high-level Python wrappers
from .core import Peptide, Spectrum, Precursor
from .indexing import IndexedDatabase
from .scoring import Scorer

__all__ = [
    # Low-level bindings
    "PyPeptide",
    "PyIndexedDatabase",
    "PyProcessedSpectrum",
    "PyScorer",
    "PyTolerance",
    "PyPrecursor",
    "PyKind",
    # High-level Python API
    "Peptide",
    "IndexedDatabase",
    "Spectrum",
    "Precursor",
    "Scorer",
]

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("peppy_sage")
except PackageNotFoundError:
    # Package is not installed (e.g. running from source tree)
    __version__ = "0.0.0"
