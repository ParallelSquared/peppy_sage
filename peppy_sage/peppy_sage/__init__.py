# peppy_sage/__init__.py

# Import the compiled PyO3 extension module
# (the same name as your PyO3 crate, built by maturin)
from .peppy_sage import (
    PyPeptide,
    PyIndexedDatabase,
    PyProcessedSpectrum,
    PyScorer,
    PyTolerance,
    PyPrecursor,
    PyKind,
)

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