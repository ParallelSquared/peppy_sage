# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`peppy_sage` is a hybrid Rust/Python library for peptide scoring in mass spectrometry proteomics. It uses PyO3 to expose Rust implementations of high-performance peptide database indexing and spectrum scoring to Python, leveraging the `sage_core` library for core proteomics algorithms.

## Build and Development Commands

### Building the Project

```bash
# Build the Rust extension and install locally
maturin develop

# Build for release (optimized)
maturin develop --release

# Build wheels for distribution
maturin build --release
```

### Testing

```bash
# Run Rust unit tests
cargo test

# Run Python integration tests
python tests/python_test_search.py
python tests/python_library_test_search.py
```

### Development Workflow

```bash
# After making Rust changes, rebuild and reinstall
maturin develop

# Check Rust code
cargo check
cargo clippy
```

## Architecture

### Hybrid Rust/Python Design

The codebase follows a two-layer architecture:

1. **Rust Core Layer** (`src/`):
   - `lib.rs`: PyO3 bindings, wrapper structs (`PyPeptide`, `PyIndexedDatabase`, `PyProcessedSpectrum`, `PyScorer`, etc.)
   - `index_logic.rs`: Database indexing, fragment generation, peptide sorting
   - `scoring_logic.rs`: Spectrum scoring using sage_core's `Scorer`

2. **Python API Layer** (`python/peppy_sage/`):
   - `__init__.py`: Exports both low-level PyO3 types and high-level wrappers
   - `core.py`: High-level `Peptide`, `Spectrum`, `Precursor` classes
   - `indexing.py`: `IndexedDatabase` wrapper
   - `scoring.py`: `Scorer` wrapper

### Key Architectural Patterns

**GIL Management**: The codebase uses `Arc<T>` and `py.allow_threads()` to release the GIL during CPU-intensive operations (indexing, scoring). Rust data is cloned cheaply via `Arc::clone()` before releasing GIL.

**Zero-Copy Data Transfer**: Uses `pyo3-polars` for zero-copy DataFrame transfer between Python and Rust. Fragment results are returned as columnar `PyFeatureArrays` for efficient conversion to Polars/Pandas.

**Parallel Processing**: Heavy use of Rayon for parallel iteration in:
- Fragment generation (`par_iter`, `flat_map`)
- Peptide sorting and deduplication (`par_sort_unstable_by`, `par_chunks_mut`)
- Spectrum scoring (`par_iter` over spectra batches)

### Database Indexing

Two indexing modes supported:

1. **From Peptides** (`build_indexed_database`): Traditional FASTA-based workflow
   - Generates theoretical fragments from peptide sequences
   - Optional decoy generation by sequence reversal
   - Filters fragments by `min_ion_index` and mass range

2. **From Library** (`build_indexed_database_from_library`): Spectral library workflow
   - Reads Polars DataFrame with columns: `StrippedPeptide`, `ModifiedPeptide`, `FragmentMz`, `RelativeIntensity`, etc.
   - Aggregates fragments across precursor charge states
   - Stores library-specific metadata in `IndexedDatabase.library_frags`

**Fragment Storage**: Fragments are sorted by m/z and bucketed for efficient binary search during scoring. Each bucket is sorted by peptide index.

### Scoring Pipeline

1. **Input**: `PyProcessedSpectrum` (contains precursor info and MS2 peaks)
2. **Scorer Config**: `ScorerConfig` holds tolerances, charge ranges, scoring algorithm
3. **Execution**: `run_scoring()` creates a native `sage_core::Scorer` and calls `scorer.score(spectrum)`
4. **Output**: `PyFeatureArrays` with 38+ columnar fields (scores, masses, fragments, q-values)

### Modification Handling

Modifications are stored as `Vec<f32>` with layout: `[N-term, per-residue..., C-term]`

For library mode, UniMod strings in `ModifiedPeptide` (e.g., `AC(UniMod:4)DEFK`) are parsed via `unimod_table()` to numeric masses.

## Important Implementation Details

### Cargo Profile

The release profile in `Cargo.toml` has `debug = true` and `strip = "none"` to enable debugging while maintaining optimizations. **Do not remove these settings.**

### Module System

The compiled Rust extension is named `_peppy_sage` (with leading underscore) and lives inside the `peppy_sage` Python package. Python wrappers re-export functionality under cleaner names.

### Fragment Ordinal Filtering

`min_ion_index` skips N-terminal fragments (b1, b2) and C-terminal fragments (yn-1, yn-2) which are typically low-information. Logic differs by ion type:
- b/a/c ions: `(ion_idx + 1) > min_ion_index`
- y/x/z ions: `peptide.len() - 1 - ion_idx > min_ion_index`

### Peptide Deduplication

`sort_and_dedup()` merges peptides with identical sequence+mods+nterm+cterm, aggregating protein lists and preserving target status if any source was a target.

### Library Fragment Reordering

After building from library, `reorder_library_mode()` sorts peptides by monoisotopic mass, remaps all `fragment.peptide_index` values, and keeps `library_frags` aligned with the new peptide order.

## Dependencies

- **sage_core**: Sibling crate at `../sage_core`, provides core proteomics types (`Peptide`, `IndexedDatabase`, `Scorer`, etc.)
- **pyo3**: Python bindings (requires `abi3-py311` for stable ABI)
- **pyo3-polars**: Zero-copy DataFrame bridge
- **rayon**: Data parallelism
- **dashmap**: Concurrent hash set for decoy generation
- **polars**: DataFrame library

## Common Pitfalls

- **Don't forget `maturin develop`**: After changing Rust code, Python won't see changes until you rebuild
- **GIL context**: Never hold Python references (`Py<T>`) while releasing GIL via `allow_threads()`
- **Fragment sorting**: Fragments MUST be sorted by m/z before bucketing, or scoring will fail
- **Test data location**: Integration tests expect data files in `tests/` directory
