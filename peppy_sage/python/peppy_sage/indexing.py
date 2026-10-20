# indexed_database.py

from typing import List, Optional, Tuple, Union
from .core import Peptide
from . import _rust
import polars as pl

ION_KIND_MAP = {
    #"a": _rust.PyKind.A, #TODO other ion types
    "b": _rust.PyKind.B,
    #"c": _rust.PyKind.C,
    #"x": _rust.PyKind.X,
    "y": _rust.PyKind.Y,
    #"z": _rust.PyKind.Z,
}

def _to_rust_peptide(x: Union[Peptide, _rust.PyPeptide, str]) -> _rust.PyPeptide:
    if isinstance(x, _rust.PyPeptide):
        return x
    if isinstance(x, Peptide):
        # assuming your wrapper keeps the raw object on ._inner (or similar)
        return x._inner
    if isinstance(x, str):
        return _rust.PyPeptide(x)
    raise TypeError(
        f"Expected Peptide | PyPeptide | str, got {type(x).__name__}"
    )

class IndexedDatabase:
    """
    Pythonic wrapper around the Rust-backed PyIndexedDatabase.

    Example:
        >>> db = IndexedDatabase.from_peptides(
        ...     peptides=my_peptides,
        ...     ion_kinds=["b", "y"],
        ...     bucket_size=50,
        ...     generate_decoys=True,
        ...     decoy_tag="DECOY_"
        ... )
        >>> print(len(db.peptides))
        >>> print(db.fragment_count)
    """

    def __init__(self, _inner: _rust.PyIndexedDatabase):
        self._inner = _inner

    # -------------------------------------------------------------------------
    # Constructors
    # -------------------------------------------------------------------------
    @classmethod
    def from_peptides(
            cls,
            peptides: List[Union[Peptide, "_rust.PyPeptide"]],
            ion_kinds: List[str],
            bucket_size: int = 8192,
            min_ion_index: int = 1,
            generate_decoys: bool = False,
            decoy_tag: str = "rev_",
            peptide_min_mass: float = 500,
            peptide_max_mass: float = 5000.0,
    ) -> "IndexedDatabase":
        """
        Build an indexed database from a list of peptides and configuration.

        Args:
            peptides: List of PyPeptide objects (from Rust layer).
            ion_kinds: List of ion types to include, e.g. ["b", "y"].
            bucket_size: Number of fragments per mass bucket.
            min_ion_index: Minimum ion index to include.
            generate_decoys: Whether to generate decoy peptides.
            decoy_tag: Prefix used to label decoy peptides.
            peptide_min_mass: Minimum peptide monoisotopic mass to include.
            peptide_max_mass: Maximum peptide monoisotopic mass to include.
        """

        try:
            ion_enum_list = [ION_KIND_MAP[k.lower()]() for k in ion_kinds]
        except KeyError as e:
            raise ValueError(f"Invalid ion kind: {e.args[0]}. Must be one of {list(ION_KIND_MAP)}")

        rust_peptides = [_to_rust_peptide(p) for p in peptides]

        inner = _rust.PyIndexedDatabase.from_peptides_and_config(
            peptides=rust_peptides,
            bucket_size=bucket_size,
            ion_kinds=ion_enum_list,
            min_ion_index=min_ion_index,
            generate_decoys=generate_decoys,
            decoy_tag=decoy_tag,
            peptide_min_mass=peptide_min_mass,
            peptide_max_mass=peptide_max_mass,
        )

        return cls(inner)


    @classmethod
    def from_library(
            cls,
            library: pl.DataFrame,
            ion_kinds: List[str],
            bucket_size: int = 8192,
            min_ion_index: int = 1,
            generate_decoys: bool = False,
            decoy_tag: str = "rev_",
            peptide_min_mass: float = 500.0,
            peptide_max_mass: float = 5000.0,
    ) -> "IndexedDatabase":
        """
        Build an indexed database directly from a DIA-NN-style library.

        Args:
            library: Polars DataFrame with columns like:
                - "StrippedPeptide", "ModifiedPeptide"
                - "FragmentMz", "FragmentType", "FragmentCharge", "FragmentSeriesNumber"
                - etc.
            ion_kinds: List of ion types to include, e.g. ["b", "y"].
            Other args mirror `from_peptides`.
        """
        if not isinstance(library, pl.DataFrame):
            raise TypeError(
                f"`library` must be a polars.DataFrame, got {type(library)!r}"
            )

        try:
            ion_enum_list = [ION_KIND_MAP[k.lower()]() for k in ion_kinds]
        except KeyError as e:
            raise ValueError(
                f"Invalid ion kind: {e.args[0]}. Must be one of {list(ION_KIND_MAP)}"
            )

        # This call is where we keep things zero-copy: `library` is passed
        # directly to Rust as a Polars / Arrow DataFrame (via pyo3-polars).
        inner = _rust.PyIndexedDatabase.from_library_and_config(
            library_df=library,
            bucket_size=bucket_size,
            ion_kinds=ion_enum_list,
            min_ion_index=min_ion_index,
            generate_decoys=generate_decoys,
            decoy_tag=decoy_tag,
            peptide_min_mass=peptide_min_mass,
            peptide_max_mass=peptide_max_mass,
        )

        return cls(inner)

    @classmethod
    def from_columnar(
            cls,
            seq,
            mod_seq,
            mods_values,
            mods_offsets,
            frag_mz,
            frag_intensity,
            frag_codes,
            frag_offsets,
            frag_lengths,
            ion_kinds: List[str],
            bucket_size: int = 8192,
            min_ion_index: int = 1,
            generate_decoys: bool = False,
            decoy_tag: str = "rev_",
            peptide_min_mass: float = 500.0,
            peptide_max_mass: float = 5000.0,
    ) -> "IndexedDatabase":
        """Build from columnar (CSR) arrays, without a one-row-per-fragment frame.

        An alternate entry point to `from_library`, which remains the general
        DataFrame route. Use this when you already hold a columnar library: it
        avoids materializing the long-form intermediate, where every precursor
        field is repeated across all of that precursor's fragments. On a large
        library that intermediate is the peak-memory event, not the ingest.

        Produces a database identical to `from_library` on the equivalent frame.

        Args:
            seq: length-N sequence of stripped peptide strings.
            mod_seq: length-N sequence of modified peptide strings. Entries
                sharing a modified peptide are merged into one peptide, with
                their fragments averaged across charge states.
            mods_values: float32 CSR values -- per-residue modification masses,
                concatenated. Entry i spans mods_offsets[i]:mods_offsets[i+1]
                and must be len(seq[i]) + 2 long (N-term, residues..., C-term).
            mods_offsets: int64, length N+1.
            frag_mz: float32, length M.
            frag_intensity: float32, length M.
            frag_codes: int32, length M. Packed as ion_type 3 bits at 0,
                ordinal 8 bits at 3, loss 3 bits at 11, charge-1 3 bits at 14,
                isotope 4 bits at 17.
            frag_offsets: int64, length N. Entry i's fragments start here.
            frag_lengths: int32, length N. Offsets are explicit, so the
                fragment arrays may contain gaps.

        Arrays are borrowed, not copied, so they must be C-contiguous and have
        exactly the dtypes above.
        """
        try:
            ion_enum_list = [ION_KIND_MAP[k.lower()]() for k in ion_kinds]
        except KeyError as e:
            raise ValueError(
                f"Invalid ion kind: {e.args[0]}. Must be one of {list(ION_KIND_MAP)}"
            )

        import numpy as np

        def _arr(a, dtype, name):
            out = np.ascontiguousarray(a, dtype=dtype)
            if out.ndim != 1:
                raise ValueError(f"{name} must be 1-dimensional, got shape {out.shape}")
            return out

        inner = _rust.PyIndexedDatabase.from_columnar_and_config(
            seq=list(seq),
            mod_seq=list(mod_seq),
            mods_values=_arr(mods_values, np.float32, "mods_values"),
            mods_offsets=_arr(mods_offsets, np.int64, "mods_offsets"),
            frag_mz=_arr(frag_mz, np.float32, "frag_mz"),
            frag_intensity=_arr(frag_intensity, np.float32, "frag_intensity"),
            frag_codes=_arr(frag_codes, np.int32, "frag_codes"),
            frag_offsets=_arr(frag_offsets, np.int64, "frag_offsets"),
            frag_lengths=_arr(frag_lengths, np.int32, "frag_lengths"),
            bucket_size=bucket_size,
            ion_kinds=ion_enum_list,
            min_ion_index=min_ion_index,
            generate_decoys=generate_decoys,
            decoy_tag=decoy_tag,
            peptide_min_mass=peptide_min_mass,
            peptide_max_mass=peptide_max_mass,
        )

        return cls(inner)

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------
    @property
    def peptides(self) -> List["_rust.PyPeptide"]:
        """Return the list of peptides stored in this indexed database."""
        return self._inner.peptides

    @property
    def fragment_count(self) -> int:
        """Return the total number of fragments indexed."""
        return self._inner.fragment_count()

    # -------------------------------------------------------------------------
    # Debug / Inspection
    # -------------------------------------------------------------------------
    def debug_fragment_summary(self) -> List[Tuple[float, int]]:
        """
        Return a summary of (fragment_mz, peptide_index) tuples for debugging.
        """
        return self._inner.debug_fragment_summary()

    def summary(self, limit: Optional[int] = 10) -> None:
        """Print a short summary of fragments for inspection."""
        summary = self.debug_fragment_summary()
        print(f"Indexed fragments: {self.fragment_count}")
        print("Example fragments:")
        for mz, pep_idx in summary[:limit]:
            print(f"  - m/z={mz:.4f}, peptide_index={pep_idx}")
