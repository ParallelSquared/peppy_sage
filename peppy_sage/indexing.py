# indexed_database.py

from typing import List, Optional, Tuple
import peppy_sage as _rust  # your compiled PyO3 module

ION_KIND_MAP = {
    "a": _rust.PyKind.A,
    "b": _rust.PyKind.B,
    "c": _rust.PyKind.C,
    "x": _rust.PyKind.X,
    "y": _rust.PyKind.Y,
    "z": _rust.PyKind.Z,
}

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
            peptides: List["_rust.PyPeptide"],
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
            ion_enum_list = [
                ION_KIND_MAP[k.lower()]
                for k in ion_kinds
            ]
        except KeyError as e:
            raise ValueError(f"Invalid ion kind: {e.args[0]}. Must be one of {list(ION_KIND_MAP)}")

        inner = _rust.PyIndexedDatabase.from_peptides_and_config(
            peptides=peptides,
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
