from typing import List, Optional, Tuple, Union
import numpy as np
from pyteomics import mass
from . import _rust

PROTON_MASS = 1.0072764


class Precursor:
    """
    Lightweight Python wrapper for PyPrecursor (MS1 parent ion).
    """

    def __init__(
            self,
            mz: float,
            charge: int,
            isolation_window: Union[_rust.PyTolerance, Tuple[float, float]] = (0.5, 0.5),
    ):
        """
        Args:
            mz: The m/z of the precursor ion.
            charge: The precursor charge state.
            isolation_window: Either a PyTolerance object or a (low, high) Da tuple.
        """
        if isinstance(isolation_window, tuple):
            isolation_window = _rust.PyTolerance.Da(*isolation_window)

        self._inner = _rust.PyPrecursor(
            mz=mz,
            charge=charge,
            isolation_window=isolation_window,
        )

    @property
    def mz(self) -> float:
        return self._inner.mz

    @property
    def charge(self) -> int:
        return self._inner.charge

    @property
    def isolation_window(self) -> _rust.PyTolerance:
        return self._inner.isolation_window

    def __repr__(self) -> str:
        return f"<Precursor mz={self.mz:.4f}, z={self.charge}>"


class Spectrum:
    """
    High-level Python wrapper around Rust's PyProcessedSpectrum.
    """

    def __init__(
            self,
            id: str,
            file_id: int,
            scan_start_time: float,
            mz_array: Union[List[float], np.ndarray],
            intensity_array: Union[List[float], np.ndarray],
            precursors: List[Precursor],
            total_ion_current: Optional[float] = None,
    ):
        """
        Create a spectrum from arrays of m/z and intensities.

        Args:
            id: Spectrum identifier (e.g. 'Scan_1234')
            file_id: Numeric file index (if multiple files are loaded)
            scan_start_time: Retention time (in minutes)
            mz_array: List or numpy array of fragment m/z values
            intensity_array: List or numpy array of corresponding intensities
            precursors: List of Precursor objects
            total_ion_current: Optional total ion current; computed automatically if None
        """

        mz_array = np.asarray(mz_array, dtype=float)
        intensity_array = np.asarray(intensity_array, dtype=float)

        if mz_array.shape != intensity_array.shape:
            raise ValueError("mz_array and intensity_array must be the same length.")

        if total_ion_current is None:
            total_ion_current = float(np.sum(intensity_array))

        mz_array   = np.asarray(mz_array, dtype=np.float32, order="C")
        intensity_array = np.asarray(intensity_array, dtype=np.float32, order="C")

        self._inner = _rust.PyProcessedSpectrum(
            id,
            file_id,
            scan_start_time,
            mz_array,
            intensity_array,
            [p._inner for p in precursors],
            total_ion_current,
        )

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------
    @property
    def id(self) -> str:
        return self._inner.id

    @property
    def scan_start_time(self) -> float:
        return self._inner.scan_start_time

    @property
    def precursors(self) -> List[Precursor]:
        """Return a list of precursor ions."""
        return [Precursor(p.mz, p.charge, p.isolation_window) for p in self._inner.precursors]

    @property
    def peaks(self) -> List[Tuple[float, float]]:
        """Return list of (m/z, intensity) tuples."""
        return [(m, i) for (m, i) in self._inner.peaks]

    @property
    def total_ion_current(self) -> float:
        return self._inner.total_ion_current

    # -------------------------------------------------------------------------
    # Utilities
    # -------------------------------------------------------------------------
    def sort_peaks(self, inplace: bool = True) -> "Spectrum":
        """Sort peaks by m/z value."""
        peaks = sorted(self.peaks, key=lambda x: x[0])
        if inplace:
            mz, ints = zip(*peaks)
            self._inner = _rust.PyProcessedSpectrum(
                self.id,
                self._inner.file_id,
                self.scan_start_time,
                list(mz),
                list(ints),
                [p._inner for p in self.precursors],
                self.total_ion_current,
            )
            return self
        else:
            mz, ints = zip(*peaks)
            return Spectrum(
                id=self.id,
                file_id=self._inner.file_id,
                scan_start_time=self.scan_start_time,
                mz_array=mz,
                intensity_array=ints,
                precursors=self.precursors,
                total_ion_current=self.total_ion_current,
            )

    def __repr__(self):
        n_peaks = len(self._inner.peaks)
        return f"<Spectrum id={self.id!r}, peaks={n_peaks}, precursors={len(self.precursors)}>"


class Peptide:
    """
    High-level Python wrapper for a Rust PyPeptide.
    Handles sequence, modifications, and mass calculation.
    """

    def __init__(self, sequence: str, mods: List[float] | None = None): #TODO handle n-term and c-term
        """
        Parameters
        ----------
        sequence : str
            Peptide amino acid sequence (e.g., "PEPTIDEK")
        mods : List[float], optional
            Array of modification masses (length of sequence + 2).
            Default is zeros (no mods).
        """
        self.sequence = sequence

        if mods is None:
            mods = [0.0] * (len(sequence) + 2)
        elif len(mods) != (len(sequence) + 2):
            raise ValueError(
                f"Modification array must have length sequence ({len(sequence)}) + 2."
            )

        self.mods = mods

        # --- Calculate precursor mass using pyteomics ---
        self.monoisotopic_mass = mass.calculate_mass(sequence=sequence)
        self.monoisotopic_mass += sum(mods)

        # --- Create the Rust-side PyPeptide object ---
        self._inner = _rust.PyPeptide(sequence, self.monoisotopic_mass, mods[1:-1], mods[0], mods[-1])

    @classmethod
    def from_rust(cls, rust_peptide: _rust.PyPeptide) -> "Peptide":
        """Construct a Python Peptide from a Rust PyPeptide object."""
        sequence = rust_peptide.sequence
        # Rust mods are floats
        mods = rust_peptide.modifications.copy()
        return cls(sequence, mods)

    def calculate_theoretical_mz(self, z=1):
        return (self.monoisotopic_mass + (PROTON_MASS * z)) / z

    @staticmethod
    def calculate_theoretical_mz(seq, mods, z) -> float:
        """
        Compute the theoretical m/z for a peptide.

        Parameters
        ----------
        monoisotopic_mass : float
            Mass of the peptide (including modifications if any)
        charge : int
            Charge state

        Returns
        -------
        float
            Theoretical m/z
        """
        monoisotopic_mass = mass.calculate_mass(sequence=seq)
        monoisotopic_mass += sum(mods)

        return (monoisotopic_mass + PROTON_MASS * z) / z

    def __repr__(self):
        mods_str = ", ".join(f"{m:.2f}" for m in self.mods)
        return (
            f"<Peptide seq='{self.sequence}' "
            f"mass={self.monoisotopic_mass:.4f} "
            f"mods=[{mods_str}]>"
        )