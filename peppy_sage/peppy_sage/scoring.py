from peppy_sage import PyScorer, PyTolerance, Peptide

class Scorer:
    """
    High-level scoring interface for peppy_sage.
    Wraps the underlying Rust-based PyScorer for convenience.
    """

    def __init__(
            self,
            precursor_tol_da: tuple[float, float] = (-1.0, 1.0),
            fragment_tol_ppm: tuple[float, float] = (-5.0, 5.0),
            wide_window: bool = False,
            chimera: bool = True,
            report_psms: int = 10,
            min_isotope_err: int = -1,
            max_isotope_err: int = 3,
            min_precursor_charge: int = 1,
            max_precursor_charge: int = 4,
            min_matched_peaks: int = 0,
            annotate_matches: bool = True,
            max_fragment_charge: int = 1,
    ):
        """Initialize the high-level scoring object."""
        precursor_tol = PyTolerance.Da(*precursor_tol_da)
        fragment_tol = PyTolerance.Ppm(*fragment_tol_ppm)

        self._scorer = PyScorer(
            precursor_tol,
            fragment_tol,
            wide_window,
            chimera,
            report_psms,
            min_isotope_err,
            max_isotope_err,
            min_precursor_charge,
            max_precursor_charge,
            min_matched_peaks,
            annotate_matches,
            max_fragment_charge,
        )

    def score(self, db, spectrum):
        """
        Score a single spectrum (or a list of spectra) against the given database.
        Returns a list of features (PSMs).
        """
        features = self._scorer.score_spectra(db._inner, spectrum._inner)

        return features

    def score_many(self, db, spectra):
        """
        Score many spectra and convert Rust peptides to Python Peptide.
        """
        rust_spectra = [s._inner for s in spectra]
        all_hits = self._scorer.score_many_spectra(db._inner, rust_spectra)

        return all_hits