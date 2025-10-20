from peppy_sage import PyScorer, PyTolerance

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
        features = self._scorer.score_spectra(db, spectrum)

        print(f"\n--- Scoring Results ---")
        print(f"Features returned: {len(features)}")

        for i, f in enumerate(features):
            pep = db.peptides[f.peptide_idx]
            print(f"\n→ Peptide: {pep.sequence}")
            if hasattr(f, "fragments") and f.fragments is not None:
                print("Matched fragments:")
                print(f.fragments.to_dict())
            print(
                f"Feature {i}: hyperscore={f.hyperscore:.2f}, "
                f"delta_mass={f.delta_mass:.4f}, matched_peaks={f.matched_peaks}"
            )

        return features