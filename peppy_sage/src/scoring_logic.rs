use sage_core::mass::Tolerance;
use sage_core::scoring::{Scorer, ScoreType, Feature};
use sage_core::database::IndexedDatabase;
use sage_core::spectrum::ProcessedSpectrum;
use sage_core::spectrum::Peak;

/// Configuration Struct: Holds all necessary scoring parameters.
/// Mirrors the native `Scorer` struct but decoupled for Python configuration.
#[derive(Clone)]
pub struct ScorerConfig {
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub wide_window: bool,
    pub min_matched_peaks: u16,
    pub min_isotope_err: i8,
    pub max_isotope_err: i8,
    pub min_precursor_charge: u8,
    pub max_precursor_charge: u8,
    pub report_psms: usize,
    pub score_type: ScoreType,
    pub override_precursor_charge: bool,
    pub max_fragment_charge: Option<u8>,
    pub chimera: bool,
    pub annotate_matches: bool,
}

/// Main scoring function: builds a native Scorer and runs scoring on a spectrum.
///
/// This is called by your PyO3 wrapper (`PyScorer.score_spectra`).
pub fn run_scoring(
    config: &ScorerConfig,
    db: &IndexedDatabase,
    spectrum: &ProcessedSpectrum<Peak>,
) -> Vec<Feature> {
    // Construct the native Scorer
    let native_scorer = Scorer {
        db,
        precursor_tol: config.precursor_tol,
        fragment_tol: config.fragment_tol,
        min_matched_peaks: config.min_matched_peaks,
        min_isotope_err: config.min_isotope_err,
        max_isotope_err: config.max_isotope_err,
        min_precursor_charge: config.min_precursor_charge,
        max_precursor_charge: config.max_precursor_charge,
        wide_window: config.wide_window,
        chimera: config.chimera,
        report_psms: config.report_psms,
        score_type: config.score_type,
        override_precursor_charge: config.override_precursor_charge,
        max_fragment_charge: config.max_fragment_charge,
        annotate_matches: config.annotate_matches,
    };

    // Run the scoring (this executes sage_core's native logic)
    native_scorer.score(spectrum)
}