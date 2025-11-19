// peppysage/src/lib.rs

// --- MODULES ---
pub mod index_logic;
pub mod scoring_logic;

use crate::index_logic::{build_indexed_database, IndexingConfig};
use crate::scoring_logic::{ScorerConfig, run_scoring};

// --- CORE AND PY03 IMPORTS ---
use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use pyo3::types::PyString;
use pyo3::types::PyDict;
use pyo3::buffer::PyBuffer;
use pyo3::types::PyAny;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::sync::Arc;

use std::fs::File;


// Import ONLY the native Rust types from sage_core
use sage_core::database::{IndexedDatabase, PeptideIx, Theoretical};
use sage_core::enzyme::Position;
use sage_core::peptide::Peptide;
use sage_core::ion_series::Kind;
use sage_core::spectrum::{ProcessedSpectrum, Peak, Precursor};
use sage_core::mass::{Tolerance};
use sage_core::scoring::{Scorer, ScoreType, Feature, Fragments}; // Required for PyScorer/PyFeature
use sage_core::modification::ModificationSpecificity;

// --- 1. DEFINE WRAPPER STRUCTS ---

// Define PyPeptide, PyKind, PyModificationSpecificity, etc., first.

#[pyclass(module = "peppy_sage")]
pub struct PyIndexedDatabase {
    pub inner: Arc<IndexedDatabase>,
}

#[pyclass(module = "peppy_sage")]
#[derive(Clone)]
pub struct PyPeptide {
    pub inner: Peptide,
}

#[pyclass(module = "peppy_sage")]
#[derive(Clone, Copy)] // Clone/Copy is usually required if the native Kind enum is Copy
pub struct PyKind {
    pub inner: Kind, // Wraps sage_core::ion_series::Kind
}

#[pyclass(module = "peppy_sage")]
#[derive(Copy, Clone)]
pub struct PyTolerance {
    pub inner: Tolerance,
}

#[pyclass(module = "peppy_sage")]
#[derive(Clone)]
pub struct PyFragments {
    pub inner: Fragments,
}

#[pyclass(module = "peppy_sage")]
#[derive(Clone)]
pub struct PyProcessedSpectrum {
    pub inner: Arc<ProcessedSpectrum<Peak>>,
}

#[pyclass(module = "peppy_sage")]
#[derive(Clone)]
pub struct PyPrecursor {
    pub inner: Precursor,
}

#[pyclass(module = "peppy_sage")]
pub struct PyScorer {
    // Note: Scorer requires a lifetime db, but PyO3 often simplifies this
    // for objects passed temporarily. We wrap the configuration/parameters here.
    // For now, we only store the configuration needed to instantiate the actual Scorer later.
    pub config: ScorerConfig,
}

#[pyclass(module = "peppy_sage")]
#[derive(Clone)]
pub struct PyFeature {
    pub inner: Feature,
    pub peptide: Option<PyPeptide>,
}

/// A columnar representation of many PyFeature objects, optimized for zero-copy
/// transfer to Python's data analysis libraries (e.g., Polars or Pandas).
#[pyclass]
#[derive(Default, Clone, Debug)]
pub struct PyFeatureArrays {
    // Basic identifiers
    #[pyo3(get)]
    pub file_id: Vec<usize>,
    #[pyo3(get)]
    pub spec_id: Vec<String>,
    #[pyo3(get)]
    pub psm_id: Vec<usize>,
    #[pyo3(get)]
    pub rank: Vec<u32>,

    // Peptide info (extracted)
    #[pyo3(get)]
    pub sequence: Vec<String>,
    // Modifications are complex (list of floats). Using a list of lists (Vec<Vec<f32>>)
    // is the most direct way to transfer this, though Polars might prefer
    // to map this column to a List<Float32> type.
    #[pyo3(get)]
    pub modifications: Vec<Vec<f32>>,
    #[pyo3(get)]
    pub peptide_len: Vec<usize>,
    #[pyo3(get)]
    pub missed_cleavages: Vec<u8>,

    // Mass and error metrics
    #[pyo3(get)]
    pub expmass: Vec<f32>,
    #[pyo3(get)]
    pub calcmass: Vec<f32>,
    #[pyo3(get)]
    pub charge: Vec<u8>,
    #[pyo3(get)]
    pub delta_mass: Vec<f32>,
    #[pyo3(get)]
    pub isotope_error: Vec<f32>,
    #[pyo3(get)]
    pub average_ppm: Vec<f32>,

    // Score metrics
    #[pyo3(get)]
    pub hyperscore: Vec<f64>,
    #[pyo3(get)]
    pub delta_next: Vec<f64>,
    #[pyo3(get)]
    pub delta_best: Vec<f64>,
    #[pyo3(get)]
    pub poisson: Vec<f64>,
    #[pyo3(get)]
    pub discriminant_score: Vec<f32>,

    // Spectrum matching metrics
    #[pyo3(get)]
    pub matched_peaks: Vec<u32>,
    #[pyo3(get)]
    pub matched_intensity_pct: Vec<f32>,
    #[pyo3(get)]
    pub longest_b: Vec<u32>,
    #[pyo3(get)]
    pub longest_y: Vec<u32>,
    #[pyo3(get)]
    pub longest_y_pct: Vec<f32>,
    #[pyo3(get)]
    pub scored_candidates: Vec<u32>,
    #[pyo3(get)]
    pub ms2_intensity: Vec<f32>,

    // RT/IMS metrics
    #[pyo3(get)]
    pub rt: Vec<f32>,
    #[pyo3(get)]
    pub aligned_rt: Vec<f32>,
    #[pyo3(get)]
    pub predicted_rt: Vec<f32>,
    #[pyo3(get)]
    pub delta_rt_model: Vec<f32>,
    #[pyo3(get)]
    pub ims: Vec<f32>,
    #[pyo3(get)]
    pub predicted_ims: Vec<f32>,
    #[pyo3(get)]
    pub delta_ims_model: Vec<f32>,

    // Q-values and label
    #[pyo3(get)]
    pub label: Vec<i32>,
    #[pyo3(get)]
    pub spectrum_q: Vec<f32>,
    #[pyo3(get)]
    pub peptide_q: Vec<f32>,
    #[pyo3(get)]
    pub protein_q: Vec<f32>,
    #[pyo3(get)]
    pub posterior_error: Vec<f32>,

    //Fragment-level metrics (one list per feature)
    #[pyo3(get)]
    pub frag_charges: Vec<Vec<i32>>,
    #[pyo3(get)]
    pub frag_kinds: Vec<Vec<String>>,
    #[pyo3(get)]
    pub frag_fragment_ordinals: Vec<Vec<i32>>,
    #[pyo3(get)]
    pub frag_intensities: Vec<Vec<f32>>,
    #[pyo3(get)]
    pub frag_mz_calculated: Vec<Vec<f32>>,
    #[pyo3(get)]
    pub frag_mz_experimental: Vec<Vec<f32>>,
}

/// Iterates over the raw features, looks up the associated Peptide data,
/// and collects all fields into columnar vectors.
fn features_to_arrays(
    // The result of score_many_spectra, flattened from Vec<Vec<Feature>> to Vec<Feature>
    all_features: Vec<Feature>,
    db_arc: Arc<IndexedDatabase>,
) -> PyFeatureArrays {
    let mut arrays = PyFeatureArrays::default();
    let total_features = all_features.len();

    // We reserve capacity up front for better performance (all 38 fields)
    arrays.file_id.reserve(total_features);
    arrays.spec_id.reserve(total_features);
    arrays.psm_id.reserve(total_features);
    arrays.rank.reserve(total_features);
    arrays.sequence.reserve(total_features);
    arrays.modifications.reserve(total_features);
    arrays.peptide_len.reserve(total_features);
    arrays.missed_cleavages.reserve(total_features);
    arrays.expmass.reserve(total_features);
    arrays.calcmass.reserve(total_features);
    arrays.charge.reserve(total_features);
    arrays.delta_mass.reserve(total_features);
    arrays.isotope_error.reserve(total_features);
    arrays.average_ppm.reserve(total_features);
    arrays.hyperscore.reserve(total_features);
    arrays.delta_next.reserve(total_features);
    arrays.delta_best.reserve(total_features);
    arrays.poisson.reserve(total_features);
    arrays.discriminant_score.reserve(total_features);
    arrays.matched_peaks.reserve(total_features);
    arrays.matched_intensity_pct.reserve(total_features);
    arrays.longest_b.reserve(total_features);
    arrays.longest_y.reserve(total_features);
    arrays.longest_y_pct.reserve(total_features);
    arrays.scored_candidates.reserve(total_features);
    arrays.ms2_intensity.reserve(total_features);
    arrays.rt.reserve(total_features);
    arrays.aligned_rt.reserve(total_features);
    arrays.predicted_rt.reserve(total_features);
    arrays.delta_rt_model.reserve(total_features);
    arrays.ims.reserve(total_features);
    arrays.predicted_ims.reserve(total_features);
    arrays.delta_ims_model.reserve(total_features);
    arrays.label.reserve(total_features);
    arrays.spectrum_q.reserve(total_features);
    arrays.peptide_q.reserve(total_features);
    arrays.protein_q.reserve(total_features);
    arrays.posterior_error.reserve(total_features);
    arrays.frag_charges.reserve(total_features);
    arrays.frag_kinds.reserve(total_features);
    arrays.frag_fragment_ordinals.reserve(total_features);
    arrays.frag_intensities.reserve(total_features);
    arrays.frag_mz_calculated.reserve(total_features);
    arrays.frag_mz_experimental.reserve(total_features);


    // Iterate through the features and extract data
    for feat in all_features.into_iter() {
        // Fetch the associated Peptide object from the database once
        let pep = peptide_by_idx(&*db_arc, &feat.peptide_idx);

        // --- Basic Identifiers ---
        arrays.file_id.push(feat.file_id);
        arrays.spec_id.push(feat.spec_id);
        arrays.psm_id.push(feat.psm_id);
        arrays.rank.push(feat.rank);

        // --- Peptide Info ---
        if let Some(p) = pep {
            arrays.peptide_len.push(p.inner.sequence.len());
            arrays.missed_cleavages.push(feat.missed_cleavages);

            // Sequence
            arrays.sequence.push(String::from_utf8_lossy(&p.inner.sequence).into_owned());

            // Modifications: Assuming the underlying Peptide struct has a method
            // equivalent to PyPeptide::mod_array() (e.g., `p.mod_array()`)
            // that returns the N-term, residue, and C-term mods as Vec<f32>.
            arrays.modifications.push(p.mod_array());
        } else {
            // Handle features without an associated peptide
            arrays.peptide_len.push(0);
            arrays.missed_cleavages.push(0);
            arrays.sequence.push(String::new());
            arrays.modifications.push(vec![]);
        }

        // --- Mass and error metrics ---
        arrays.expmass.push(feat.expmass);
        arrays.calcmass.push(feat.calcmass);
        arrays.charge.push(feat.charge);
        arrays.delta_mass.push(feat.delta_mass);
        arrays.isotope_error.push(feat.isotope_error);
        arrays.average_ppm.push(feat.average_ppm);

        // --- Score metrics ---
        arrays.hyperscore.push(feat.hyperscore);
        arrays.delta_next.push(feat.delta_next);
        arrays.delta_best.push(feat.delta_best);
        arrays.poisson.push(feat.poisson);
        arrays.discriminant_score.push(feat.discriminant_score);

        // --- Spectrum matching metrics ---
        arrays.matched_peaks.push(feat.matched_peaks);
        arrays.matched_intensity_pct.push(feat.matched_intensity_pct);
        arrays.longest_b.push(feat.longest_b);
        arrays.longest_y.push(feat.longest_y);
        arrays.longest_y_pct.push(feat.longest_y_pct);
        arrays.scored_candidates.push(feat.scored_candidates);
        arrays.ms2_intensity.push(feat.ms2_intensity);

        // --- RT/IMS metrics ---
        arrays.rt.push(feat.rt);
        arrays.aligned_rt.push(feat.aligned_rt);
        arrays.predicted_rt.push(feat.predicted_rt);
        arrays.delta_rt_model.push(feat.delta_rt_model);
        arrays.ims.push(feat.ims);
        arrays.predicted_ims.push(feat.predicted_ims);
        arrays.delta_ims_model.push(feat.delta_ims_model);

        // --- Q-values and label ---
        arrays.label.push(feat.label);
        arrays.spectrum_q.push(feat.spectrum_q);
        arrays.peptide_q.push(feat.peptide_q);
        arrays.protein_q.push(feat.protein_q);
        arrays.posterior_error.push(feat.posterior_error);

        //--- Fragments ---
        if let Some(ref frags) = feat.fragments {
            // charges: Vec<i32>
            arrays.frag_charges.push(frags.charges.clone());

            // kinds: convert Kind enum -> String, same as PyFragments::kinds()
            arrays.frag_kinds.push(
                frags.kinds.iter().map(|k| format!("{:?}", k)).collect()
            );

            arrays.frag_fragment_ordinals.push(frags.fragment_ordinals.clone());
            arrays.frag_intensities.push(frags.intensities.clone());
            arrays.frag_mz_calculated.push(frags.mz_calculated.clone());
            arrays.frag_mz_experimental.push(frags.mz_experimental.clone());
        } else {
            // No fragments: store empty lists so the column lengths still match
            arrays.frag_charges.push(Vec::new());
            arrays.frag_kinds.push(Vec::new());
            arrays.frag_fragment_ordinals.push(Vec::new());
            arrays.frag_intensities.push(Vec::new());
            arrays.frag_mz_calculated.push(Vec::new());
            arrays.frag_mz_experimental.push(Vec::new());
        }
    }

    arrays
}

// --- 2. IMPLEMENT METHODS FOR WRAPPER STRUCTS ---

// Now that the struct is defined, you can write the methods.
#[pymethods]
impl PyIndexedDatabase {
    // ... all your static methods, getters, and other methods go here ...

    #[getter]
    pub fn peptides(&self) -> Vec<PyPeptide> {
        self.inner
            .peptides
            .iter()
            .map(|p| PyPeptide { inner: p.clone() })
            .collect()
    }

    #[staticmethod]
    #[pyo3(signature = (peptides, bucket_size, ion_kinds, min_ion_index, generate_decoys, decoy_tag, peptide_min_mass, peptide_max_mass))]
    pub fn from_peptides_and_config(
        peptides: Vec<PyPeptide>,
        bucket_size: usize,
        ion_kinds: Vec<PyKind>,
        min_ion_index: usize,
        generate_decoys: bool,
        decoy_tag: String,
        // Assuming PyModificationSpecificity is defined and correct
        //potential_mods: Vec<(PyModificationSpecificity, f32)>,
        peptide_min_mass: f32,
        peptide_max_mass: f32,
    ) -> PyResult<Self> {
        // 1. Convert PyO3 wrappers to native Rust types
        let native_peptides = peptides.into_iter().map(|p| p.inner).collect();
        let native_ion_kinds = ion_kinds.into_iter().map(|k| k.inner).collect();
        //let native_potential_mods = potential_mods
        //    .into_iter()
        //    .map(|(k, v)| (k.inner, v))
        //    .collect();

        // 2. Create the IndexingConfig struct
        let config = IndexingConfig {
            bucket_size,
            ion_kinds: native_ion_kinds,
            min_ion_index,
            generate_decoys,
            decoy_tag,
            //potential_mods: native_potential_mods,
            peptide_min_mass,
            peptide_max_mass,
        };

        // 3. Call the core Rust logic, handling panics
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            build_indexed_database(config, native_peptides)
        }));

        match result {
            Ok(db) => Ok(PyIndexedDatabase { inner: Arc::new(db) }), // ✅ wrap once here
            Err(_) => Err(PyRuntimeError::new_err(
                "Rust panic occurred during indexed database generation.",
            )),
        }
    }

    pub fn sequence_for(&self, pep_index: u32) -> Option<String> {
        self.inner.peptides
            .get(pep_index as usize)
            .map(|p| String::from_utf8_lossy(&p.sequence).into_owned())
    }

    /// Debug helper: return a summary of all indexed fragments
    pub fn debug_fragment_summary(&self) -> PyResult<Vec<(f32, usize)>> {
        let mut summary = Vec::new();

        for frag in &self.inner.fragments {
            let mass = frag.fragment_mz;
            let pep_idx = frag.peptide_index.0 as usize;
            summary.push((mass, pep_idx));
        }

        Ok(summary)
    }

    /// Quick count of total fragments
    pub fn fragment_count(&self) -> usize {
        self.inner.fragments.len()
    }
}

fn peptide_by_idx(db: &IndexedDatabase, idx: &PeptideIx) -> Option<PyPeptide> {
    db.peptides
        .get(idx.0 as usize)
        .map(|p| PyPeptide { inner: p.clone() })
}

#[pymethods]
impl PyScorer {
    /// Initializes the Scorer configuration. Defaults are set for a robust DIA/Chimeric search.
    #[new]
    #[pyo3(signature = (
        precursor_tol,
        fragment_tol,
        wide_window=true,
        chimera=true,
        report_psms=5,
        min_isotope_err=-1,
        max_isotope_err=3,
        min_precursor_charge=2,
        max_precursor_charge=4,
        min_matched_peaks=3,
        // score_type=PyScoreType::SageHyperScore, REMOVED and set as default
        annotate_matches=false,
        max_fragment_charge=2,
    ))]
    pub fn new(
        precursor_tol: PyTolerance,
        fragment_tol: PyTolerance,
        wide_window: bool,
        chimera: bool,
        report_psms: usize,
        min_isotope_err: i8,
        max_isotope_err: i8,
        min_precursor_charge: u8,
        max_precursor_charge: u8,
        min_matched_peaks: u16,
        // score_type: PyScoreType, REMOVED and set SageHyperScore as default
        annotate_matches: bool,
        max_fragment_charge: Option<u8>
    ) -> Self {
        PyScorer {
            config: ScorerConfig {
                precursor_tol: precursor_tol.inner,
                fragment_tol: fragment_tol.inner,
                wide_window,
                chimera,
                report_psms,
                min_matched_peaks,
                min_isotope_err,
                max_isotope_err,
                min_precursor_charge,
                max_precursor_charge,
                score_type: ScoreType::SageHyperScore,

                // Defaults for native Scorer fields
                override_precursor_charge: false, // May need to change these defaults TODO
                annotate_matches,
                max_fragment_charge,
            }
        }
    }

    pub fn score_many_spectra(
        &self,
        py: pyo3::Python<'_>,
        db: &PyIndexedDatabase,
        spectra: Vec<Py<PyProcessedSpectrum>>,   // <— faster extraction
    ) -> PyResult<PyFeatureArrays> {
        //let guard = ProfilerGuard::new(100).ok(); // 100 Hz sampler

        // 1) Snapshot Rust-owned data while holding the GIL
        let mut rust_specs: Vec<Arc<ProcessedSpectrum<Peak>>> = Vec::with_capacity(spectra.len());
        for handle in &spectra {
            let s = handle.borrow(py);
            rust_specs.push(Arc::clone(&s.inner)); // O(1) Arc clone
        }
        let db_arc = Arc::clone(&db.inner);        // O(1) Arc clone
        let cfg = self.config.clone();             // cheap struct clone

        // 2) Build a custom Rayon thread pool if requested
        let mut num_threads = None; //TODO this should be a parameter
        let results = py.allow_threads(|| {
            if let Some(n_threads) = num_threads {
                ThreadPoolBuilder::new()
                    .num_threads(n_threads)
                    .build()
                    .expect("Failed to build Rayon thread pool")
                    .install(|| {
                        rust_specs
                            .par_iter()
                            .map(|spec| crate::scoring_logic::run_scoring(&cfg, &db_arc, spec))
                            .collect::<Vec<_>>()
                    })
            } else {
                rust_specs
                    .par_iter()
                    .map(|spec| crate::scoring_logic::run_scoring(&cfg, &db_arc, spec))
                    .collect::<Vec<_>>()
            }
        });

        // The first flat_map flattens the results from Vec<Vec<Feature>> to Vec<Feature>
        let all_features: Vec<Feature> = results
            .into_iter()
            .flat_map(|v| v.into_iter())
            .collect();

        // Perform the final data extraction and formatting, all in fast Rust code
        let feature_arrays = py.allow_threads(|| {
            features_to_arrays(all_features, db_arc)
        });

        Ok(feature_arrays)
    }

    /// Executes the scoring logic against the built database and a single spectrum.
    pub fn score_spectra(
        &self,
        py: pyo3::Python<'_>,
        db: &PyIndexedDatabase,
        spectrum: &PyProcessedSpectrum
    ) -> PyResult<Vec<PyFeature>> {

        // 1. Delegate execution to the external logic function
        // run native Rust without the GIL
        let native_features = py.allow_threads(|| {
            run_scoring(&self.config, &db.inner, &spectrum.inner)
        });

        println!("Scoring returned {} features", native_features.len());

        // 2. Convert the native Rust output back to PyO3 wrappers
        let db_arc = Arc::clone(&db.inner);
        Ok(native_features.into_iter().map(|f| {
            // If db is Arc<IndexedDatabase>, pass &*db_arc to get &IndexedDatabase
            let pep = peptide_by_idx(&*db_arc, &f.peptide_idx);
            PyFeature { inner: f, peptide: pep }
        }).collect())
    }
}

#[pymethods]
impl PyPeptide {
    // --- CONSTRUCTOR: Allows creation from Python (e.g., PyPeptide("ABC", 1200.5)) ---
    // The signature attribute defines the Python API, allowing optional/default arguments.
    #[new]
    #[pyo3(signature = (sequence, monoisotopic, modifications, nterm=None, cterm=None))]
    pub fn new(
        sequence: String,
        monoisotopic: f32,
        modifications: Option<Vec<f32>>,
        nterm: Option<f32>,
        cterm: Option<f32>,
    ) -> Self {
        // WARNING: This is a minimal constructor. In a real project,
        // you'd take all relevant fields (proteins, mods, etc.) as arguments.
        let seq_len = sequence.len();

        // 1. Create the Box<[u8]> from the string
        let boxed_sequence: Box<[u8]> = sequence.into_bytes().into_boxed_slice();

        // 2. Wrap the Boxed sequence in an Arc<[u8]> to satisfy the Peptide struct
        //    We convert the Box<[u8]> to an Arc<[u8]>
        let arc_sequence: Arc<[u8]> = boxed_sequence.into();

        PyPeptide {
            inner: Peptide {
                // Convert Python string to the native Rust Box<[u8]> sequence type
                sequence: arc_sequence,
                monoisotopic,
                proteins: vec![Arc::from("USER_PROVIDED_PEPPYSAGE".to_string())],
                decoy: false,
                // Mods vector must match sequence length
                modifications: modifications.unwrap_or(vec![0.0; seq_len]),
                nterm: nterm,
                cterm: cterm,
                missed_cleavages: 0,
                semi_enzymatic: false,
                position: Position::default(),
            },
        }
    }

    #[getter]
    pub fn sequence(&self) -> PyResult<String> {
        // Convert the native Box<[u8]> sequence back to a Python string
        String::from_utf8(self.inner.sequence.to_vec())
            .map_err(|e| pyo3::exceptions::PyUnicodeError::new_err(e.to_string()))
    }

    #[getter]
    pub fn modifications(&self) -> PyResult<Vec<f32>> {
        Ok(self.mod_array())
    }

    /// Return combined mod array in the convention: [N-term, per-res..., C-term]
    #[getter]
    pub fn mod_array(&self) -> Vec<f32> {
        let n = self.inner.nterm.unwrap_or(0.0);
        let c = self.inner.cterm.unwrap_or(0.0);
        let mut out = Vec::with_capacity(self.inner.modifications.len() + 2);
        out.push(n);
        out.extend_from_slice(&self.inner.modifications);
        out.push(c);
        out
    }

    /// Accept a combined mod array [N-term, per-res..., C-term] and split it
    #[setter]
    pub fn set_mod_array(&mut self, mods: Vec<f32>) -> PyResult<()> {
        let seq_len = self.inner.modifications.len();
        if mods.len() != seq_len + 2 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "mod_array length {} must equal sequence length + 2 (= {})",
                mods.len(),
                seq_len + 2
            )));
        }
        self.inner.nterm = Some(mods[0]);
        self.inner.modifications.copy_from_slice(&mods[1..(1 + seq_len)]);
        self.inner.cterm = Some(mods[seq_len + 1]);
        Ok(())
    }
}


#[pymethods]
impl PyKind {
    // Implement a simple static constructor for the B-ion kind
    #[staticmethod]
    pub fn B() -> Self {
        Self { inner: Kind::B }
    }

    // Implement a simple static constructor for the Y-ion kind
    #[staticmethod]
    pub fn Y() -> Self {
        Self { inner: Kind::Y }
    }
}

// And your other impl blocks (impl PyPeptide, impl PyKind, etc.)
#[pymethods]
impl PyTolerance {
    #[staticmethod]
    pub fn Da(minus: f32, plus: f32) -> Self {
        Self { inner: Tolerance::Da(minus, plus) }
    }

    #[staticmethod]
    pub fn Ppm(minus: f32, plus: f32) -> Self {
        Self { inner: Tolerance::Ppm(minus, plus) }
    }
}

#[pymethods]
impl PyFragments {
    #[getter]
    pub fn charges(&self) -> Vec<i32> {
        self.inner.charges.clone()
    }

    #[getter]
    pub fn kinds(&self) -> Vec<String> {
        // Convert `Kind` enum variants to strings for readability in Python
        self.inner.kinds.iter().map(|k| format!("{:?}", k)).collect()
    }

    #[getter]
    pub fn fragment_ordinals(&self) -> Vec<i32> {
        self.inner.fragment_ordinals.clone()
    }

    #[getter]
    pub fn intensities(&self) -> Vec<f32> {
        self.inner.intensities.clone()
    }

    #[getter]
    pub fn mz_calculated(&self) -> Vec<f32> {
        self.inner.mz_calculated.clone()
    }

    #[getter]
    pub fn mz_experimental(&self) -> Vec<f32> {
        self.inner.mz_experimental.clone()
    }

    /// Optional helper: return as Python dict for easier pandas conversion
    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new_bound(py);
        dict.set_item("charges", self.charges())?;
        dict.set_item("kinds", self.kinds())?;
        dict.set_item("fragment_ordinals", self.fragment_ordinals())?;
        dict.set_item("intensities", self.intensities())?;
        dict.set_item("mz_calculated", self.mz_calculated())?;
        dict.set_item("mz_experimental", self.mz_experimental())?;
        Ok(dict)
    }
}

#[pymethods]
impl PyFeature {
    #[getter]
    pub fn peptide(&self) -> Option<PyPeptide> {
        self.peptide.clone()
    }

    #[setter]
    pub fn set_peptide(&mut self, peptide: Option<PyPeptide>) {
        self.peptide = peptide;
    }

    #[getter]
    pub fn sequence(&self) -> Option<String> {
        self.peptide.as_ref()
            .map(|p| String::from_utf8_lossy(&p.inner.sequence).into_owned())
    }

    #[getter]
    pub fn modifications(&self) -> Option<Vec<f32>> {
        self.peptide.as_ref().map(|p| p.mod_array())
    }

    #[getter]
    pub fn hyperscore(&self) -> f64 {
        self.inner.hyperscore
    }

    #[getter]
    pub fn delta_mass(&self) -> f32 {
        self.inner.delta_mass
    }

    #[getter]
    pub fn rank(&self) -> u32 {
        self.inner.rank
    }

    #[getter]
    pub fn matched_peaks(&self) -> u32 {
        self.inner.matched_peaks
    }

    #[getter]
    pub fn peptide_idx(&self) -> u32 {
        self.inner.peptide_idx.0
    }

    #[getter]
    pub fn psm_id(&self) -> usize {
        self.inner.psm_id
    }

    #[getter]
    pub fn peptide_len(&self) -> usize {
        self.inner.peptide_len
    }

    #[getter]
    pub fn spec_id(&self) -> &str {
        &self.inner.spec_id
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[getter]
    pub fn label(&self) -> i32 {
        self.inner.label
    }

    #[getter]
    pub fn expmass(&self) -> f32 {
        self.inner.expmass
    }

    #[getter]
    pub fn calcmass(&self) -> f32 {
        self.inner.calcmass
    }

    #[getter]
    pub fn charge(&self) -> u8 {
        self.inner.charge
    }

    #[getter]
    pub fn rt(&self) -> f32 {
        self.inner.rt
    }

    #[getter]
    pub fn aligned_rt(&self) -> f32 {
        self.inner.aligned_rt
    }

    #[getter]
    pub fn predicted_rt(&self) -> f32 {
        self.inner.predicted_rt
    }

    #[getter]
    pub fn delta_rt_model(&self) -> f32 {
        self.inner.delta_rt_model
    }

    #[getter]
    pub fn ims(&self) -> f32 {
        self.inner.ims
    }

    #[getter]
    pub fn predicted_ims(&self) -> f32 {
        self.inner.predicted_ims
    }

    #[getter]
    pub fn delta_ims_model(&self) -> f32 {
        self.inner.delta_ims_model
    }

    #[getter]
    pub fn isotope_error(&self) -> f32 {
        self.inner.isotope_error
    }

    #[getter]
    pub fn average_ppm(&self) -> f32 {
        self.inner.average_ppm
    }

    #[getter]
    pub fn delta_next(&self) -> f64 {
        self.inner.delta_next
    }

    #[getter]
    pub fn delta_best(&self) -> f64 {
        self.inner.delta_best
    }

    #[getter]
    pub fn longest_b(&self) -> u32 {
        self.inner.longest_b
    }

    #[getter]
    pub fn longest_y(&self) -> u32 {
        self.inner.longest_y
    }

    #[getter]
    pub fn longest_y_pct(&self) -> f32 {
        self.inner.longest_y_pct
    }

    #[getter]
    pub fn missed_cleavages(&self) -> u8 {
        self.inner.missed_cleavages
    }

    #[getter]
    pub fn matched_intensity_pct(&self) -> f32 {
        self.inner.matched_intensity_pct
    }

    #[getter]
    pub fn scored_candidates(&self) -> u32 {
        self.inner.scored_candidates
    }

    #[getter]
    pub fn poisson(&self) -> f64 {
        self.inner.poisson
    }

    #[getter]
    pub fn discriminant_score(&self) -> f32 {
        self.inner.discriminant_score
    }

    #[getter]
    pub fn posterior_error(&self) -> f32 {
        self.inner.posterior_error
    }

    #[getter]
    pub fn spectrum_q(&self) -> f32 {
        self.inner.spectrum_q
    }

    #[getter]
    pub fn peptide_q(&self) -> f32 {
        self.inner.peptide_q
    }

    #[getter]
    pub fn protein_q(&self) -> f32 {
        self.inner.protein_q
    }

    #[getter]
    pub fn ms2_intensity(&self) -> f32 {
        self.inner.ms2_intensity
    }

    #[getter]
    pub fn fragments(&self) -> Option<PyFragments> {
        self.inner.fragments.clone().map(|f| PyFragments { inner: f })
    }

    /// Helpful repr for quick printing
    pub fn __repr__(&self) -> PyResult<String> {
        let (seq, mods) = if let Some(p) = &self.peptide {
            (
                String::from_utf8_lossy(&p.inner.sequence).into_owned(),
                format!("{:?}", p.mod_array()),
            )
        } else {
            ("<None>".to_string(), "<None>".to_string())
        };

        Ok(format!(
            "Feature(spec={}, peptide_idx={}, rank={}, seq={}, mods={}, hyperscore={:.3}, \
                delta_mass={:.4}, \n matched_peaks={}, isotope error={}, average_ppm={:.3}, poisson={:.4})",
            self.inner.spec_id,
            self.inner.peptide_idx.0,
            self.inner.rank,
            seq,
            mods,
            self.inner.hyperscore,
            self.inner.delta_mass,
            self.inner.matched_peaks,
            self.inner.isotope_error,
            self.inner.average_ppm,
            self.inner.poisson
        ))
    }

    pub fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        let d = PyDict::new_bound(py);
        d.set_item("file_id", self.inner.file_id)?;
        d.set_item("spec_id", self.inner.spec_id.clone())?;
        d.set_item("psm_id", self.inner.psm_id)?;
        d.set_item("rank", self.inner.rank)?;
        if let Some(p) = &self.peptide {
            d.set_item("sequence", String::from_utf8_lossy(&p.inner.sequence).into_owned())?;
            d.set_item("modifications", p.mod_array())?;
        } else {
            d.set_item("sequence", py.None())?;
            d.set_item("modifications", py.None())?;
        }
        d.set_item("label", self.inner.label)?;
        d.set_item("hyperscore", self.inner.hyperscore)?;
        d.set_item("delta_mass", self.inner.delta_mass)?;
        d.set_item("matched_peaks", self.inner.matched_peaks)?;
        d.set_item("peptide_len", self.inner.peptide_len)?;
        d.set_item("expmass", self.inner.expmass)?;
        d.set_item("calcmass", self.inner.calcmass)?;
        d.set_item("charge", self.inner.charge)?;
        d.set_item("rt", self.inner.rt)?;
        d.set_item("aligned_rt", self.inner.aligned_rt)?;
        d.set_item("predicted_rt", self.inner.predicted_rt)?;
        d.set_item("delta_rt_model", self.inner.delta_rt_model)?;
        d.set_item("ims", self.inner.ims)?;
        d.set_item("predicted_ims", self.inner.predicted_ims)?;
        d.set_item("delta_ims_model", self.inner.delta_ims_model)?;
        d.set_item("isotope_error", self.inner.isotope_error)?;
        d.set_item("average_ppm", self.inner.average_ppm)?;
        d.set_item("delta_next", self.inner.delta_next)?;
        d.set_item("delta_best", self.inner.delta_best)?;
        d.set_item("longest_b", self.inner.longest_b)?;
        d.set_item("longest_y", self.inner.longest_y)?;
        d.set_item("longest_y_pct", self.inner.longest_y_pct)?;
        d.set_item("missed_cleavages", self.inner.missed_cleavages)?;
        d.set_item("matched_intensity_pct", self.inner.matched_intensity_pct)?;
        d.set_item("scored_candidates", self.inner.scored_candidates)?;
        d.set_item("poisson", self.inner.poisson)?;
        d.set_item("discriminant_score", self.inner.discriminant_score)?;
        d.set_item("posterior_error", self.inner.posterior_error)?;
        d.set_item("spectrum_q", self.inner.spectrum_q)?;
        d.set_item("peptide_q", self.inner.peptide_q)?;
        d.set_item("protein_q", self.inner.protein_q)?;
        d.set_item("ms2_intensity", self.inner.ms2_intensity)?;

        // Optionally include fragment info, if present
        if let Some(f) = &self.inner.fragments {
            d.set_item("fragments", PyFragments { inner: f.clone() }.to_dict(py)?)?;
        } else {
            d.set_item("fragments", py.None())?;
        }

        Ok(d.into())
    }
}

#[pymethods]
impl PyFeatureArrays {
    /// Extends this PyFeatureArrays object with the contents of another,
    /// performing the vector concatenation efficiently in Rust.
    pub fn extend(&mut self, other: PyFeatureArrays) {
        self.file_id.extend(other.file_id);
        self.spec_id.extend(other.spec_id);
        self.psm_id.extend(other.psm_id);
        self.rank.extend(other.rank);

        // Peptide info
        self.sequence.extend(other.sequence);
        self.modifications.extend(other.modifications);
        self.peptide_len.extend(other.peptide_len);
        self.missed_cleavages.extend(other.missed_cleavages);

        // Mass and error metrics
        self.expmass.extend(other.expmass);
        self.calcmass.extend(other.calcmass);
        self.charge.extend(other.charge);
        self.delta_mass.extend(other.delta_mass);
        self.isotope_error.extend(other.isotope_error);
        self.average_ppm.extend(other.average_ppm);

        // Score metrics
        self.hyperscore.extend(other.hyperscore);
        self.delta_next.extend(other.delta_next);
        self.delta_best.extend(other.delta_best);
        self.poisson.extend(other.poisson);
        self.discriminant_score.extend(other.discriminant_score);

        // Spectrum matching metrics
        self.matched_peaks.extend(other.matched_peaks);
        self.matched_intensity_pct.extend(other.matched_intensity_pct);
        self.longest_b.extend(other.longest_b);
        self.longest_y.extend(other.longest_y);
        self.longest_y_pct.extend(other.longest_y_pct);
        self.scored_candidates.extend(other.scored_candidates);
        self.ms2_intensity.extend(other.ms2_intensity);

        // RT/IMS metrics
        self.rt.extend(other.rt);
        self.aligned_rt.extend(other.aligned_rt);
        self.predicted_rt.extend(other.predicted_rt);
        self.delta_rt_model.extend(other.delta_rt_model);
        self.ims.extend(other.ims);
        self.predicted_ims.extend(other.predicted_ims);
        self.delta_ims_model.extend(other.delta_ims_model);

        // Q-values and label
        self.label.extend(other.label);
        self.spectrum_q.extend(other.spectrum_q);
        self.peptide_q.extend(other.peptide_q);
        self.protein_q.extend(other.protein_q);
        self.posterior_error.extend(other.posterior_error);

        // Fragment columns
        self.frag_charges.extend(other.frag_charges);
        self.frag_kinds.extend(other.frag_kinds);
        self.frag_fragment_ordinals.extend(other.frag_fragment_ordinals);
        self.frag_intensities.extend(other.frag_intensities);
        self.frag_mz_calculated.extend(other.frag_mz_calculated);
        self.frag_mz_experimental.extend(other.frag_mz_experimental);
    }

    /// Provides a debugging representation for the Python object.
    pub fn __repr__(&self) -> PyResult<String> {
        let n = self.file_id.len();
        let sample_size = n.min(3);

        let seq_sample: Vec<String> = self.sequence.iter().take(sample_size).cloned().collect();
        let score_sample: Vec<f64> = self.hyperscore.iter().take(sample_size).copied().collect();
        let charge_sample: Vec<u8> = self.charge.iter().take(sample_size).copied().collect();

        Ok(format!(
            "PyFeatureArrays(\n\
             \tFeatures: {n},\n\
             \tSample Sequence: {:?},\n\
             \tSample Hyperscore: {:?},\n\
             \tSample Charge: {:?}\n\
             )",
            seq_sample,
            score_sample,
            charge_sample,
        ))
    }

    /// Returns a list of all column names (field names) available in this struct.
    pub fn get_column_names(&self) -> Vec<String> {
        let debug_string = format!("{:?}", self);
        // This is a common trick: the Debug output for a struct lists all field names.
        // We parse this string to get the names.

        let struct_name = "PyFeatureArrays {";
        if let Some(start) = debug_string.find(struct_name) {
            let fields_string = &debug_string[start + struct_name.len()..];
            let mut names = Vec::new();

            for field in fields_string.split(", ") {
                if let Some(name_end) = field.find(':') {
                    let name = field[..name_end].trim().to_string();
                    if !name.is_empty() {
                        names.push(name);
                    }
                }
            }
            // Remove the trailing '}' from the last field name if present
            if let Some(last) = names.last_mut() {
                if last.ends_with('}') {
                    last.pop();
                }
            }
            return names;
        }
        Vec::new()
    }
}

#[pymethods]
impl PyProcessedSpectrum {
    #[new]
    #[pyo3(signature = (id, file_id, scan_start_time, mz_array, intensity_array, precursors, total_ion_current))]
    pub fn new<'py>(
        py: Python<'py>,                        // get a GIL handle in the constructor
        id: String,
        file_id: usize,
        scan_start_time: f32,
        // accept *any* buffer-exporting object (NumPy array, memoryview, array('f'), etc.)
        mz_array: Bound<'py, PyAny>,
        intensity_array: Bound<'py, PyAny>,
        precursors: Vec<PyPrecursor>,
        total_ion_current: f32,
    ) -> PyResult<Self> {
        // 1) Acquire typed buffers (float32) using the *bound* API
        let buf_mz:  PyBuffer<f32> = PyBuffer::get_bound(&mz_array)?;
        let buf_int: PyBuffer<f32> = PyBuffer::get_bound(&intensity_array)?;

        // 2) Validate shape/layout
        if !buf_mz.is_c_contiguous() || !buf_int.is_c_contiguous() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "mz_array and intensity_array must be C-contiguous float32 buffers",
            ));
        }
        if buf_mz.item_count() != buf_int.item_count() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "mz_array and intensity_array must be the same length",
            ));
        }

        // 3) Borrow as slices (zero Python-float boxing), then copy once to Vecs
        //    SAFETY: dtype verified by PyBuffer<f32> and contiguity checked above
        let mz_cells  = unsafe { buf_mz.as_slice(py) }
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("Failed to get mz_array slice"))?;
        let int_cells = unsafe { buf_int.as_slice(py) }
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("Failed to get intensity_array slice"))?;

        // 4) Pack into Peaks (owned, mutable on Rust side)
        let peaks: Vec<Peak> = mz_cells.iter()
            .zip(int_cells.iter())
            .map(|(m, i)| Peak { mass: m.get(), intensity: i.get() })
            .collect();

        Ok(PyProcessedSpectrum {
            inner: Arc::new(ProcessedSpectrum {
                level: 2,
                id,
                file_id,
                scan_start_time,
                ion_injection_time: 0.0,
                precursors: precursors.into_iter().map(|p| p.inner).collect(),
                peaks,
                total_ion_current,
            }),
        })
    }

    #[getter]
    pub fn id(&self) -> PyResult<String> {
        Ok(self.inner.id.clone())
    }

    #[getter]
    pub fn peaks(&self) -> PyResult<Vec<(f32, f32)>> {
        // return Vec of (mz, intensity)
        Ok(self
            .inner
            .peaks
            .iter()
            .map(|p| (p.mass, p.intensity))
            .collect())
    }

    #[getter]
    pub fn precursors(&self) -> PyResult<Vec<PyPrecursor>> {
        Ok(self
            .inner
            .precursors
            .iter()
            .cloned()
            .map(|p| PyPrecursor { inner: p })
            .collect())
    }
}

#[pymethods]
impl PyPrecursor {
    #[new]
    #[pyo3(signature = (mz, charge=None, intensity=None, isolation_window=None, inverse_ion_mobility=None))]
    pub fn new(
        mz: f32,
        charge: Option<u8>,
        intensity: Option<f32>,
        isolation_window: Option<PyTolerance>,
        inverse_ion_mobility: Option<f32>,
    ) -> PyResult<Self> {
        Ok(PyPrecursor {
            inner: Precursor {
                mz,
                intensity,
                charge,
                spectrum_ref: None,
                isolation_window: isolation_window.map(|t| t.inner),
                inverse_ion_mobility,
                },
            }
        )}
    #[getter]
    pub fn mz(&self) -> PyResult<f32> {
        Ok(self.inner.mz)
    }
    #[getter]
    pub fn charge(&self) -> PyResult<Option<u8>> {
        Ok(self.inner.charge)
    }
    #[getter]
    pub fn intensity(&self) -> PyResult<Option<f32>> {
        Ok(self.inner.intensity)
    }
    #[getter]
    pub fn isolation_window(&self) -> PyResult<Option<PyTolerance>> {
        Ok(self.inner.isolation_window.map(|t| PyTolerance { inner: t }))
    }
}

// --- 3. PYMODULE EXPORT ---

// The #[pymodule] function goes last.
#[cfg(test)]
mod tests {
    use super::*;
    use sage_core::ion_series::Kind;

    // PyO3 utility
    use pyo3::Python;

    /// Helper to create a basic PyPeptide for input.
    fn create_test_peptide(sequence: &str, monoisotopic: f32, modifications: Vec<f32>) -> PyPeptide {
        // Use the PyPeptide constructor we just implemented
        PyPeptide::new(
            sequence.to_string(),
            monoisotopic,
            Some(modifications),
            Some(0.0), // nterm
            Some(0.0)  // cterm
        )
    }

    #[test]
    fn test_build_from_peptides_fasta_free() {
        // 0. Initialize PyO3 environment
        pyo3::prepare_freethreaded_python();

        Python::with_gil(|_py| {
            // 1. Setup Input Peptides
            let sequence_a = "ABCDEFK";
            let sequence_b = "XYZPQKR";
            let peptides = vec![
                // Note: Monoisotopic masses must be somewhat realistic
                create_test_peptide(sequence_a, 700.0, vec![0.0; sequence_a.len()]),
                create_test_peptide(sequence_b, 950.0, vec![0.0; sequence_b.len()]),
            ];
            let num_peptides = peptides.len();

            // 2. Call the new PyO3 constructor method
            let py_db_result = PyIndexedDatabase::from_peptides_and_config(
                peptides,
                128, // bucket_size
                vec![PyKind { inner: Kind::B }, PyKind { inner: Kind::Y }], // ion_kinds
                2, // min_ion_index (skips b1, b2, y1, y2)
                true, // generate_decoys
                "rev_".to_string(), // decoy_tag
                500.0, // peptide_min_mass
                1500.0 // peptide_max_mass
            );

            // 3. Assertions
            assert!(py_db_result.is_ok(),
                    "Failed to build PyIndexedDatabase: {:?}", py_db_result.err());

            let py_db = py_db_result.unwrap();

            // a. Verify peptide count (targets + decoys)
            // Expect 4 total: A (target), A (decoy), B (target), B (decoy)
            assert_eq!(py_db.inner.peptides.len(), 4,
                       "Should have 2 targets and 2 decoys.");

            // b. Verify fragments were generated
            // A 7-residue peptide generates (7-1)*2 = 12 fragments. After min_ion_index=2 filter, roughly 8 per kind * 2 kinds = 16.
            // Total fragments should be around 60+.
            assert!(py_db.inner.fragments.len() > 50,
                    "Should have generated a significant number of fragments.");

            // c. Verify a getter works on the resulting objects
            let first_peptide = py_db.peptides().into_iter().next().unwrap();
            assert!(first_peptide.sequence().unwrap().contains(sequence_a) ||
                        first_peptide.sequence().unwrap().contains(sequence_b),
                    "First peptide sequence does not match input.");
        });
    }
}

// Build Test
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyfunction]
fn hello_py() -> PyResult<&'static str> {
    Ok("Hello from Rust + PyO3 🎉")
}

#[pymodule]
fn _peppy_sage(py: Python, m: Bound<'_, PyModule>) -> PyResult<()> {
    // Register ALL classes used for I/O and configuration
    m.add_class::<PyIndexedDatabase>()?;
    m.add_class::<PyPeptide>()?;
    m.add_class::<PyKind>()?;
    m.add_class::<PyTolerance>()?; // New
    m.add_class::<PyPrecursor>()?; // New
    m.add_class::<PyProcessedSpectrum>()?; // New
    m.add_class::<PyScorer>()?; // New
    m.add_class::<PyFeature>()?; // New

    m.add_function(wrap_pyfunction!(hello_py, m.clone())?)?;

    Ok(())
}