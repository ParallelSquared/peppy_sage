use sage_core::database::{IndexedDatabase, PeptideIx, Theoretical};
use sage_core::peptide::Peptide;
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::modification::ModificationSpecificity;
use sage_core::enzyme::Position;
use std::sync::Arc;
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use rayon::prelude::*;
use pyo3_polars::PyDataFrame;
use polars_core::frame::DataFrame;
use polars_core::prelude::*;
// Add imports for types used in the original logic (e.g., DashSet, FnvBuildHasher, etc. if you copied that part)

// 1. Define the minimal configuration struct
pub struct IndexingConfig {
    pub bucket_size: usize,
    pub ion_kinds: Vec<Kind>,
    pub min_ion_index: usize,
    pub generate_decoys: bool,
    pub decoy_tag: String,
    //pub potential_mods: Vec<(ModificationSpecificity, f32)>,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
}

// 2. Define the core indexing function (pasting the adapted logic)
pub fn build_indexed_database(config: IndexingConfig, targets: Vec<Peptide>) -> IndexedDatabase {
    // NOTE: This is where you paste the full, adapted logic from
    // the original Parameters::build_from_peptides method.

    let peptides = if config.generate_decoys {
        // gather target sequences for collision checks
        let target_seqs: DashSet<Arc<[u8]>, FnvBuildHasher> = DashSet::default();
        targets
            .par_iter()
            .filter(|p| !p.decoy)
            .for_each(|p| {
                target_seqs.insert(p.sequence.clone());
            });

        // build output = targets (+ valid reversed decoys)
        let out: Vec<Peptide> = targets
            .into_par_iter()
            .flat_map_iter(|p| {
                // always keep the target
                let mut v = Vec::with_capacity(2);
                v.push(p.clone());

                // reversed decoy; skip if reversed equals a real target sequence
                let rev = p.reverse();
                if !target_seqs.contains(&rev.sequence) {
                    v.push(rev);
                }
                v.into_iter()
            })
            .collect(); // CHANGE: collect the parallel stream directly
        out
    } else {
        targets
    };

    // --- Adapted Fragmentation Logic ---
    let mut fragments = peptides
        .par_iter()
        .enumerate()
        .flat_map(|(idx, peptide)| {
            // ... (mass filter logic remains the same) ...

            // Generate ions (using config fields)
            config.ion_kinds
                .iter()
                .flat_map(|kind| IonSeries::new(peptide, *kind).enumerate())
                .filter(|(ion_idx, ion)| {
                    // ... (filtering logic remains the same) ...
                    match ion.kind {
                        Kind::A | Kind::B | Kind::C => (ion_idx + 1) > config.min_ion_index,
                        Kind::X | Kind::Y | Kind::Z => {
                            peptide.sequence.len().saturating_sub(1) - ion_idx > config.min_ion_index
                        }
                    }
                })
                .map(move |(_, ion)| Theoretical {
                    peptide_index: PeptideIx(idx as u32),
                    fragment_mz: ion.monoisotopic_mass,
                })
                // FIX: Collect into a Vec, and immediately convert that Vec
                // into a parallel iterator (`IntoParallelIterator`) for Rayon.
                .collect::<Vec<Theoretical>>()
                .into_par_iter() // <--- THIS IS THE KEY CHANGE
        })
        // Final collect is now correct because flat_map returns a ParallelIterator
        .collect::<Vec<_>>();

    // --- Sorting and Bucketing Logic ---
    fragments.par_sort_unstable_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

    let min_value = fragments
        .par_chunks_mut(config.bucket_size)
        .map(|chunk| {
            let min = chunk[0].fragment_mz;
            chunk.par_sort_unstable_by(|a, b| a.peptide_index.cmp(&b.peptide_index));
            min
        })
        .collect::<Vec<_>>();

    // 3. Final Struct Return
    IndexedDatabase {
        peptides: peptides,
        fragments,
        min_value,
        bucket_size: config.bucket_size,
        ion_kinds: config.ion_kinds,
        generate_decoys: config.generate_decoys,
        potential_mods: Vec::new(), // inserting dummy values here since they won't be used
        decoy_tag: config.decoy_tag,
    }
}

/// Local sort+dedup that mirrors `Parameters::reorder_peptides`.
fn sort_and_dedup(peptides: &mut Vec<Peptide>) {
    peptides.par_sort_unstable_by(|a, b| {
        a.monoisotopic
            .total_cmp(&b.monoisotopic)
            .then_with(|| a.initial_sort(b))
    });

    peptides.dedup_by(|remove, keep| {
        if remove.monoisotopic == keep.monoisotopic
            && remove.sequence == keep.sequence
            && remove.modifications == keep.modifications
            && remove.nterm == keep.nterm
            && remove.cterm == keep.cterm
        {
            keep.proteins.extend(remove.proteins.iter().cloned());
            // if any source was target, the merged stays target
            keep.decoy &= remove.decoy;
            true
        } else {
            false
        }
    });

    peptides
        .par_iter_mut()
        .for_each(|p| p.proteins.sort_unstable());
}

pub fn build_indexed_database_from_library(
    config: IndexingConfig,
    df: DataFrame,
) -> IndexedDatabase {
    let n_rows = df.height();

    // Required columns
    let seq_col = df
        .column("sequence")
        .expect("library missing required column 'sequence'");
    let mods_col = df
        .column("modifications")
        .expect("library missing required column 'modifications'");
    let frag_mz_col = df
        .column("fragment_mz")
        .expect("library missing required column 'fragment_mz'");
    let frag_int_col = df
        .column("fragment_intensity")
        .expect("library missing required column 'fragment_intensity'");

    let mut peptides: Vec<Peptide> = Vec::with_capacity(n_rows);
    let mut fragments: Vec<Theoretical> = Vec::new();

    for idx in 0..n_rows {
        // -------------------- sequence --------------------
        let seq = seq_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'sequence' at row {}", idx))
            .to_string();

        let seq_bytes: Arc<[u8]> = seq.clone().into_bytes().into_boxed_slice().into();

        // -------------------- modifications --------------------
        let mods_val = mods_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'modifications' at row {}", idx));

        let mods_series = match mods_val {
            AnyValue::List(s) => s,
            other => panic!(
                "'modifications' must be List(Float32) per row, got {:?} at row {}",
                other, idx
            ),
        };

        let mut mods_vec: Vec<f32> = Vec::with_capacity(mods_series.len());

        for v in mods_series.iter() {
            match v {
                AnyValue::Float32(x) => mods_vec.push(x),
                AnyValue::Float64(x) => mods_vec.push(x as f32),
                AnyValue::Int32(x) => mods_vec.push(x as f32),
                AnyValue::Int64(x) => mods_vec.push(x as f32),
                AnyValue::Null => mods_vec.push(0.0), // or panic! if you prefer
                other => panic!(
                    "'modifications' inner values must be numeric, got {:?} at row {}",
                    other, idx
                ),
            }
        }

        // mods = [N-term, per-res..., C-term]
        if mods_vec.len() != seq.len() + 2 {
            panic!(
                "modifications length {} != sequence length + 2 ({}) at row {} (seq='{}')",
                mods_vec.len(),
                seq.len() + 2,
                idx,
                seq
            );
        }

        let nterm = mods_vec[0];
        let cterm = mods_vec[mods_vec.len() - 1];
        let per_res_mods = mods_vec[1..mods_vec.len() - 1].to_vec();

        // -------------------- build peptide --------------------
        let peptide = Peptide {
            sequence: seq_bytes.clone(),
            monoisotopic: 0.0, // you can fill from precursor mass later if you want
            proteins: vec![Arc::from("LIBRARY".to_string())],
            decoy: false,
            modifications: per_res_mods,
            nterm: Some(nterm),
            cterm: Some(cterm),
            missed_cleavages: 0,
            semi_enzymatic: false,
            position: Position::default(),
        };

        let pep_ix = PeptideIx(idx as u32);

        // -------------------- fragments --------------------
        let frag_mz_val = frag_mz_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'fragment_mz' at row {}", idx));
        let frag_int_val = frag_int_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'fragment_intensity' at row {}", idx));

        let frag_mz_series = match frag_mz_val {
            AnyValue::List(s) => s,
            other => panic!(
                "'fragment_mz' must be List(Float32) per row, got {:?} at row {}",
                other, idx
            ),
        };

        let frag_int_series = match frag_int_val {
            AnyValue::List(s) => s,
            other => panic!(
                "'fragment_intensity' must be List(Float32) per row, got {:?} at row {}",
                other, idx
            ),
        };

        let mut mz_vec: Vec<f32> = Vec::with_capacity(frag_mz_series.len());
        let mut int_vec: Vec<f32> = Vec::with_capacity(frag_int_series.len());

        for v in frag_mz_series.iter() {
            match v {
                AnyValue::Float32(x) => mz_vec.push(x),
                AnyValue::Float64(x) => mz_vec.push(x as f32),
                AnyValue::Int32(x) => mz_vec.push(x as f32),
                AnyValue::Int64(x) => mz_vec.push(x as f32),
                AnyValue::Null => mz_vec.push(0.0),
                other => panic!(
                    "'fragment_mz' inner values must be numeric, got {:?} at row {}",
                    other, idx
                ),
            }
        }

        for v in frag_int_series.iter() {
            match v {
                AnyValue::Float32(x) => int_vec.push(x),
                AnyValue::Float64(x) => int_vec.push(x as f32),
                AnyValue::Int32(x) => int_vec.push(x as f32),
                AnyValue::Int64(x) => int_vec.push(x as f32),
                AnyValue::Null => int_vec.push(0.0),
                other => panic!(
                    "'fragment_intensity' inner values must be numeric, got {:?} at row {}",
                    other, idx
                ),
            }
        }

        if mz_vec.len() != int_vec.len() {
            panic!(
                "fragment_mz len {} != fragment_intensity len {} at row {}",
                mz_vec.len(),
                int_vec.len(),
                idx
            );
        }

        // For now we only care about fragment_mz for indexing
        for mz in mz_vec {
            fragments.push(Theoretical {
                peptide_index: pep_ix,
                fragment_mz: mz,
            });
        }

        peptides.push(peptide);
    }

    // -------------------- decoys --------------------
    // For now, ignore config.generate_decoys as you suggested.
    // If you want to enforce that it's false:
    if config.generate_decoys {
        eprintln!(
            "Warning: generate_decoys is not yet supported in build_indexed_database_from_library; returning targets only."
        );
    }
    let peptides = peptides;

    // -------------------- sort + bucket fragments --------------------
    fragments.par_sort_unstable_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

    let min_value = if fragments.is_empty() {
        Vec::new()
    } else {
        fragments
            .par_chunks_mut(config.bucket_size)
            .map(|chunk| {
                let min = chunk[0].fragment_mz;
                chunk.par_sort_unstable_by(|a, b| a.peptide_index.cmp(&b.peptide_index));
                min
            })
            .collect::<Vec<_>>()
    };

    IndexedDatabase {
        peptides,
        fragments,
        min_value,
        bucket_size: config.bucket_size,
        ion_kinds: config.ion_kinds.clone(), // unused but kept for compat
        generate_decoys: config.generate_decoys,
        potential_mods: Vec::new(),
        decoy_tag: config.decoy_tag.clone(),
    }
}