use sage_core::database::{IndexedDatabase, PeptideIx, Theoretical};
use sage_core::peptide::Peptide;
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::modification::ModificationSpecificity;
use rayon::prelude::*;
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
pub fn build_indexed_database(config: IndexingConfig, target_decoys: Vec<Peptide>) -> IndexedDatabase {
    // NOTE: This is where you paste the full, adapted logic from
    // the original Parameters::build_from_peptides method.

    // TODO DEBUG LOGIC
    println!(
        "\n[DEBUG] Building IndexedDatabase with peptide mass window: {:.4} - {:.4}",
        config.peptide_min_mass, config.peptide_max_mass
    );

    for (i, pep) in target_decoys.iter().enumerate() {
        let mass = pep.monoisotopic;
        let passes = mass >= config.peptide_min_mass && mass <= config.peptide_max_mass;
        let seq = String::from_utf8_lossy(&pep.sequence);
        println!(
            "[DEBUG] Peptide {:>2}: {:<10}  mass={:.4}  within_range={} ",
            i, seq, mass, passes
        );
    }

    // --- Adapted Fragmentation Logic ---
    let mut fragments = target_decoys
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
        peptides: target_decoys,
        fragments,
        min_value,
        bucket_size: config.bucket_size,
        ion_kinds: config.ion_kinds,
        generate_decoys: config.generate_decoys,
        potential_mods: Vec::new(), // inserting dummy values here since they won't be used
        decoy_tag: config.decoy_tag,
    }
}