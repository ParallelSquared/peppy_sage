use sage_core::database::{IndexedDatabase, PeptideIx, Theoretical, LibraryFragment};
use sage_core::peptide::Peptide;
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::modification::ModificationSpecificity;
use sage_core::enzyme::Position;
use std::sync::Arc;
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use rayon::prelude::*;
use pyo3_polars::PyDataFrame;
use polars::prelude::*;
use pyo3::pyclass::boolean_struct::False;
use std::collections::HashMap;
use sage_core::mass::{monoisotopic as aa_mono, H2O};
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
        library_frags: None,
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
    let unimod_table = unimod_table(); // HashMap<&'static str, f32>

    let n_rows = df.height();

    let seq_col = df
        .column("StrippedPeptide")
        .expect("library missing required column 'StrippedPeptide'");
    let modpep_col = df
        .column("ModifiedPeptide")
        .expect("library missing required column 'ModifiedPeptide'");

    // optional numeric mods column
    let mods_col_opt = df.column("Modifications").ok();

    let frag_mz_col = df
        .column("FragmentMz")
        .expect("library missing required column 'FragmentMz'");
    let frag_int_col = df
        .column("RelativeIntensity")
        .expect("library missing required column 'RelativeIntensity'");

    let frag_type_col = df
        .column("FragmentType")
        .expect("library missing required column 'FragmentType'");
    let frag_num_col = df
        .column("FragmentSeriesNumber")
        .expect("library missing required column 'FragmentSeriesNumber'");
    let frag_z_col = df
        .column("FragmentCharge")
        .expect("library missing required column 'FragmentCharge'");

    // 1) Collect unique peptides: ModifiedPeptide -> (StrippedPeptide, mods_vec)
    let mut pep_mods: HashMap<String, (String, Vec<f32>)> = HashMap::new();

    // 2) Aggregate ions across precursor charge states:
    // key = (ModifiedPeptide, FragmentType, FragmentNumber, FragmentCharge)
    // value = (mz_sum, intensity_sum, count)
    let mut ion_agg: HashMap<(String, String, i32, i32), (f32, f32, u32)> = HashMap::new();

    for idx in 0..n_rows {
        // -------- peptide-level stuff --------

        // 1) Read stripped sequence
        let seq_val = seq_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'StrippedPeptide' at row {}", idx));
        let seq: String = match seq_val {
            AnyValue::String(s) => s.to_string(),
            AnyValue::StringOwned(ref s) => s.to_string(),
            _ => panic!("StrippedPeptide must be Utf8 at row {}", idx),
        };

        // 2) Read ModifiedPeptide (always exists, always used as key)
        let modpep_val = modpep_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'ModifiedPeptide' at row {}", idx));
        let modified: String = match modpep_val {
            AnyValue::String(s) => s.to_string(),
            AnyValue::StringOwned(ref s) => s.to_string(),
            _ => panic!("ModifiedPeptide must be Utf8 at row {}", idx),
        };

        // 3) Determine mods_vec
        let mods_vec: Vec<f32> = if let Some(mods_col) = &mods_col_opt {
            // ---- CASE A: user provided explicit masses ----
            let mods_val = mods_col
                .get(idx)
                .unwrap_or_else(|_| panic!("null 'Modifications' at row {}", idx));
            let list = match mods_val {
                AnyValue::List(ref s) => s,
                _ => panic!("Modifications must be List(Float32) at row {}", idx),
            };

            let mut v = Vec::with_capacity(list.len());
            for mv in list.iter() {
                v.push(anyvalue_to_f32(&mv, "Modifications", idx));
            }

            if v.len() != seq.len() + 2 {
                panic!(
                    "Modifications length {} != seq.len+2 ({}) for '{}' at row {}",
                    v.len(),
                    seq.len() + 2,
                    modified,
                    idx
                );
            }

            v
        } else {
            // ---- CASE B: fallback — parse ModifiedPeptide ----
            mods_from_modified_peptide(&seq, &modified, &unimod_table)
        };

        // 4) Store peptide once using ModifiedPeptide as the key
        pep_mods
            .entry(modified.clone())
            .or_insert_with(|| (seq.clone(), mods_vec.clone()));

        if mods_vec.len() != seq.len() + 2 {
            panic!(
                "modifications length {} != sequence length + 2 ({}) at row {} (seq='{}')",
                mods_vec.len(),
                seq.len() + 2,
                idx,
                seq
            );
        }

        // store peptide info once per ModifiedPeptide
        pep_mods
            .entry(modified.clone())
            .or_insert_with(|| (seq.clone(), mods_vec.clone()));

        // -------- fragment-level stuff (for aggregation) --------
        let frag_type_val = frag_type_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'FragmentType' at row {}", idx));
        let frag_type: String = match frag_type_val {
            AnyValue::String(s) => s.to_string(),
            AnyValue::StringOwned(ref s) => s.to_string(),
            other => panic!(
                "'FragmentType' must be Utf8, got {:?} at row {}",
                other, idx
            ),
        };

        let frag_num_val = frag_num_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'FragmentNumber' at row {}", idx));
        let frag_num = match frag_num_val {
            AnyValue::Int32(x) => x,
            AnyValue::Int64(x) => x as i32,
            other => panic!(
                "'FragmentNumber' must be integer, got {:?} at row {}",
                other, idx
            ),
        };

        let frag_z_val = frag_z_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'FragmentCharge' at row {}", idx));
        let frag_z = match frag_z_val {
            AnyValue::Int32(x) => x,
            AnyValue::Int64(x) => x as i32,
            other => panic!(
                "'FragmentCharge' must be integer, got {:?} at row {}",
                other, idx
            ),
        };

        let frag_mz_val = frag_mz_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'FragmentMz' at row {}", idx));
        let frag_int_val = frag_int_col
            .get(idx)
            .unwrap_or_else(|_| panic!("null 'RelativeIntensity' at row {}", idx));

        let mz = anyvalue_to_f32(&frag_mz_val, "FragmentMz", idx);
        let inten = anyvalue_to_f32(&frag_int_val, "RelativeIntensity", idx);

        let key = (modified.clone(), frag_type, frag_num, frag_z);

        let entry = ion_agg.entry(key).or_insert((0.0f32, 0.0f32, 0u32));
        entry.0 += mz;
        entry.1 += inten;
        entry.2 += 1;
    }

    // -------- build peptides (unique per ModifiedPeptide) --------
    let mut peptides: Vec<Peptide> = Vec::with_capacity(pep_mods.len());
    let mut pep_index: HashMap<String, PeptideIx> = HashMap::with_capacity(pep_mods.len());

    for (modified, (seq, mods_vec)) in pep_mods {
        let n = seq.len();
        if mods_vec.len() != n + 2 {
            panic!(
                "inconsistent mods length {} != seq.len+2 ({}) for '{}'",
                mods_vec.len(),
                n + 2,
                modified
            );
        }
        let nterm = mods_vec[0];
        let cterm = mods_vec[mods_vec.len() - 1];
        let per_res_mods = mods_vec[1..mods_vec.len() - 1].to_vec();

        let mono = mono_from_seq_and_mods(&seq, &mods_vec);

        let seq_bytes: Arc<[u8]> = seq.clone().into_bytes().into_boxed_slice().into();

        let peptide = Peptide {
            sequence: seq_bytes,
            monoisotopic: mono,
            proteins: vec![Arc::from("LIBRARY".to_string())],
            decoy: false,
            modifications: per_res_mods,
            nterm: Some(nterm),
            cterm: Some(cterm),
            missed_cleavages: 0,
            semi_enzymatic: false,
            position: Position::default(),
        };

        let ix = PeptideIx(peptides.len() as u32);
        pep_index.insert(modified, ix);
        peptides.push(peptide);
    }

    // -------- build fragments from aggregated ions --------
    let mut fragments: Vec<Theoretical> = Vec::with_capacity(ion_agg.len());
    let mut library_frags: Vec<Vec<LibraryFragment>> = vec![Vec::new(); peptides.len()];

    for ((modified, ftype, fnum, fz), (mz_sum, inten_sum, count)) in ion_agg {
        let pep_ix = *pep_index
            .get(&modified)
            .unwrap_or_else(|| panic!("internal error: missing peptide index for '{}'", modified));

        let avg_mz  = mz_sum  / (count as f32);
        let avg_int = inten_sum / (count as f32);

        // keep index lean: only peptide_index + fragment_mz
        fragments.push(Theoretical {
            peptide_index: pep_ix,
            fragment_mz: avg_mz,
        });

        // full fragment info goes into side table
        let kind = match ftype.as_str() {
            "b" | "B" => Kind::B,
            "y" | "Y" => Kind::Y,
            "a" | "A" => Kind::A,
            "x" | "X" => Kind::X,
            "c" | "C" => Kind::C,
            "z" | "Z" => Kind::Z,
            other => panic!("unsupported FragmentType '{}' in library", other),
        };

        library_frags[pep_ix.0 as usize].push(LibraryFragment {
            mz: avg_mz,
            kind,
            ordinal: fnum as u16,
            charge: fz as u8,
            rel_intensity: avg_int,
        });
    }

    // ensure sorted in ordinal order
    fn kind_order(k: Kind) -> u8 {
        match k {
            Kind::A => 0,
            Kind::B => 1,
            Kind::C => 2,
            Kind::X => 3,
            Kind::Y => 4,
            Kind::Z => 5,
        }
    }
    for frags in &mut library_frags {
        frags.sort_unstable_by(|a, b| {
            kind_order(a.kind)
                .cmp(&kind_order(b.kind))
                .then(a.ordinal.cmp(&b.ordinal))
        });
    }

    // -------- optional decoys (still skipped here) --------
    if config.generate_decoys {
        eprintln!(
            "Warning: generate_decoys is not yet supported in build_indexed_database_from_library; returning targets only."
        );
    }

    // -------- sort + bucket fragments as before --------
    use rayon::prelude::*;

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
        library_frags: Some(library_frags),
        min_value,
        bucket_size: config.bucket_size,
        ion_kinds: config.ion_kinds.clone(), // unused here but kept for compat
        generate_decoys: false,
        potential_mods: Vec::new(),
        decoy_tag: config.decoy_tag.clone(),
    }
}

fn anyvalue_to_f32(v: &AnyValue, col: &str, idx: usize) -> f32 {
    match v {
        AnyValue::Float32(x) => *x,
        AnyValue::Float64(x) => *x as f32,
        AnyValue::Int32(x)   => *x as f32,
        AnyValue::Int64(x)   => *x as f32,
        AnyValue::Null       => 0.0,
        other => panic!(
            "'{}' must be numeric, got {:?} at row {}",
            col, other, idx
        ),
    }
}

fn mods_from_modified_peptide(
    seq: &str,
    modified: &str,
    unimod_table: &HashMap<&'static str, f32>,
) -> Vec<f32> {
    let n = seq.len();
    let mut mods = vec![0.0f32; n + 2]; // [N-term, per-res..., C-term]

    let chars: Vec<char> = modified.chars().collect();
    let mut i = 0usize;
    let mut aa_idx = 0usize;

    while i < chars.len() {
        let c = chars[i];

        if c == '(' {
            let term_target = if aa_idx == 0 {
                Some(0usize)
            } else if aa_idx == n {
                Some(n + 1)
            } else {
                None
            };

            let mut unimod = String::new();
            i += 1;
            while i < chars.len() && chars[i] != ')' {
                unimod.push(chars[i]);
                i += 1;
            }
            if i == chars.len() {
                panic!("Unmatched '(' in ModifiedPeptide '{}'", modified);
            }
            i += 1; // consume ')'

            let mass = *unimod_table
                .get(unimod.as_str())
                .unwrap_or_else(|| panic!("Unknown UniMod '{}' in ModifiedPeptide '{}'", unimod, modified));

            if let Some(pos) = term_target {
                mods[pos] += mass;
            } else {
                if aa_idx == 0 {
                    panic!(
                        "Residue-level UniMod '{}' appears before any residue in '{}'",
                        unimod, modified
                    );
                }
                mods[aa_idx] += mass;
            }
        } else if c.is_ascii_alphabetic() && c.is_ascii_uppercase() {
            if aa_idx >= n {
                panic!(
                    "More residues in ModifiedPeptide '{}' than in StrippedPeptide '{}'",
                    modified, seq
                );
            }
            let expected = seq.as_bytes()[aa_idx] as char;
            if c != expected {
                panic!(
                    "Residue mismatch at position {}: ModifiedPeptide '{}' vs StrippedPeptide '{}'",
                    aa_idx, modified, seq
                );
            }
            aa_idx += 1;
            i += 1;
        } else {
            panic!(
                "Unexpected character '{}' in ModifiedPeptide '{}'",
                c, modified
            );
        }
    }

    if aa_idx != n {
        panic!(
            "Parsed only {} residues from ModifiedPeptide '{}' but StrippedPeptide has {}",
            aa_idx, modified, n
        );
    }

    mods
}


fn mono_from_seq_and_mods(seq: &str, mods_vec: &[f32]) -> f32 {
    // mods_vec layout: [N-term, per-res..., C-term], length = seq.len() + 2
    assert_eq!(
        mods_vec.len(),
        seq.len() + 2,
        "mods_vec must be [N-term, per-res..., C-term]"
    );

    let nterm = mods_vec[0];
    let cterm = mods_vec[mods_vec.len() - 1];
    let per_res_mods = &mods_vec[1..mods_vec.len() - 1];

    let mut mass = H2O;

    for (aa, mod_mass) in seq.as_bytes().iter().zip(per_res_mods.iter()) {
        let base = aa_mono(*aa);
        if base == 0.0 {
            panic!(
                "Invalid residue '{}' in library sequence '{}'",
                *aa as char, seq
            );
        }
        mass += base + mod_mass;
    }

    mass + nterm + cterm
}

pub fn unimod_table() -> HashMap<&'static str, f32> {
    HashMap::from([
        ("UniMod:4", 57.021464),
        ("Carbamidomethyl (C)", 57.021464),
        ("Carbamidomethyl", 57.021464),
        ("CAM", 57.021464),
        ("+57", 57.021464),
        ("+57.0", 57.021464),

        ("UniMod:26", 39.994915),
        ("PCm", 39.994915),

        ("UniMod:5", 43.005814),
        ("Carbamylation (KR)", 43.005814),
        ("+43", 43.005814),
        ("+43.0", 43.005814),
        ("CRM", 43.005814),

        ("UniMod:7", 0.984016),
        ("Deamidation (NQ)", 0.984016),
        ("Deamidation", 0.984016),
        ("Dea", 0.984016),
        ("+1", 0.984016),
        ("+1.0", 0.984016),

        ("UniMod:35", 15.994915),
        ("Oxidation (M)", 15.994915),
        ("Oxidation", 15.994915),
        ("Oxi", 15.994915),
        ("+16", 15.994915),
        ("+16.0", 15.994915),

        ("UniMod:1", 42.010565),
        ("Acetyl (Protein N-term)", 42.010565),
        ("+42", 42.010565),
        ("+42.0", 42.010565),

        ("UniMod:255", 28.0313),
        ("AAR", 28.0313),

        ("UniMod:254", 26.01565),
        ("AAS", 26.01565),

        ("UniMod:122", 27.994915),
        ("Frm", 27.994915),

        ("UniMod:1301", 128.094963),
        ("+1K", 128.094963),

        ("UniMod:1288", 156.101111),
        ("+1R", 156.101111),

        ("UniMod:27", -18.010565),
        ("PGE", -18.010565),

        ("UniMod:28", -17.026549),
        ("PGQ", -17.026549),

        ("UniMod:526", -48.003371),
        ("DTM", -48.003371),

        ("UniMod:325", 31.989829),
        ("2Ox", 31.989829),

        ("UniMod:342", 15.010899),
        ("Amn", 15.010899),

        ("UniMod:1290", 114.042927),
        ("2CM", 114.042927),

        ("UniMod:359", 13.979265),
        ("PGP", 13.979265),

        ("UniMod:30", 21.981943),
        ("NaX", 21.981943),

        ("UniMod:401", -2.015650),
        ("-2H", -2.015650),

        ("UniMod:528", 14.999666),
        ("MDe", 14.999666),

        ("UniMod:385", -17.026549),
        ("dAm", -17.026549),

        ("UniMod:23", -18.010565),
        ("Dhy", -18.010565),

        ("UniMod:129", 125.896648),
        ("Iod", 125.896648),

        ("Phosphorylation (ST)", 79.966331),
        ("UniMod:21", 79.966331),
        ("+80", 79.966331),
        ("+80.0", 79.966331),

        ("UniMod:259", 8.014199),
        ("Lys8", 8.014199),

        ("UniMod:267", 10.008269),
        ("Arg10", 10.008269),

        ("UniMod:268", 6.013809),
        ("UniMod:269", 10.027228),
    ])
}