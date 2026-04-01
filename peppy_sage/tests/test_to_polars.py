import numpy as np
from pyteomics import mass

import peppy_sage as ps
import polars as pl


def test_database_build():
    print("\n--- Database Build ---")

    df = pl.read_csv('trimmed_lib.tsv', separator='\t')
    print(df)
    # Create indexed database
    db = ps.IndexedDatabase.from_library(
        library=df,
        bucket_size=128,
        ion_kinds=["b", "y"],
        min_ion_index=0,
        generate_decoys=False,
        decoy_tag="rev_",
        peptide_min_mass=0.0,
        peptide_max_mass=5000.0,
    )

    print("Peptide count:", len(db.peptides))
    for p in db.peptides:
        print("  →", p.sequence)

    print("Fragment summary:", db.debug_fragment_summary())

    return db


def test_spectrum_build():
    print("\n--- Spectrum Build ---")

    proton_mass = 1.0072764
    # 1. Precursor setup
    #precursor = ps.Precursor(mz=(956.55240+proton_mass)/2, charge=0, isolation_window=(-2.4, 2.4))
    precursor = ps.Precursor(mz=(956.55240+proton_mass)/2, charge=0, isolation_window=(-100002.4, 100002.4))

    # 2. Peak data
    mz_arr = [72.04444,
     143.08155,
     214.11866,
     285.15578,
     356.19289,
     427.23000,
     498.26712,
     611.35118,
     739.40976,
     810.44687,
     938.54183,
     956.55240,
     885.51528,
     814.47817,
     743.44106,
     672.40394,
     601.36683,
     530.32972,
     459.29260,
     346.20854,
     218.14996,
     147.11285]

    mz_arr = sorted(mz_arr)
    int_arr = mz_arr  # dummy intensities for test

    # 3. Build processed spectrum
    spectrum = ps.Spectrum(
        id="Scan_100",
        file_id=1,
        scan_start_time=10.5,
        mz_array=mz_arr,
        intensity_array=int_arr,
        precursors=[precursor],
    )

    print(f"Spectrum ID: {spectrum.id}")
    print(f"Peak count: {len(spectrum.peaks)}")
    print(f"Peaks: {spectrum.peaks}")
    print(f"Precursor MZ: {spectrum.precursors[0].mz:.5f}")

    return spectrum


def test_scoring_to_polars(db, spectrum):
    """Test scoring with to_polars() method for fast DataFrame conversion."""
    print("\n--- Scoring Test (to_polars) ---")

    scorer = ps.Scorer(
        precursor_tol_da=(-1, 1),
        fragment_tol_ppm=(-10, 10),
        min_isotope_err=0,
        max_isotope_err=0,
        wide_window=True,
        chimera=False,
        annotate_matches=True,
        report_psms=100
    )

    spectra = [spectrum]

    BATCH_SIZE = 10000
    feature_arrays = None

    for i in range(0, len(spectra), BATCH_SIZE):
        batch_spectra = spectra[i:i + BATCH_SIZE]
        batch_arrays = scorer.score_many(db, batch_spectra)

        if feature_arrays is None:
            feature_arrays = batch_arrays
        else:
            feature_arrays.extend(batch_arrays)

    # Use the new to_polars() method - DataFrame built entirely in Rust
    print("\nConverting to Polars DataFrame using to_polars()...")
    df = feature_arrays.to_polars()

    print(f"\nDataFrame shape: {df.shape}")
    print(f"Columns: {df.columns}")
    print(f"\nFirst few rows:")
    print(df.head())

    # Verify key columns exist and have correct types
    print(f"\nColumn dtypes:")
    print(df.schema)

    # Check specific columns
    print(f"\nexpmass: {df['expmass']}")
    print(f"calcmass: {df['calcmass']}")
    print(f"modified_peptide: {df['modified_peptide']}")

    # Check list columns (fragments)
    print(f"\nfrag_charges (list column): {df['frag_charges']}")
    print(f"frag_kinds (list column): {df['frag_kinds']}")

    return df


def test_compare_methods(db, spectrum):
    """Compare old method vs to_polars() for correctness."""
    print("\n--- Comparing DataFrame creation methods ---")

    import time

    scorer = ps.Scorer(
        precursor_tol_da=(-1, 1),
        fragment_tol_ppm=(-10, 10),
        min_isotope_err=0,
        max_isotope_err=0,
        wide_window=True,
        chimera=False,
        annotate_matches=True,
        report_psms=100
    )

    spectra = [spectrum]
    feature_arrays = scorer.score_many(db, spectra)

    # Method 1: Old way (getattr loop)
    start = time.perf_counter()
    col_names = feature_arrays.get_column_names()
    df_old = pl.DataFrame({name: getattr(feature_arrays, name) for name in col_names})
    time_old = time.perf_counter() - start
    print(f"Old method time: {time_old*1000:.2f} ms")

    # Method 2: New to_polars()
    start = time.perf_counter()
    df_new = feature_arrays.to_polars()
    time_new = time.perf_counter() - start
    print(f"to_polars() time: {time_new*1000:.2f} ms")

    # Compare results
    print(f"\nOld DataFrame shape: {df_old.shape}")
    print(f"New DataFrame shape: {df_new.shape}")

    # Check that scalar columns match
    scalar_cols = ['file_id', 'spec_id', 'hyperscore', 'charge', 'sequence']
    for col in scalar_cols:
        if col in df_old.columns and col in df_new.columns:
            old_vals = df_old[col].to_list()
            new_vals = df_new[col].to_list()
            match = old_vals == new_vals
            print(f"Column '{col}' matches: {match}")

    print(f"\nSpeedup: {time_old/time_new:.2f}x")


if __name__ == "__main__":
    db = test_database_build()
    spectrum = test_spectrum_build()

    # Test the new to_polars method
    df = test_scoring_to_polars(db, spectrum)

    # Compare old vs new method
    test_compare_methods(db, spectrum)
