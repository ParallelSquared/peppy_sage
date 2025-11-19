import numpy as np
from pyteomics import mass

import peppy_sage as ps
#from peppy_sage import Precursor, Spectrum, Scorer, IndexedDatabase



def test_database_build():
    print("\n--- Database Build ---")

    # Create peptide(s)
    seq = "PEPTIDEK"
    #peptide = ps.Peptide(seq, mods=[5.0] + ([0.0] * (len(seq)) + [5.0]))  # calculates mass and creates PyPeptide internally
    peptide = ps.Peptide(seq, mods=[0.0, 5.0] + ([0.0] * (len(seq)-2) + [0.0, 5.0]))

    # Create indexed database
    db = ps.IndexedDatabase.from_peptides(
        peptides=[peptide],
        bucket_size=128,
        ion_kinds=["b", "y"],
        min_ion_index=0,
        generate_decoys=True,
        decoy_tag="rev_",
        peptide_min_mass=0.0,
        peptide_max_mass=5000.0,
    )


    print("Peptide count:", len(db.peptides))
    for p in db.peptides:
        print("  →", p.sequence)

    print("Fragment summary:", db.debug_fragment_summary())

    expected_mass = mass.calculate_mass(sequence="PEPTIDEK")
    print("Expected neutral mass (pyteomics):", expected_mass)
    print("Precursor neutral mass (from m/z):", (464.73478 - 1.007276466812) * 2)

    return db


def test_spectrum_build():
    print("\n--- Spectrum Build ---")

    # 1. Precursor setup
    precursor = ps.Precursor(mz=310.15896, charge=0, isolation_window=(-1000, 1000))

    # 2. Peak data
    proton_mass = 1.0072764
    mz_arr = sorted([
        98.06009, 227.10268, 324.15544, 425.20312, 538.28718, 653.31413,
        782.35672, 910.45168, 928.46225, 831.40948, 702.36689, 605.31413,
        504.26645, 391.18238, 276.15544, 147.11285
    ])
    mz_arr = [mz - proton_mass for mz in mz_arr]  # Convert to neutral masses

    # Testing mz_arr with dummy mod of 1000
    mz_arr = [mz + 5.0 for mz in mz_arr]

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


def test_scoring(db, spectrum):
    print("\n--- Scoring Test ---")

    scorer = ps.Scorer(
        precursor_tol_da=(-1, 1),
        fragment_tol_ppm=(-5, 5),
        annotate_matches=True,
        report_psms=10,
        chimera=False
    )

    scorer = ps.Scorer(
        precursor_tol_da=(-1,1), # TODO placeholder
        fragment_tol_ppm=(-10,10),
        min_isotope_err=0,
        max_isotope_err=0,
        wide_window=True,
        chimera=False,
        annotate_matches=True,
        report_psms=1
    )

    spectra = [spectrum]

    BATCH_SIZE = 10000
    feature_arrays = None # Start with None or PyFeatureArrays()

    for i in range(0, len(spectra), BATCH_SIZE):
        batch_spectra = spectra[i:i + BATCH_SIZE]

        # 1. Get the columnar results for the current batch
        batch_arrays = scorer.score_many(db, batch_spectra)

        if feature_arrays is None:
            # Initialize the total container with the first batch
            feature_arrays = batch_arrays
        else:
            # 2. Append the new batch's data using the highly efficient Rust-side extend
            feature_arrays.extend(batch_arrays)

    import pandas as pd

    col_names = feature_arrays.get_column_names()
    print(feature_arrays.frag_intensities)
    df = pd.DataFrame({name: getattr(feature_arrays, name) for name in col_names})

    print(df.loc[0, 'frag_intensities'])

if __name__ == "__main__":
    db = test_database_build()
    spectrum = test_spectrum_build()
    test_scoring(db, spectrum)