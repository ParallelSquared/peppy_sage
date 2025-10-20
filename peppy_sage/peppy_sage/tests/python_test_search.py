import numpy as np
from pyteomics import mass

import peppy_sage as ps


def test_database_build():
    print("\n--- Database Build ---")

    # Create peptide(s)
    seq = "PEPTIDEK"
    peptide = ps.core.Peptide(seq)  # calculates mass and creates PyPeptide internally

    # Create indexed database
    db = Database.from_peptides(
        peptides=[peptide],
        bucket_size=128,
        ion_kinds=["b", "y"],
        min_ion_index=2,
        generate_decoys=False,
        decoy_tag="rev_",
        peptide_min_mass=0.0,
        peptide_max_mass=1500.0,
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
    precursor = Precursor(mz=310.15896, charge=3, isolation_window=(1.5, 1.5))

    # 2. Peak data
    proton_mass = 1.0072764
    mz_arr = sorted([
        98.06009, 227.10268, 324.15544, 425.20312, 538.28718, 653.31413,
        782.35672, 910.45168, 928.46225, 831.40948, 702.36689, 605.31413,
        504.26645, 391.18238, 276.15544, 147.11285
    ])
    mz_arr = [mz - proton_mass for mz in mz_arr]  # Convert to neutral masses
    int_arr = mz_arr  # dummy intensities for test

    # 3. Build processed spectrum
    spectrum = Spectrum(
        id="Scan_100",
        file_id=1,
        scan_start_time=10.5,
        mz_array=mz_arr,
        intensity_array=int_arr,
        precursors=[precursor],
    )

    print(f"Spectrum ID: {spectrum.id}")
    print(f"Peak count: {len(spectrum.peaks)}")
    print(f"Precursor MZ: {spectrum.precursors[0].mz:.5f}")

    return spectrum


def test_scoring(db, spectrum):
    print("\n--- Scoring Test ---")

    scorer = Scorer(
        precursor_tol_da=(-1, 1),
        fragment_tol_ppm=(-5, 5),
        annotate_matches=True,
        report_psms=10,
    )

    features = scorer.score(db, spectrum)
    print(f"Total features returned: {len(features)}")


if __name__ == "__main__":
    db = test_database_build()
    spectrum = test_spectrum_build()
    test_scoring(db, spectrum)