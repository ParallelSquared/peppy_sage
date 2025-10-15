import peppy_sage
import numpy as np
from pyteomics import mass

def test_database_build():
    seq1 = "ACDEFGK"
    seq2 = "PEPTIDEK"
    seq3 = "PEPTIQEK"
    seq4 = "PEPTQDEK"
    seq5 = "PEPQIDEK"
    seq6 = "PEQTIDEK"
    mass1 = mass.calculate_mass(sequence = seq1)
    mass2 = mass.calculate_mass(sequence = seq2)
    mass3 = mass.calculate_mass(sequence = seq3)
    mass4 = mass.calculate_mass(sequence = seq4)
    mass5 = mass.calculate_mass(sequence = seq5)
    mass6 = mass.calculate_mass(sequence = seq6)
    pep1 = peppy_sage.PyPeptide(seq1, mass1, [0.0]*7)
    pep2 = peppy_sage.PyPeptide(seq2, mass2, [0.0]*8)
    pep3 = peppy_sage.PyPeptide(seq3, mass3, [0.0]*8)
    pep4 = peppy_sage.PyPeptide(seq4, mass4, [0.0]*8)
    pep5 = peppy_sage.PyPeptide(seq5, mass5, [0.0]*8)
    pep6 = peppy_sage.PyPeptide(seq6, mass6, [0.0]*8)
    b = peppy_sage.PyKind.B()
    y = peppy_sage.PyKind.Y()

    db = peppy_sage.PyIndexedDatabase.from_peptides_and_config(
        [pep1, pep2, pep3, pep4, pep5, pep6],
        128,
        [b, y],
        2,
        False,
        "rev_",
        500.0,
        1500.0
    )

    print("Peptide count:", len(db.peptides))
    for p in db.peptides:
        print("  →", p.sequence)

    print(db.debug_fragment_summary())

    print("Expected neutral mass (from pyteomics):", mass.calculate_mass(sequence="PEPTIDEK"))
    print("Precursor neutral mass (from m/z):", (464.73478 - 1.007276466812) * 2)

    return db

def test_spectrum_build():
    # 1. Create a minimal Precursor object (Input for the spectrum)
    # This precursor simulates a parent ion with a mass of 1000.0 m/z and charge 2.
    precursor = peppy_sage.PyPrecursor(
        mz=310.15896,
        charge=3,
        isolation_window=peppy_sage.PyTolerance.Da(1.5, 1.5) # Example: 3.0 Da window
    )

    # 2. Define the fragment ions (peaks) in the MS2 scan
    # Note: These masses are usually proton-removed (monoisotopic mass) in sage_core Peak struct
    #mz_arr = np.array([98.06009, 227.10268, 324.15544, 425.20312, 538.28718, 653.31413, 782.35672, 910.45168])
    proton_mass = 1.007276466812

    mz_arr = [98.06009, 227.10268, 324.15544, 425.20312, 538.28718, 653.31413, 782.35672, 910.45168,
              928.46225, 831.40948, 702.36689, 605.31413, 504.26645, 391.18238, 276.15544, 147.11285]
    mz_arr = [mz - proton_mass for mz in mz_arr]
    print(mz_arr)
    #int_arr = [50.0, 100.0, 75.0, 200.0]
    int_arr = [mz for mz in mz_arr]


    # 3. Create the ProcessedSpectrum object
    # PROCESSED SPECTRUM
    spectrum = peppy_sage.PyProcessedSpectrum(
        "Scan_100",         # id
        1,                  # file_id
        10.5,               # scan_start_time
        mz_arr,             # mz_array
        int_arr,            # intensity_array
        [precursor],        # precursors
        sum(int_arr)        # total_ion_current
    )

    print("\n--- Spectrum Construction Test ---")
    print("Spectrum ID:", spectrum.id)
    print("Peak count:", len(spectrum.peaks)) # Assuming you implement a getter for peaks
    print("Precursor MZ:", spectrum.precursors[0].mz)

    return spectrum

def test_scoring(db, spectrum):
    print("\n--- Scoring Test ---")

    #db = test_database_build()
    #spectrum = test_spectrum_build()

    precursor_tol = peppy_sage.PyTolerance.Da(-100, 100.0)
    fragment_tol = peppy_sage.PyTolerance.Ppm(-100, 100)

    scorer = peppy_sage.PyScorer(
        precursor_tol,
        fragment_tol,
        wide_window=False,
        chimera=False,
        report_psms=100,
        min_isotope_err=-1,
        max_isotope_err=3,
        min_precursor_charge=1,
        max_precursor_charge=4,
        min_matched_peaks=0,
        annotate_matches=False
    )

    features = scorer.score_spectra(db, spectrum)

    print(f"Features returned: {len(features)}")

    for i, f in enumerate(features):
        print(db.peptides[f.peptide_idx].sequence)
        print(f.fragments.to_dict())
        print(f" Feature {i}: hyperscore={f.hyperscore:.2f}, delta_mass={f.delta_mass:.4f}, matched_peaks={f.matched_peaks}")
        print(f.__repr__())


if __name__ == "__main__":
    db = test_database_build()
    s = test_spectrum_build()
    test_scoring(db, s)