"""Test that library-mode charge state tracking works correctly.

When a peptide appears in the library at multiple charge states (e.g., z=2 and z=3),
each charge state should be indexed separately. A spectrum with a known precursor
charge should only match the library entry at that charge state.
"""

import peppy_sage as ps
import polars as pl
from pyteomics import mass as pyteomass


PROTON = 1.0072764
SEQUENCE = "PEPTIDEK"


def _neutral_mass():
    return pyteomass.calculate_mass(sequence=SEQUENCE)


def _precursor_mz(neutral_mass, charge):
    """Standard precursor m/z: (M + z*H) / z"""
    return (neutral_mass + charge * PROTON) / charge


def _b_ion_mzs(sequence):
    """Compute singly-charged b-ion m/z values for a sequence."""
    mzs = []
    cumulative = 0.0
    for i in range(1, len(sequence)):
        cumulative += pyteomass.calculate_mass(sequence=sequence[i - 1:i])
        # b-ion = cumulative residue mass + proton
        mzs.append(cumulative + PROTON)
    return mzs


def _y_ion_mzs(sequence):
    """Compute singly-charged y-ion m/z values for a sequence."""
    mzs = []
    for i in range(1, len(sequence)):
        suffix = sequence[len(sequence) - i:]
        mzs.append(pyteomass.calculate_mass(sequence=suffix) + PROTON)
    return mzs


def _build_library_df(charges):
    """Build a minimal library DataFrame with PEPTIDEK at the given charge states.

    Each entry gets synthetic b and y ions so the scorer has fragments to match.
    """
    n = len(SEQUENCE)
    mods = [0.0] * (n + 2)  # no modifications

    neutral = _neutral_mass()

    rows = []
    for z in charges:
        prec_mz = _precursor_mz(neutral, z)
        b_mzs = _b_ion_mzs(SEQUENCE)
        y_mzs = _y_ion_mzs(SEQUENCE)

        for i, mz in enumerate(b_mzs):
            rows.append({
                "StrippedPeptide": SEQUENCE,
                "ModifiedPeptide": SEQUENCE,
                "PrecursorCharge": z,
                "PrecursorMz": prec_mz,
                "FragmentType": "b",
                "FragmentSeriesNumber": i + 1,
                "FragmentCharge": 1,
                "FragmentMz": mz,
                "RelativeIntensity": 100.0,
            })

        for i, mz in enumerate(y_mzs):
            rows.append({
                "StrippedPeptide": SEQUENCE,
                "ModifiedPeptide": SEQUENCE,
                "PrecursorCharge": z,
                "PrecursorMz": prec_mz,
                "FragmentType": "y",
                "FragmentSeriesNumber": i + 1,
                "FragmentCharge": 1,
                "FragmentMz": mz,
                "RelativeIntensity": 100.0,
            })

    df = pl.DataFrame(rows)
    # Add explicit Modifications column
    df = df.with_columns(
        pl.Series("Modifications", [mods] * len(rows), dtype=pl.List(pl.Float32))
    )
    return df


def _build_db(charges):
    df = _build_library_df(charges)
    return ps.IndexedDatabase.from_library(
        library=df,
        bucket_size=128,
        ion_kinds=["b", "y"],
        min_ion_index=0,
        generate_decoys=False,
        decoy_tag="rev_",
        peptide_min_mass=0.0,
        peptide_max_mass=5000.0,
    )


def _build_spectrum(charge):
    """Build a synthetic spectrum for PEPTIDEK at the given charge state."""
    neutral = _neutral_mass()
    mz = _precursor_mz(neutral, charge)

    precursor = ps.Precursor(mz=mz, charge=charge, isolation_window=(-2.4, 2.4))

    # Combine b and y ions as the observed peaks
    peak_mzs = sorted(_b_ion_mzs(SEQUENCE) + _y_ion_mzs(SEQUENCE))
    peak_intensities = [100.0] * len(peak_mzs)

    return ps.Spectrum(
        id=f"Scan_z{charge}",
        file_id=1,
        scan_start_time=10.0,
        mz_array=peak_mzs,
        intensity_array=peak_intensities,
        precursors=[precursor],
    )


def _score(db, spectrum):
    scorer = ps.Scorer(
        precursor_tol_da=(-0.5, 0.5),
        fragment_tol_ppm=(-10, 10),
        min_isotope_err=0,
        max_isotope_err=0,
        min_precursor_charge=1,
        max_precursor_charge=5,
        annotate_matches=True,
        report_psms=10,
        chimera=False,
    )
    spectra = [spectrum]
    return scorer.score_many(db, spectra).to_polars()


def test_separate_charge_state_entries():
    """Same peptide at z=2 and z=3 should produce two separate index entries."""
    db = _build_db(charges=[2, 3])
    print(f"\n  Peptide count: {len(db.peptides)}")
    for i, p in enumerate(db.peptides):
        print(f"    Peptide {i}: seq={p.sequence}")
    assert len(db.peptides) == 2, f"Expected 2 peptide entries, got {len(db.peptides)}"


def test_z2_spectrum_matches_z2_entry_only():
    """A z=2 spectrum should only match the z=2 library entry, not z=3."""
    db = _build_db(charges=[2, 3])
    spectrum_z2 = _build_spectrum(charge=2)

    neutral = _neutral_mass()
    mz = _precursor_mz(neutral, 2)
    print(f"\n  Neutral mass: {neutral:.5f}")
    print(f"  Precursor m/z (z=2): {mz:.5f}")
    print(f"  DB peptides: {[p.sequence for p in db.peptides]}")

    df = _score(db, spectrum_z2)
    print(f"  Matches: {len(df)}")
    if len(df) > 0:
        print(f"  charges: {df['charge'].to_list()}")
        print(f"  calcmass: {df['calcmass'].to_list()}")
        print(f"  expmass: {df['expmass'].to_list()}")
        print(f"  hyperscore: {df['hyperscore'].to_list()}")
        print(f"  matched_peaks: {df['matched_peaks'].to_list()}")

    assert len(df) >= 1, "Expected at least one match for z=2 spectrum"
    charges = df["charge"].to_list()
    assert all(z == 2 for z in charges), (
        f"z=2 spectrum matched non-z=2 entries: charges={charges}"
    )


def test_z3_spectrum_matches_z3_entry_only():
    """A z=3 spectrum should only match the z=3 library entry, not z=2."""
    db = _build_db(charges=[2, 3])
    spectrum_z3 = _build_spectrum(charge=3)

    neutral = _neutral_mass()
    mz = _precursor_mz(neutral, 3)
    print(f"\n  Neutral mass: {neutral:.5f}")
    print(f"  Precursor m/z (z=3): {mz:.5f}")
    print(f"  DB peptides: {[p.sequence for p in db.peptides]}")

    df = _score(db, spectrum_z3)
    print(f"  Matches: {len(df)}")
    if len(df) > 0:
        print(f"  charges: {df['charge'].to_list()}")
        print(f"  calcmass: {df['calcmass'].to_list()}")
        print(f"  hyperscore: {df['hyperscore'].to_list()}")

    assert len(df) >= 1, "Expected at least one match for z=3 spectrum"
    charges = df["charge"].to_list()
    assert all(z == 3 for z in charges), (
        f"z=3 spectrum matched non-z=3 entries: charges={charges}"
    )


def test_calcmass_is_neutral_mass():
    """The reported calcmass should be the neutral mass, not precursor m/z."""
    db = _build_db(charges=[2])
    spectrum = _build_spectrum(charge=2)
    df = _score(db, spectrum)

    print(f"\n  Matches: {len(df)}")
    if len(df) > 0:
        print(f"  calcmass: {df['calcmass'].to_list()}")
        print(f"  expected neutral: {_neutral_mass():.5f}")

    assert len(df) >= 1, "Expected at least one match"
    expected_neutral = _neutral_mass()
    reported_calcmass = df["calcmass"][0]
    # Allow 0.01 Da tolerance for float precision
    assert abs(reported_calcmass - expected_neutral) < 0.01, (
        f"calcmass {reported_calcmass:.4f} doesn't match expected neutral mass "
        f"{expected_neutral:.4f}"
    )


if __name__ == "__main__":
    test_separate_charge_state_entries()
    print("PASS: test_separate_charge_state_entries")

    test_z2_spectrum_matches_z2_entry_only()
    print("PASS: test_z2_spectrum_matches_z2_entry_only")

    test_z3_spectrum_matches_z3_entry_only()
    print("PASS: test_z3_spectrum_matches_z3_entry_only")

    test_calcmass_is_neutral_mass()
    print("PASS: test_calcmass_is_neutral_mass")

    print("\nAll charge state tracking tests passed.")
