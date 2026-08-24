"""Test reading the combined_JMod_lib.tsv spectral library.

This library uses the JMod/DIA-NN column naming convention which differs from
the schema `from_library` expects:

    JMod column         -> expected column
    -----------------------------------------
    PeptideSequence     -> StrippedPeptide
    ProductMz           -> FragmentMz
    LibraryIntensity    -> RelativeIntensity

We rename the columns, build a `Modifications` list[float32] column from
`ModifiedPeptide` (UniMod tags), and then attempt to build an IndexedDatabase
so that any crash inside `from_library` is reproduced under pytest.
"""

from pathlib import Path

import polars as pl
import pytest

import peppy_sage as ps


LIB_PATH = Path(__file__).parent / "combined_JMod_lib.tsv"

# Minimal UniMod -> monoisotopic mass shift table. Extend as needed when new
# tags appear in the library. Anything not listed here will raise so we notice.
UNIMOD_MASSES = {
    "UniMod:1": 42.010565,    # Acetyl
    "UniMod:4": 57.021464,    # Carbamidomethyl
    "UniMod:5": 43.005814,    # Carbamyl
    "UniMod:7": 0.984016,     # Deamidated
    "UniMod:21": 79.966331,   # Phospho
    "UniMod:23": -18.010565,  # Dehydrated
    "UniMod:26": 39.994915,   # Pyro-carbamidomethyl
    "UniMod:27": -18.010565,  # Glu->pyro-Glu
    "UniMod:28": -17.026549,  # Gln->pyro-Glu
    "UniMod:35": 15.994915,   # Oxidation
    "UniMod:121": 114.042927, # GG (ubiquitin remnant)
    "UniMod:259": 8.014199,   # Label:13C(6)15N(2)
    "UniMod:267": 10.008269,  # Label:13C(6)15N(4)
    "UniMod:737": 229.162932, # TMT6plex
}


def mods_for_row(modified_peptide: str, stripped: str) -> list[float]:
    """Parse `ModifiedPeptide` into a [N-term, res1..resN, C-term] mass array."""
    n = len(stripped)
    mods = [0.0] * (n + 2)

    aa_idx = 0
    i = 0
    chars = modified_peptide
    while i < len(chars):
        c = chars[i]
        if c == "(":
            j = chars.find(")", i + 1)
            if j == -1:
                raise ValueError(
                    f"unmatched '(' in ModifiedPeptide {modified_peptide!r}"
                )
            tag = chars[i + 1:j]
            mass = UNIMOD_MASSES.get(tag)
            if mass is None:
                raise ValueError(
                    f"unknown UniMod tag {tag!r} in {modified_peptide!r}"
                )
            if aa_idx == 0:
                mods[0] += mass
            elif aa_idx >= n:
                mods[n + 1] += mass
            else:
                mods[aa_idx] += mass
            i = j + 1
        elif c.isalpha() and c.isupper():
            aa_idx += 1
            i += 1
        else:
            i += 1
    return mods


def _load_library() -> pl.DataFrame:
    df = pl.read_csv(str(LIB_PATH), separator="\t")

    df = df.rename({
        "PeptideSequence": "StrippedPeptide",
        "ProductMz": "FragmentMz",
        "LibraryIntensity": "RelativeIntensity",
    })

    # The JMod TSV stores PrecursorCharge etc. as floats (e.g. "3.0"), but
    # `from_library` panics if they aren't integers. Cast explicitly.
    df = df.with_columns([
        pl.col("PrecursorCharge").cast(pl.Int32),
        pl.col("FragmentCharge").cast(pl.Int32),
        pl.col("FragmentSeriesNumber").cast(pl.Int32),
    ])

    mod_lists = [
        mods_for_row(mp, sp)
        for mp, sp in zip(
            df["ModifiedPeptide"].to_list(),
            df["StrippedPeptide"].to_list(),
        )
    ]
    df = df.with_columns(
        pl.Series("Modifications", mod_lists, dtype=pl.List(pl.Float32))
    )
    return df


@pytest.mark.skipif(not LIB_PATH.exists(), reason="combined_JMod_lib.tsv missing")
def test_combined_jmod_lib_loads():
    """Smoke test: rename + load the JMod library without crashing."""
    df = _load_library()

    assert df.height > 0
    assert "StrippedPeptide" in df.columns
    assert "FragmentMz" in df.columns
    assert "RelativeIntensity" in df.columns
    assert "Modifications" in df.columns

    print(f"rows: {df.height}, unique peptides: {df['ModifiedPeptide'].n_unique()}")


@pytest.mark.skipif(not LIB_PATH.exists(), reason="combined_JMod_lib.tsv missing")
def test_combined_jmod_lib_build_database():
    """Try to build IndexedDatabase from the JMod library — reproduces the crash."""
    df = _load_library()

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

    print(f"peptide count: {len(db.peptides)}")
    print(f"fragment summary: {db.debug_fragment_summary()}")
    assert len(db.peptides) > 0


if __name__ == "__main__":
    test_combined_jmod_lib_loads()
    test_combined_jmod_lib_build_database()
