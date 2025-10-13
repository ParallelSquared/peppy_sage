import peppy_sage
from pyteomics import mass

def test_database_build():
    seq1 = "ACDEFGK"
    seq2 = "PEPTIDEK"
    mass1 = mass.calculate_mass(sequence = seq1)
    mass2 = mass.calculate_mass(sequence = seq2)
    pep1 = peppy_sage.PyPeptide(seq1, mass1, [0.0]*7)
    pep2 = peppy_sage.PyPeptide(seq2, mass2, [0.0]*7)
    b = peppy_sage.PyKind.B()
    y = peppy_sage.PyKind.Y()

    db = peppy_sage.PyIndexedDatabase.from_peptides_and_config(
        [pep1, pep2],
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

if __name__ == "__main__":
    test_database_build()