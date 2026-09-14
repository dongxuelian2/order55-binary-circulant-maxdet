from maxdet.mu3 import (
    Eisenstein,
    OMEGA,
    ONE,
    ZERO,
    dephase,
    determinant,
    exponent_matrix,
    gram_from_exponents,
    inner_product_values,
)


def test_eisenstein_ring_identities_and_exact_division() -> None:
    assert OMEGA * OMEGA + OMEGA + ONE == ZERO
    z = Eisenstein(17, -9)
    assert z * z.conjugate() == Eisenstein(z.norm())
    assert (z * Eisenstein(4, 7)).divexact(Eisenstein(4, 7)) == z


def test_fraction_free_determinant() -> None:
    matrix = exponent_matrix([[0, 0, 0], [0, 1, 2], [0, 2, 1]])
    value = determinant(matrix)
    assert value.norm() == 27
    assert determinant(gram_from_exponents([[0, 0, 0], [0, 1, 2], [0, 2, 1]])) == Eisenstein(27)


def test_dephase_preserves_gram_determinant_norm() -> None:
    exponents = [[1, 2, 0], [2, 2, 1], [0, 1, 1]]
    normalized = dephase(exponents)
    assert normalized[0] == [0, 0, 0]
    assert [row[0] for row in normalized] == [0, 0, 0]
    assert determinant(exponent_matrix(exponents)).norm() == determinant(exponent_matrix(normalized)).norm()


def test_order_15_inner_product_catalogue() -> None:
    values = inner_product_values(15)
    assert len(values) == 136
    assert values[0]["value"] == ZERO
    assert values[0]["counts"] == ((5, 5, 5),)
    positive_norms = sorted({entry["norm"] for entry in values if entry["norm"]})
    assert positive_norms[:4] == [3, 9, 12, 21]
