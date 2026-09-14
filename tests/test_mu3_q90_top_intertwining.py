from order15_mu3.scripts.z3_q90_top_commuting import intertwining_certificate


def test_top_q90_orbit_pairs_have_no_mu3_intertwiner() -> None:
    result = intertwining_certificate(timeout_seconds=60)
    assert result["monomial_orbits"] == 2
    assert all(row["result"] == "unsat" for row in result["pair_results"])
