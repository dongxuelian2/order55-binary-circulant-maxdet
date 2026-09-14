from maxdet.mu3 import inner_product_values
from order15_mu3.scripts.verify_color_15_low_energy import (
    connected_types,
    has_full_support_kernel,
    q54_k4_final_obstruction,
    q63_five_vertex_final_obstruction,
    q72_no_isolate_obstruction,
    q9_layer,
    weighted_layer,
)
from maxdet.mu3 import Eisenstein


def test_all_same_color_layers_27_and_36_have_no_norm_survivor() -> None:
    norm9 = tuple(entry["value"] for entry in inner_product_values(15) if entry["norm"] == 9)
    types = connected_types(4, norm9)
    assert q9_layer(types, 3)["arithmetic_survivors"] == []
    assert q9_layer(types, 4)["arithmetic_survivors"] == []


def test_full_support_kernel_exact_field_elimination() -> None:
    zero = Eisenstein()
    one = Eisenstein(1)
    edge = [[zero, one], [one, zero]]
    cycle4 = [
        [zero, one, zero, one],
        [one, zero, one, zero],
        [zero, one, zero, one],
        [one, zero, one, zero],
    ]
    assert not has_full_support_kernel(edge)
    assert has_full_support_kernel(cycle4)


def test_q45_survivors_fail_intertwining_kernel_condition() -> None:
    norm9 = tuple(entry["value"] for entry in inner_product_values(15) if entry["norm"] == 9)
    result = q9_layer(connected_types(5, norm9), 5)
    assert len(result["arithmetic_survivors"]) == 8
    assert result["intertwining_arithmetic_survivors"] == []


def test_q54_last_k4_profile_has_negative_residual_gram() -> None:
    norm9 = tuple(entry["value"] for entry in inner_product_values(15) if entry["norm"] == 9)
    result = q9_layer(connected_types(6, norm9), 6)
    assert len(result["arithmetic_survivors"]) == 18
    certificate = q54_k4_final_obstruction(result)
    assert certificate["all_ones_eigenvalue"] == -29


def test_q63_last_five_vertex_profile_has_negative_residual_gram() -> None:
    norm9 = tuple(entry["value"] for entry in inner_product_values(15) if entry["norm"] == 9)
    result = q9_layer(connected_types(7, norm9, kernel_only=True), 7)
    assert result["intertwining_arithmetic_survivors"] == [325793054443359375]
    certificate = q63_five_vertex_final_obstruction(result)
    assert certificate["all_ones_eigenvalue"] == -35


def test_q72_unique_no_isolate_support_fails_norm_test() -> None:
    certificate = q72_no_isolate_obstruction()
    assert certificate["eisenstein_norm"] is False


def test_q72_last_mixed_k4_profile_has_negative_residual_gram() -> None:
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    types = connected_types(6, norm9, kernel_only=True)
    result = weighted_layer(types, (3, 1, 1, 1, 1, 1), catalogue, kernel_only=True)
    assert result["intertwining_arithmetic_survivors"] == [325793054443359375]
    assert q54_k4_final_obstruction(result)["all_ones_eigenvalue"] == -29


def test_all_same_color_energy54_mixed_partitions_are_enumerated() -> None:
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    types = connected_types(4, norm9)
    results = {weights: weighted_layer(types, weights, catalogue) for weights in ((3, 3), (4, 1, 1), (3, 1, 1, 1))}
    assert results[(3, 3)]["profiles"] == []
    assert all(result["arithmetic_survivors"] == [] for result in results.values())
