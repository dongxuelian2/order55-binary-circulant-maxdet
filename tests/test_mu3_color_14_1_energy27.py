from order15_mu3.scripts.verify_color_14_1_energy27 import orbit_representatives
from maxdet.mu3 import Eisenstein, inner_product_values


def test_norm9_has_two_cube_root_orbits() -> None:
    values = tuple(entry["value"] for entry in inner_product_values(15) if entry["norm"] == 9)
    assert len(orbit_representatives(values)) == 2
    assert Eisenstein(-3) in values and Eisenstein(3) in values
