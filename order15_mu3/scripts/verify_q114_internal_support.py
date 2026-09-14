"""Exact 14-row internal determinant maxima used in the pre-Q126 audit."""

from __future__ import annotations

import json
from functools import lru_cache

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy90 import weight_partitions
    from order15_mu3.scripts.verify_color_15_low_energy import connected_types, weighted_layer
except ModuleNotFoundError:
    from verify_color_15_energy90 import weight_partitions
    from verify_color_15_low_energy import connected_types, weighted_layer


EXPECTED_MAXIMA = {
    45: 24288514453125000,
    54: 23914845000000000,
    63: 22958251200000000,
    72: 22121231625000000,
}


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    types = connected_types(8, norm9)
    rows = []
    for energy, expected in EXPECTED_MAXIMA.items():
        candidates = {}
        for weights in weight_partitions(energy // 9):
            layer = weighted_layer(types, weights, catalogue, order=14)
            candidates[str(weights)] = layer["maximum_determinant"]
        maximum = max(candidates.values())
        assert maximum == expected
        rows.append({
            "internal_energy": energy,
            "partition_maxima": candidates,
            "maximum_determinant": maximum,
        })
    return {
        "order": 14,
        "connected_component_types": len(types),
        "energy_rows": rows,
        "theorem": "The displayed values are the exact internal-block determinant maxima.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
