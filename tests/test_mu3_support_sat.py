from order15_mu3.scripts.sat_support_graph import SHAPES, build


def test_q45_support_models_build_with_exact_order() -> None:
    for shape, (edges, isolates) in SHAPES.items():
        support_vertices = max(max(edge) for edge in edges) + 1
        assert isolates + support_vertices == 15, shape
        cnf, pool, _x, shifted_edges, representatives, selectors = build(shape)
        assert cnf.clauses
        assert pool.top > 0
        assert len(shifted_edges) == 5
        assert len(representatives) == len(selectors) == 5
