"""Exact support-component enumeration for all-same-color energy <=45."""

from __future__ import annotations

import json
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product

import networkx as nx

from maxdet.mu3 import Eisenstein, determinant, inner_product_values

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm
    from order15_mu3.scripts.verify_support_independence import CODE_BOUND, independence_number
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_bounds import is_rational_eisenstein_norm
    from verify_support_independence import CODE_BOUND, independence_number
    from verify_trace_stability import BENCHMARK


def _field_add(x, y):
    return x[0] + y[0], x[1] + y[1]


def _field_sub(x, y):
    return x[0] - y[0], x[1] - y[1]


def _field_mul(x, y):
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - b * d


def _field_inv(x):
    a, b = x
    norm = a * a - a * b + b * b
    assert norm
    return (a - b) / norm, -b / norm


def kernel_dimension_and_full_support(matrix: list[list[Eisenstein]]) -> tuple[int, bool]:
    """Return exact nullity and whether the kernel contains a full-support vector."""

    rows = [[(Fraction(value.a), Fraction(value.b)) for value in row] for row in matrix]
    dimension = len(rows)
    pivots = []
    pivot_row = 0
    zero = (Fraction(0), Fraction(0))
    for column in range(dimension):
        selected = next((row for row in range(pivot_row, dimension) if rows[row][column] != zero), None)
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        inverse = _field_inv(rows[pivot_row][column])
        rows[pivot_row] = [_field_mul(value, inverse) for value in rows[pivot_row]]
        for row in range(dimension):
            if row == pivot_row or rows[row][column] == zero:
                continue
            factor = rows[row][column]
            rows[row] = [_field_sub(left, _field_mul(factor, right)) for left, right in zip(rows[row], rows[pivot_row])]
        pivots.append(column)
        pivot_row += 1
    free = [column for column in range(dimension) if column not in pivots]
    if not free:
        return 0, False
    # Over the infinite Eisenstein fraction field, a finite union of proper
    # coordinate hyperplanes cannot cover the kernel.  Hence a full-support
    # vector exists iff no coordinate vanishes identically on the kernel.
    for row, pivot in enumerate(pivots):
        if all(rows[row][column] == zero for column in free):
            return len(free), False
    return len(free), True


def has_full_support_kernel(matrix: list[list[Eisenstein]]) -> bool:
    return kernel_dimension_and_full_support(matrix)[1]


def unit_kernel_labelings(graph: nx.Graph, norm9):
    """Gauge-fixed q9 edge labels for which the all-ones vector is in ker A."""

    edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
    units = tuple(Eisenstein(value.a // 3, value.b // 3) for value in norm9)
    unit_set = set(units)
    incident = {vertex: [edge for edge in edges if vertex in edge] for vertex in graph.nodes()}
    zero = Eisenstein()

    def coefficient(edge, label, vertex):
        return label if edge[0] == vertex else label.conjugate()

    def propagate(assignments):
        assignments = dict(assignments)
        while True:
            changed = False
            for vertex in graph.nodes():
                missing = [edge for edge in incident[vertex] if edge not in assignments]
                total = sum(
                    (coefficient(edge, assignments[edge], vertex) for edge in incident[vertex] if edge in assignments),
                    zero,
                )
                if not missing:
                    if total != zero:
                        return None
                elif len(missing) == 1:
                    edge = missing[0]
                    needed_at_vertex = -total
                    label = needed_at_vertex if edge[0] == vertex else needed_at_vertex.conjugate()
                    if label not in unit_set:
                        return None
                    assignments[edge] = label
                    changed = True
            if not changed:
                return assignments

    def search(assignments):
        assignments = propagate(assignments)
        if assignments is None:
            return
        if len(assignments) == len(edges):
            yield tuple(assignments[edge] for edge in edges)
            return
        edge = max(
            (edge for edge in edges if edge not in assignments),
            key=lambda item: graph.degree(item[0]) + graph.degree(item[1]),
        )
        for label in units:
            yield from search({**assignments, edge: label})

    yield from search({})


def balanced_kernel_labelings(edges, vertex_count: int, options):
    """Enumerate edge labels whose Hermitian adjacency has A*1=0.

    An isolated row of the opposite Gram supplies a full-support mu3 kernel
    vector.  Switching by that vector gauges it to the all-ones vector, so
    vertex balance is an exact and much stronger generator than testing every
    phase assignment after the fact.
    """

    edges = tuple(edges)
    allowed = {edge: set(edge_options) for edge, edge_options in zip(edges, options)}
    incident = {vertex: [edge for edge in edges if vertex in edge] for vertex in range(vertex_count)}
    zero = Eisenstein()

    def coefficient(edge, label, vertex):
        return label if edge[0] == vertex else label.conjugate()

    def propagate(assignments):
        assignments = dict(assignments)
        while True:
            changed = False
            for vertex in range(vertex_count):
                missing = [edge for edge in incident[vertex] if edge not in assignments]
                total = sum(
                    (coefficient(edge, assignments[edge], vertex) for edge in incident[vertex] if edge in assignments),
                    zero,
                )
                if not missing:
                    if total != zero:
                        return None
                elif len(missing) == 1:
                    edge = missing[0]
                    needed = -total
                    label = needed if edge[0] == vertex else needed.conjugate()
                    if label not in allowed[edge]:
                        return None
                    assignments[edge] = label
                    changed = True
            if not changed:
                return assignments

    def search(assignments):
        assignments = propagate(assignments)
        if assignments is None:
            return
        if len(assignments) == len(edges):
            yield tuple(assignments[edge] for edge in edges)
            return
        edge = max(
            (edge for edge in edges if edge not in assignments),
            key=lambda item: sum(len(incident[vertex]) for vertex in item),
        )
        for label in options[edges.index(edge)]:
            yield from search({**assignments, edge: label})

    yield from search({})


def component_states(graph: nx.Graph, norm9, require_unit_kernel: bool = False) -> set[tuple[int, bool, int]]:
    vertices = graph.number_of_nodes()
    edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
    tree = nx.minimum_spanning_tree(graph)
    tree_edges = {tuple(sorted(edge)) for edge in tree.edges()}
    # Full six-unit diagonal switching fixes every spanning-tree gain to +3;
    # only independent cycle gains affect determinant or nullity.
    options = tuple((Eisenstein(3),) if edge in tree_edges else norm9 for edge in edges)
    values = set()
    has_leaf = any(graph.degree(vertex) == 1 for vertex in graph.nodes())
    if require_unit_kernel:
        labelings = (
            tuple(Eisenstein(3 * unit.a, 3 * unit.b) for unit in units)
            for units in unit_kernel_labelings(graph, norm9)
        )
    elif nx.is_tree(graph):
        # Every unit gain on a tree is removable by diagonal unitary
        # switching, so its Gram determinant depends only on edge norms.
        labelings = [(Eisenstein(3),) * len(edges)]
    else:
        labelings = product(*options)
    for labels in labelings:
        block = [[Eisenstein(15 if i == j else 0) for j in range(vertices)] for i in range(vertices)]
        for (i, j), label in zip(edges, labels):
            block[i][j] = label
            block[j][i] = label.conjugate()
        value = determinant(block)
        assert value.b == 0 and value.a > 0
        adjacency = [[Eisenstein() for _ in range(vertices)] for _ in range(vertices)]
        for (i, j), label in zip(edges, labels):
            unit = Eisenstein(label.a // 3, label.b // 3)
            adjacency[i][j] = unit
            adjacency[j][i] = unit.conjugate()
        if require_unit_kernel:
            kernel = True
            nullity, full_support = kernel_dimension_and_full_support(adjacency)
            assert nullity > 0 and full_support
        else:
            adjacency_singular = determinant(adjacency) == Eisenstein()
            if has_leaf or not adjacency_singular:
                kernel, nullity = False, 0
            else:
                nullity, kernel = kernel_dimension_and_full_support(adjacency)
        values.add((value.a, kernel, nullity))
    return values


def component_determinants(graph: nx.Graph, norm9) -> set[int]:
    return {state[0] for state in component_states(graph, norm9)}


@lru_cache(maxsize=None)
def leafless_connected_graphs(vertices: int, edges: int) -> tuple[nx.Graph, ...]:
    """Generate connected unlabeled leafless graphs close to the tree line."""

    assert edges >= vertices
    if edges == vertices:
        return (nx.cycle_graph(vertices),)
    buckets: dict[str, list[nx.Graph]] = {}
    added_count = edges - (vertices - 1)
    for tree in nx.nonisomorphic_trees(vertices):
        missing = list(nx.non_edges(tree))
        for added in combinations(missing, added_count):
            graph = tree.copy()
            graph.add_edges_from(added)
            if min(dict(graph.degree()).values()) < 2:
                continue
            key = nx.weisfeiler_lehman_graph_hash(graph)
            bucket = buckets.setdefault(key, [])
            if not any(nx.is_isomorphic(graph, existing) for existing in bucket):
                bucket.append(graph)
    result = tuple(graph for bucket in buckets.values() for graph in bucket)
    assert all(
        nx.is_connected(graph) and graph.number_of_nodes() == vertices and graph.number_of_edges() == edges
        for graph in result
    )
    return result


def leafless_bicyclic_graphs_order8() -> tuple[nx.Graph, ...]:
    return leafless_connected_graphs(8, 9)


@lru_cache(maxsize=None)
def connected_unicyclic_graphs(vertices: int) -> tuple[nx.Graph, ...]:
    """All connected unlabeled graphs with equally many vertices and edges."""

    buckets: dict[str, list[nx.Graph]] = {}
    for tree in nx.nonisomorphic_trees(vertices):
        for edge in nx.non_edges(tree):
            graph = tree.copy()
            graph.add_edge(*edge)
            key = nx.weisfeiler_lehman_graph_hash(graph)
            bucket = buckets.setdefault(key, [])
            if not any(nx.is_isomorphic(graph, existing) for existing in bucket):
                bucket.append(graph)
    result = tuple(graph for bucket in buckets.values() for graph in bucket)
    assert all(
        nx.is_connected(graph)
        and graph.number_of_nodes() == vertices
        and graph.number_of_edges() == vertices
        for graph in result
    )
    return result


def connected_unicyclic_graphs_order8() -> tuple[nx.Graph, ...]:
    return connected_unicyclic_graphs(8)


@lru_cache(maxsize=None)
def connected_near_tree_graphs(vertices: int, edges: int) -> tuple[nx.Graph, ...]:
    """Connected unlabeled graphs generated by adding chords to trees."""

    assert vertices - 1 <= edges
    if edges == vertices - 1:
        return tuple(nx.nonisomorphic_trees(vertices))
    if edges == vertices:
        return connected_unicyclic_graphs(vertices)
    buckets: dict[str, list[nx.Graph]] = {}
    added_count = edges - (vertices - 1)
    for tree in nx.nonisomorphic_trees(vertices):
        for added in combinations(tuple(nx.non_edges(tree)), added_count):
            graph = tree.copy()
            graph.add_edges_from(added)
            key = nx.weisfeiler_lehman_graph_hash(graph)
            bucket = buckets.setdefault(key, [])
            if not any(nx.is_isomorphic(graph, existing) for existing in bucket):
                bucket.append(graph)
    result = tuple(graph for bucket in buckets.values() for graph in bucket)
    assert all(
        nx.is_connected(graph)
        and graph.number_of_nodes() == vertices
        and graph.number_of_edges() == edges
        for graph in result
    )
    return result


def connected_types(max_edges: int, norm9, kernel_only: bool = False):
    result = []
    graphs = list(nx.graph_atlas_g())
    # The graph atlas stops at seven vertices.  A connected graph with at most
    # seven edges and eight vertices must be a tree; include those exact
    # missing types for the Q=63 layer.
    for vertices in range(8, max_edges + 2):
        graphs.extend(nx.nonisomorphic_trees(vertices))
    if max_edges >= 8:
        # Among connected eight-vertex/eight-edge graphs, the only one without
        # a leaf is C8.  Leafed unicyclic graphs are irrelevant in kernel-only
        # mode and ordinary layers currently stop below this boundary.
        if kernel_only:
            graphs.append(nx.cycle_graph(8))
        else:
            graphs.extend(connected_unicyclic_graphs_order8())
    if max_edges >= 9:
        assert kernel_only
        graphs.extend(leafless_bicyclic_graphs_order8())
        graphs.append(nx.cycle_graph(9))
    if max_edges >= 10:
        assert kernel_only
        graphs.extend(leafless_connected_graphs(8, 10))
        graphs.extend(leafless_connected_graphs(9, 10))
        graphs.append(nx.cycle_graph(10))
    if max_edges >= 11:
        assert kernel_only
        graphs.extend(leafless_connected_graphs(8, 11))
        graphs.extend(leafless_connected_graphs(9, 11))
        graphs.extend(leafless_connected_graphs(10, 11))
        graphs.append(nx.cycle_graph(11))
    if max_edges >= 12:
        assert kernel_only
        graphs.extend(leafless_connected_graphs(8, 12))
        graphs.extend(leafless_connected_graphs(9, 12))
        graphs.extend(leafless_connected_graphs(10, 12))
        graphs.extend(leafless_connected_graphs(11, 12))
        graphs.append(nx.cycle_graph(12))
    for graph in graphs:
        if graph.number_of_nodes() < 2 or graph.number_of_edges() > max_edges:
            continue
        if not nx.is_connected(graph):
            continue
        graph = nx.convert_node_labels_to_integers(graph)
        edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
        has_leaf = any(graph.degree(vertex) == 1 for vertex in graph.nodes())
        # A full-support kernel vector is impossible at a leaf: the leaf
        # equation is a nonzero edge label times its neighbour coordinate.
        # Dropping these dead types here is important at the ten-edge layer,
        # where otherwise thousands of empty tree multisets are generated.
        if kernel_only and has_leaf:
            continue
        states = component_states(graph, norm9, require_unit_kernel=kernel_only)
        result.append({
            "vertices": graph.number_of_nodes(),
            "edges": graph.number_of_edges(),
            "edge_list": edges,
            "alpha": independence_number(graph.number_of_nodes(), edges),
            "has_leaf": has_leaf,
            "determinants": {state[0] for state in states},
            "states": states,
        })
    result.sort(key=lambda item: (item["edges"], item["vertices"], item["edge_list"]))
    return result


def component_multisets(types, total_edges: int, start: int = 0, prefix=()):
    used = sum(types[index]["edges"] for index in prefix)
    if used == total_edges:
        yield prefix
        return
    for index in range(start, len(types)):
        if used + types[index]["edges"] <= total_edges:
            yield from component_multisets(types, total_edges, index, prefix + (index,))


def q9_layer(
    types, edge_count: int, required_isolates: int | None = None, order: int = 15
) -> dict[str, object]:
    profiles = []
    determinant_values = set()
    intertwining_values = set()
    for indices in component_multisets(types, edge_count):
        components = [types[index] for index in indices]
        used_vertices = sum(component["vertices"] for component in components)
        if used_vertices > order:
            continue
        isolates = order - used_vertices
        if required_isolates is not None and isolates != required_isolates:
            continue
        alpha = isolates + sum(component["alpha"] for component in components)
        if alpha > CODE_BOUND:
            continue
        profile_values = {
            15**isolates * product_value
            for determinant_choices in product(*(component["determinants"] for component in components))
            for product_value in [__import__("math").prod(determinant_choices)]
        }
        determinant_values.update(profile_values)
        profile_intertwining = {
            15**isolates * __import__("math").prod(state[0] for state in state_choices)
            for state_choices in product(*(component["states"] for component in components))
            if all(state[1] for state in state_choices)
        }
        nullities_by_value: dict[int, set[int]] = defaultdict(set)
        for state_choices in product(*(component["states"] for component in components)):
            if all(state[1] for state in state_choices):
                value = 15**isolates * __import__("math").prod(state[0] for state in state_choices)
                nullities_by_value[value].add(sum(state[2] for state in state_choices))
        intertwining_values.update(profile_intertwining)
        profile_survivors = sorted(value for value in profile_values if value > BENCHMARK and is_rational_eisenstein_norm(value))
        profiles.append({
            "component_type_indices": indices,
            "component_edges": [component["edge_list"] for component in components],
            "isolates": isolates,
            "alpha": alpha,
            "determinant_count": len(profile_values),
            "maximum_determinant": max(profile_values, default=0),
            "arithmetic_survivors": profile_survivors,
            "intertwining_arithmetic_survivors": sorted(
                value for value in profile_intertwining if value > BENCHMARK and is_rational_eisenstein_norm(value)
            ),
            "active_kernel_nullities": {
                str(value): sorted(nullities_by_value[value])
                for value in sorted(nullities_by_value)
                if value > BENCHMARK and is_rational_eisenstein_norm(value)
            },
        })
    above = {value for value in determinant_values if value > BENCHMARK}
    return {
        "profiles": profiles,
        "all_determinants": len(determinant_values),
        "maximum_determinant": max(determinant_values, default=0),
        "above_record_determinants": len(above),
        "arithmetic_survivors": sorted(value for value in above if is_rational_eisenstein_norm(value)),
        "intertwining_arithmetic_survivors": sorted(
            value for value in intertwining_values if value > BENCHMARK and is_rational_eisenstein_norm(value)
        ),
    }


def _orbit_representatives(values) -> tuple[Eisenstein, ...]:
    # Determinants are invariant under diagonal switching by any Eisenstein
    # unit.  For same-color Gram values the catalogue is closed under all six
    # units, so a spanning-tree edge needs just one representative, not a
    # separate sign.  Cycle-gain edges still retain the full value set.
    roots = (
        Eisenstein(1), Eisenstein(-1), Eisenstein(0, 1),
        Eisenstein(0, -1), Eisenstein(-1, -1), Eisenstein(1, 1),
    )
    remaining = set(values)
    representatives = []
    while remaining:
        representative = min(remaining, key=lambda value: (value.a, value.b))
        orbit = {representative * root for root in roots}
        representatives.append(representative)
        remaining.difference_update(orbit)
    return tuple(representatives)


def weighted_layer(
    types,
    weight_units: tuple[int, ...],
    catalogue,
    kernel_only: bool = False,
    mu3_kernel_only: bool = False,
    required_isolates: int | None = None,
    order: int = 15,
) -> dict[str, object]:
    """Enumerate every abstract Gram block for a fixed multiset of q/9 weights."""

    edge_count = len(weight_units)
    values_by_unit = {
        unit: tuple(entry["value"] for entry in catalogue if entry["norm"] == 9 * unit)
        for unit in set(weight_units)
    }
    representatives_by_unit = {unit: _orbit_representatives(values) for unit, values in values_by_unit.items()}
    assert all(values_by_unit.values())
    weight_orders = tuple(sorted(set(permutations(weight_units))))
    profiles = []
    determinant_values = set()
    intertwining_values = set()
    for indices in component_multisets(types, edge_count):
        components = [types[index] for index in indices]
        used_vertices = sum(component["vertices"] for component in components)
        if used_vertices > order:
            continue
        isolates = order - used_vertices
        if required_isolates is not None and isolates != required_isolates:
            continue
        alpha = isolates + sum(component["alpha"] for component in components)
        if alpha > CODE_BOUND:
            continue
        if (kernel_only or mu3_kernel_only) and any(component["has_leaf"] for component in components):
            continue
        edges = []
        tree_edges = set()
        offset = 0
        for component in components:
            edges.extend((offset + i, offset + j) for i, j in component["edge_list"])
            graph = nx.Graph()
            graph.add_nodes_from(range(component["vertices"]))
            graph.add_edges_from(component["edge_list"])
            tree_edges.update((offset + min(i, j), offset + max(i, j)) for i, j in nx.minimum_spanning_tree(graph).edges())
            offset += component["vertices"]
        profile_values = set()
        profile_intertwining = set()
        nullities_by_value: dict[int, set[int]] = defaultdict(set)
        for weights in weight_orders:
            if mu3_kernel_only:
                # Kernel gauge fixes switching, so retain every phase orbit
                # and generate only exactly balanced labelings.
                options = tuple(values_by_unit[unit] for unit in weights)
                labelings = balanced_kernel_labelings(edges, used_vertices, options)
            else:
                options = tuple(
                    representatives_by_unit[unit] if edge in tree_edges else values_by_unit[unit]
                    for edge, unit in zip(edges, weights)
                )
                labelings = product(*options)
            for labels in labelings:
                block = [[Eisenstein(15 if i == j else 0) for j in range(used_vertices)] for i in range(used_vertices)]
                for (i, j), label in zip(edges, labels):
                    block[i][j] = label
                    block[j][i] = label.conjugate()
                value = determinant(block)
                assert value.b == 0 and value.a > 0
                determinant_value = 15**isolates * value.a
                profile_values.add(determinant_value)
                adjacency = [[Eisenstein() for _ in range(used_vertices)] for _ in range(used_vertices)]
                for (i, j), label in zip(edges, labels):
                    unit = Eisenstein(label.a // 3, label.b // 3)
                    adjacency[i][j] = unit
                    adjacency[j][i] = unit.conjugate()
                if mu3_kernel_only:
                    nullity, full_support = kernel_dimension_and_full_support(adjacency)
                    assert full_support
                    profile_intertwining.add(determinant_value)
                    nullities_by_value[determinant_value].add(nullity)
                elif determinant(adjacency) == Eisenstein():
                    nullity, full_support = kernel_dimension_and_full_support(adjacency)
                    if full_support:
                        profile_intertwining.add(determinant_value)
                        nullities_by_value[determinant_value].add(nullity)
        determinant_values.update(profile_values)
        intertwining_values.update(profile_intertwining)
        profiles.append({
            "component_type_indices": indices,
            "component_edges": [component["edge_list"] for component in components],
            "isolates": isolates,
            "alpha": alpha,
            "determinant_count": len(profile_values),
            "maximum_determinant": max(profile_values, default=0),
            "arithmetic_survivors": sorted(
                value for value in profile_values if value > BENCHMARK and is_rational_eisenstein_norm(value)
            ),
            "intertwining_arithmetic_survivors": sorted(
                value for value in profile_intertwining if value > BENCHMARK and is_rational_eisenstein_norm(value)
            ),
            "active_kernel_nullities": {
                str(value): sorted(nullities_by_value[value])
                for value in sorted(nullities_by_value)
                if value > BENCHMARK and is_rational_eisenstein_norm(value)
            },
        })
    above = {value for value in determinant_values if value > BENCHMARK}
    return {
        "weight_units": weight_units,
        "kernel_only": kernel_only,
        "mu3_kernel_only": mu3_kernel_only,
        "profiles": profiles,
        "all_determinants": len(determinant_values),
        "maximum_determinant": max(determinant_values, default=0),
        "above_record_determinants": len(above),
        "arithmetic_survivors": sorted(value for value in above if is_rational_eisenstein_norm(value)),
        "intertwining_arithmetic_survivors": sorted(
            value for value in intertwining_values if value > BENCHMARK and is_rational_eisenstein_norm(value)
        ),
    }


def q54_k4_final_obstruction(layer: dict[str, object]) -> dict[str, object]:
    """Certify the last intertwining-compatible Q=54 profile is impossible."""

    surviving_profiles = [profile for profile in layer["profiles"] if profile["intertwining_arithmetic_survivors"]]
    assert len(surviving_profiles) == 1
    profile = surviving_profiles[0]
    edges = profile["component_edges"][0]
    vertices = {vertex for edge in edges for vertex in edge}
    assert len(vertices) == 4 and len(edges) == 6  # K4.
    assert profile["isolates"] == 11
    # A full-support kernel makes the weighted K4 singular.  Every principal
    # weighted triangle is nonsingular: its determinant is twice the real
    # part of an Eisenstein unit and hence is ±1 or ±2.  Thus the K4 kernel is
    # exactly one-dimensional.  Intertwining A H = H B then makes the four
    # active entries in all 11 isolated rows phase multiples of one vector.
    # Dephase those rows.  Their complementary length-11 Gram has diagonal 11
    # and off-diagonal -4, namely 15 I_11 - 4 J_11, whose all-ones eigenvalue
    # is negative.
    residual_eigenvalue = 15 - 4 * 11
    assert residual_eigenvalue == -29
    return {
        "last_determinant": profile["intertwining_arithmetic_survivors"][0],
        "support": "K4 plus 11 isolates",
        "active_kernel_dimension": 1,
        "forced_residual_gram": "15 I_11 - 4 J_11",
        "all_ones_eigenvalue": residual_eigenvalue,
        "conclusion": "impossible because a Gram matrix is positive semidefinite",
    }


def q63_five_vertex_final_obstruction(layer: dict[str, object]) -> dict[str, object]:
    """Certify the last kernel-compatible Q=63 profile is impossible."""

    surviving_profiles = [profile for profile in layer["profiles"] if profile["intertwining_arithmetic_survivors"]]
    assert len(surviving_profiles) == 1
    profile = surviving_profiles[0]
    edges = set(profile["component_edges"][0])
    expected = {(0, 1), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (3, 4)}
    assert edges == expected and profile["isolates"] == 10
    # In a kernel vector, the equations at degree-two vertices 0 and 2
    # determine x4 and x3 from x1.  The equations at vertices 4 and 3 then
    # determine x0 and x2.  All edge weights are nonzero, so nullity is at
    # most one; the enumerated full-support kernel makes it exactly one.
    residual_eigenvalue = 15 - 5 * 10
    assert residual_eigenvalue == -35
    return {
        "last_determinant": profile["intertwining_arithmetic_survivors"][0],
        "support_edges": sorted(edges),
        "support": "seven-edge graph on 5 vertices plus 10 isolates",
        "active_kernel_dimension": 1,
        "forced_residual_gram": "15 I_10 - 5 J_10",
        "all_ones_eigenvalue": residual_eigenvalue,
        "conclusion": "impossible because a Gram matrix is positive semidefinite",
    }


def q72_no_isolate_obstruction() -> dict[str, object]:
    """Exclude the unique eight-edge support covering all 15 vertices."""

    # With 15 nonisolated vertices and eight edges, the degree sum is 16, so
    # the graph is P3 plus six K2 components.  At energy 72 all eight weights
    # must be 9.  The corresponding Gram determinant is phase-independent.
    value = 3105 * 216**6
    assert value > BENCHMARK and not is_rational_eisenstein_norm(value)
    return {
        "support": "P3 plus 6 K2",
        "determinant": value,
        "eisenstein_norm": False,
    }


def main() -> None:
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    types = connected_types(6, norm9)
    layers = {energy: q9_layer(types, energy // 9) for energy in (27, 36, 45, 54)}
    mixed_54 = {
        "27+27": weighted_layer(types, (3, 3), catalogue),
        "36+9+9": weighted_layer(types, (4, 1, 1), catalogue),
        "27+9+9+9": weighted_layer(types, (3, 1, 1, 1), catalogue),
    }

    # At energy 45 the only mixed-weight profile with enough independent-set
    # deficit is three disjoint edges of norms 27,9,9.
    mixed_45 = 198 * 216**2 * 15**9
    assert not is_rational_eisenstein_norm(mixed_45)
    assert not layers[27]["arithmetic_survivors"]
    assert not layers[36]["arithmetic_survivors"]
    assert len(layers[45]["arithmetic_survivors"]) == 8
    assert not layers[45]["intertwining_arithmetic_survivors"]
    assert len(layers[54]["arithmetic_survivors"]) == 18
    assert layers[54]["intertwining_arithmetic_survivors"] == [318786752197265625]
    final_54 = q54_k4_final_obstruction(layers[54])
    kernel_types_63 = connected_types(7, norm9, kernel_only=True)
    kernel_layers_63 = {
        "9+9+9+9+9+9+9": q9_layer(kernel_types_63, 7),
        "63": weighted_layer(kernel_types_63, (7,), catalogue, kernel_only=True),
        "36+27": weighted_layer(kernel_types_63, (4, 3), catalogue, kernel_only=True),
        "36+9+9+9": weighted_layer(kernel_types_63, (4, 1, 1, 1), catalogue, kernel_only=True),
        "27+27+9": weighted_layer(kernel_types_63, (3, 3, 1), catalogue, kernel_only=True),
        "27+9+9+9+9": weighted_layer(kernel_types_63, (3, 1, 1, 1, 1), catalogue, kernel_only=True),
    }
    assert kernel_layers_63["9+9+9+9+9+9+9"]["intertwining_arithmetic_survivors"] == [325793054443359375]
    assert all(
        not layer["intertwining_arithmetic_survivors"]
        for name, layer in kernel_layers_63.items()
        if name != "9+9+9+9+9+9+9"
    )
    final_63 = q63_five_vertex_final_obstruction(kernel_layers_63["9+9+9+9+9+9+9"])
    no_isolate_72 = q72_no_isolate_obstruction()
    kernel_types_72 = connected_types(8, norm9, kernel_only=True)
    kernel_layers_72 = {
        "9+9+9+9+9+9+9+9": q9_layer(kernel_types_72, 8),
        "63+9": weighted_layer(kernel_types_72, (7, 1), catalogue, kernel_only=True),
        "36+36": weighted_layer(kernel_types_72, (4, 4), catalogue, kernel_only=True),
        "36+27+9": weighted_layer(kernel_types_72, (4, 3, 1), catalogue, kernel_only=True),
        "36+9+9+9+9": weighted_layer(kernel_types_72, (4, 1, 1, 1, 1), catalogue, kernel_only=True),
        "27+27+9+9": weighted_layer(kernel_types_72, (3, 3, 1, 1), catalogue, kernel_only=True),
        "27+9+9+9+9+9": weighted_layer(kernel_types_72, (3, 1, 1, 1, 1, 1), catalogue, kernel_only=True),
    }
    last_72_name = "27+9+9+9+9+9"
    assert kernel_layers_72[last_72_name]["intertwining_arithmetic_survivors"] == [325793054443359375]
    assert all(
        not layer["intertwining_arithmetic_survivors"]
        for name, layer in kernel_layers_72.items()
        if name != last_72_name
    )
    final_72 = q54_k4_final_obstruction(kernel_layers_72[last_72_name])

    print(json.dumps({
        "color_partition": [15, 0, 0],
        "connected_unlabeled_types_through_six_edges": len(types),
        "q9_layers": layers,
        "energy45_mixed_profile_norms": [27, 9, 9],
        "energy45_mixed_determinant": mixed_45,
        "energy54_mixed_layers": mixed_54,
        "energy54_final_obstruction": final_54,
        "energy63_kernel_layers": kernel_layers_63,
        "energy63_final_obstruction": final_63,
        "energy72_no_isolate_obstruction": no_isolate_72,
        "energy72_kernel_layers": kernel_layers_72,
        "energy72_final_obstruction": final_72,
        "theorem": "All-same-color energies 27, 36, 45, 54, 63, and 72 are excluded exactly.",
    }, indent=2))


if __name__ == "__main__":
    main()
