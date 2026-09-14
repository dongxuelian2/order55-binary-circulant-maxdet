"""Consolidated exact audit that every energy below Q=144 is below its envelope."""

from __future__ import annotations

import json
from fractions import Fraction
from math import isqrt

import networkx as nx

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy90 import weight_partitions
    from order15_mu3.scripts.verify_color_15_low_energy import connected_types, q9_layer, weighted_layer
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q114_internal_support import EXPECTED_MAXIMA
    from order15_mu3.scripts.verify_q126_boundary import capped_trace_candidates
    from order15_mu3.scripts.verify_q126_k6_minus_edge import certificate as q126_k6_minus_edge_certificate
    from order15_mu3.scripts.verify_q126_isolated_k4 import certificate as q126_isolated_k4_certificate
    from order15_mu3.scripts.verify_q126_k4_four_forest import certificate as q126_k4_four_forest_certificate
    from order15_mu3.scripts.verify_q126_k4_three_forest import certificate as q126_k4_three_forest_certificate
    from order15_mu3.scripts.verify_q126_k4_leaf import certificate as q126_k4_leaf_certificate
    from order15_mu3.scripts.verify_q126_pure_k4_six_forest import (
        certificate as q126_pure_k4_six_forest_certificate,
    )
    from order15_mu3.scripts.verify_q126_pure_k5 import certificate as q126_pure_k5_certificate
    from order15_mu3.scripts.verify_q126_triangle_components import certificate as q126_triangle_certificate
    from order15_mu3.scripts.verify_q126_triangle_active6 import certificate as q126_triangle_active6_certificate
    from order15_mu3.scripts.verify_q126_triangle_active5 import certificate as q126_triangle_active5_certificate
    from order15_mu3.scripts.verify_q126_triangle_active4 import certificate as q126_triangle_active4_certificate
    from order15_mu3.scripts.verify_q126_triangle_free_moment import (
        NO_K4_PURE_UPPER,
        UPPER as TRIANGLE_FREE_UPPER,
    )
    from order15_mu3.scripts.verify_q135_boundary import fixed_isolate_capped_upper
    from order15_mu3.scripts.verify_q141_boundary import certificate as q141_certificate
    from order15_mu3.scripts.verify_q87_q90_top_obstruction import (
        certificate as energy9_top_obstruction,
    )
    from order15_mu3.scripts.verify_q87_q90_second_obstruction import (
        NEXT as ENERGY9_REFINED_UPPER,
        certificate as energy9_second_obstruction,
    )
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:
    from verify_color_15_energy90 import weight_partitions
    from verify_color_15_low_energy import connected_types, q9_layer, weighted_layer
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q114_internal_support import EXPECTED_MAXIMA
    from verify_q126_boundary import capped_trace_candidates
    from verify_q126_k6_minus_edge import certificate as q126_k6_minus_edge_certificate
    from verify_q126_isolated_k4 import certificate as q126_isolated_k4_certificate
    from verify_q126_k4_four_forest import certificate as q126_k4_four_forest_certificate
    from verify_q126_k4_three_forest import certificate as q126_k4_three_forest_certificate
    from verify_q126_k4_leaf import certificate as q126_k4_leaf_certificate
    from verify_q126_pure_k4_six_forest import certificate as q126_pure_k4_six_forest_certificate
    from verify_q126_pure_k5 import certificate as q126_pure_k5_certificate
    from verify_q126_triangle_components import certificate as q126_triangle_certificate
    from verify_q126_triangle_active6 import certificate as q126_triangle_active6_certificate
    from verify_q126_triangle_active5 import certificate as q126_triangle_active5_certificate
    from verify_q126_triangle_active4 import certificate as q126_triangle_active4_certificate
    from verify_q126_triangle_free_moment import NO_K4_PURE_UPPER, UPPER as TRIANGLE_FREE_UPPER
    from verify_q135_boundary import fixed_isolate_capped_upper
    from verify_q141_boundary import certificate as q141_certificate
    from verify_q87_q90_top_obstruction import certificate as energy9_top_obstruction
    from verify_q87_q90_second_obstruction import (
        NEXT as ENERGY9_REFINED_UPPER,
        certificate as energy9_second_obstruction,
    )
    from verify_trace_stability import stationary_upper


def _capped_upper(isolates: int, cap: Fraction, variance_energy: int = 252) -> Fraction:
    if isolates:
        return fixed_isolate_capped_upper(isolates, variance_energy, cap)
    return max(row[3] for row in capped_trace_candidates(cap, variance_energy))


def _sqrt_upper(value: Fraction, denominator: int = 100) -> Fraction:
    scaled = value.numerator * denominator**2
    root = isqrt(scaled // value.denominator)
    while root**2 * value.denominator < scaled:
        root += 1
    return Fraction(root, denominator)


def _determinant(matrix: list[list[Fraction]]) -> Fraction:
    rows = [row[:] for row in matrix]
    value = Fraction(1)
    for column in range(len(rows)):
        pivot = next(row for row in range(column, len(rows)) if rows[row][column])
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            value = -value
        diagonal = rows[column][column]
        value *= diagonal
        for row in range(column + 1, len(rows)):
            factor = rows[row][column] / diagonal
            for index in range(column, len(rows)):
                rows[row][index] -= factor * rows[column][index]
    return value


def _weighted_k5_two_leaf_cap() -> Fraction:
    """Certify lambda_max<28.3 for the sole tight weighted K5 topology."""

    # Perron edge-moving puts both new leaves at one K5 vertex and the sole
    # sqrt(3)-weighted edge incident with that vertex.  Replace sqrt(3) by
    # the strict rational upper 26/15.  Sylvester's criterion then proves the
    # comparison adjacency has spectral radius below 133/30.
    assert Fraction(26, 15) ** 2 > 3
    edges = [(i, j) for i in range(5) for j in range(i + 1, 5)] + [(0, 5), (0, 6)]
    adjacency = [[Fraction(0) for _ in range(7)] for _ in range(7)]
    for i, j in edges:
        weight = Fraction(26, 15) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    threshold = Fraction(133, 30)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(7)]
        for i in range(7)
    ]
    assert all(_determinant([row[:size] for row in comparison[:size]]) > 0 for size in range(1, 8))
    return Fraction(283, 10)


def _weighted_k5_one_leaf_cap() -> Fraction:
    """Certify lambda_max<28.6 for weighted K5 with one new leaf."""

    # The eleven-edge partition has one weight-4 edge.  When every edge is
    # in the K5 component, the support is exactly K5 with one leaf.  Perron
    # edge-moving places the weight-2 edge on a clique edge incident with the
    # leaf's root.  Sylvester's criterion gives rho<68/15.
    edges = [(i, j) for i in range(5) for j in range(i + 1, 5)] + [(0, 5)]
    adjacency = [[Fraction(0) for _ in range(6)] for _ in range(6)]
    for i, j in edges:
        weight = Fraction(2) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    threshold = Fraction(68, 15)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(6)]
        for i in range(6)
    ]
    assert all(_determinant([row[:size] for row in comparison[:size]]) > 0 for size in range(1, 7))
    return Fraction(143, 5)


def _weighted_k5_two_heavy_edges_cap() -> Fraction:
    """Certify lambda_max<28.95 for K5 with two norm-3 edges."""

    # There are only two placements up to K5 symmetry.  Adjacent heavy edges
    # dominate disjoint ones by Perron edge-moving, so use the adjacent case.
    assert Fraction(26, 15) ** 2 > 3
    adjacency = [[Fraction(0) for _ in range(5)] for _ in range(5)]
    for i in range(5):
        for j in range(i + 1, 5):
            weight = Fraction(26, 15) if (i, j) in {(0, 1), (0, 2)} else Fraction(1)
            adjacency[i][j] = adjacency[j][i] = weight
    threshold = Fraction(93, 20)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(5)]
        for i in range(5)
    ]
    assert all(_determinant([row[:size] for row in comparison[:size]]) > 0 for size in range(1, 6))
    return Fraction(579, 20)


def _pure_k5_two_vertex_cap() -> Fraction:
    """Certify lambda_max<28.5 for a 14-edge K5 component on seven vertices."""

    # Four edges beyond K5 activate two vertices.  Splitting on whether the
    # new vertices are adjacent and applying Perron edge-moving leaves the
    # extremal pattern: one new vertex has three K5 neighbours and the other
    # is a leaf at one of those neighbours.
    edges = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    edges.extend([(0, 5), (1, 5), (2, 5), (0, 6)])
    adjacency = [[Fraction(0) for _ in range(7)] for _ in range(7)]
    for i, j in edges:
        adjacency[i][j] = adjacency[j][i] = Fraction(1)
    threshold = Fraction(9, 2)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(7)]
        for i in range(7)
    ]
    assert all(_determinant([row[:size] for row in comparison[:size]]) > 0 for size in range(1, 8))
    return Fraction(57, 2)


def _weighted_k4_six_leaf_cap() -> Fraction:
    """Certify lambda_max<27 for the tight ten-vertex weighted K4 branch."""

    # Six edges beyond K4 add exactly six vertices, so they form a rooted
    # forest.  Perron edge-moving makes them six leaves at one K4 vertex and
    # moves the sqrt(3) edge to an incident clique edge.  Replacing sqrt(3)
    # by 26/15, Sylvester's criterion proves rho<4, hence lambda_max<27.
    assert Fraction(26, 15) ** 2 > 3
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    edges.extend((0, vertex) for vertex in range(4, 10))
    adjacency = [[Fraction(0) for _ in range(10)] for _ in range(10)]
    for i, j in edges:
        weight = Fraction(26, 15) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    comparison = [
        [Fraction(4) * (i == j) - adjacency[i][j] for j in range(10)]
        for i in range(10)
    ]
    assert all(_determinant([row[:size] for row in comparison[:size]]) > 0 for size in range(1, 11))
    return Fraction(27)


def _weighted_k4_five_leaf_cap() -> Fraction:
    """Certify lambda_max<27.25 for K4 with five leaves and one weight-4 edge."""

    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    edges.extend((0, vertex) for vertex in range(4, 9))
    adjacency = [[Fraction(0) for _ in range(9)] for _ in range(9)]
    for i, j in edges:
        weight = Fraction(2) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    threshold = Fraction(49, 12)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(9)]
        for i in range(9)
    ]
    assert all(_determinant([row[:size] for row in comparison[:size]]) > 0 for size in range(1, 10))
    return Fraction(109, 4)


def _weighted_k4_five_vertex_cap() -> Fraction:
    """Certify lambda_max<27.25 for the tight nine-vertex K4 branch."""

    # Six edges beyond K4 add five vertices.  Perron edge-moving first makes
    # four of them leaves at a maximal K4 vertex; the remaining edge joins
    # the fifth vertex to the two largest K4 coordinates.  Moving the sole
    # sqrt(3)-weighted edge to the largest incident clique edge dominates
    # every placement.  The displayed comparison is therefore extremal.
    assert Fraction(26, 15) ** 2 > 3
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    edges.extend([(0, 4), (1, 4)])
    edges.extend((0, vertex) for vertex in range(5, 9))
    adjacency = [[Fraction(0) for _ in range(9)] for _ in range(9)]
    for i, j in edges:
        weight = Fraction(26, 15) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    threshold = Fraction(49, 12)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(9)]
        for i in range(9)
    ]
    assert all(_determinant([row[:size] for row in comparison[:size]]) > 0 for size in range(1, 10))
    return Fraction(109, 4)


def _weighted_k4_four_vertex_cap() -> Fraction:
    """Certify lambda_max<27 for K4 plus five edges on four vertices."""

    assert Fraction(26, 15) ** 2 > 3
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    edges.extend([(0, 4), (1, 4)])
    edges.extend((0, vertex) for vertex in range(5, 8))
    adjacency = [[Fraction(0) for _ in range(8)] for _ in range(8)]
    for i, j in edges:
        weight = Fraction(26, 15) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    comparison = [
        [Fraction(4) * (i == j) - adjacency[i][j] for j in range(8)]
        for i in range(8)
    ]
    assert all(_determinant([row[:size] for row in comparison[:size]]) > 0 for size in range(1, 9))
    return Fraction(27)


def _seven_isolate_k4_cap(kind: str) -> Fraction:
    """Perron/Sylvester caps for the three weighted active-order-eight rows."""

    assert kind in {"two_sqrt3_forest", "weight4_bicyclic", "sqrt3_tricyclic"}
    assert Fraction(26, 15) ** 2 > 3
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    if kind == "two_sqrt3_forest":
        edges.extend((0, vertex) for vertex in range(4, 8))
        heavy = {(0, 1), (0, 2)}
        heavy_weight = Fraction(26, 15)
        threshold = Fraction(83, 20)
    elif kind == "weight4_bicyclic":
        edges.extend([(0, 4), (1, 4)])
        edges.extend((0, vertex) for vertex in range(5, 8))
        heavy = {(0, 1)}
        heavy_weight = Fraction(2)
        threshold = Fraction(83, 20)
    else:
        edges.extend([(0, 4), (1, 4), (2, 4)])
        edges.extend((0, vertex) for vertex in range(5, 8))
        heavy = {(0, 1)}
        heavy_weight = Fraction(26, 15)
        threshold = Fraction(127, 30)
    adjacency = [[Fraction(0) for _ in range(8)] for _ in range(8)]
    for i, j in edges:
        adjacency[i][j] = adjacency[j][i] = heavy_weight if (i, j) in heavy else Fraction(1)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(8)]
        for i in range(8)
    ]
    assert all(
        _determinant([row[:size] for row in comparison[:size]]) > 0
        for size in range(1, 9)
    )
    return Fraction(15) + 3 * threshold


def _dense_pure_clique_cap(kind: str) -> Fraction:
    """Perron/Sylvester caps for the remaining dense pure clique rows."""

    assert kind in {"k4_order6", "k4_order8", "k5_order6", "k5_order8"}
    if kind == "k4_order6":
        order, clique, threshold = 6, 4, Fraction(22, 5)
        # The support itself is K6 minus two disjoint edges.
        graph = nx.complete_graph(6)
        graph.remove_edges_from([(0, 1), (2, 3)])
        edges = list(graph.edges())
    elif kind == "k4_order8":
        order, clique, threshold = 8, 4, Fraction(13, 3)
        edges = [(i, j) for i in range(clique) for j in range(i + 1, clique)]
        edges.extend([(0, 4), (1, 4), (2, 4), (0, 5), (1, 5), (2, 5), (0, 6), (0, 7)])
    elif kind == "k5_order6":
        order, clique, threshold = 6, 5, Fraction(133, 30)
        edges = [(i, j) for i in range(clique) for j in range(i + 1, clique)]
        edges.extend([(0, 5), (1, 5), (2, 5)])
    else:
        order, clique, threshold = 8, 5, Fraction(13, 3)
        edges = [(i, j) for i in range(clique) for j in range(i + 1, clique)]
        edges.extend([(0, 5), (1, 5), (0, 6), (0, 7)])
    adjacency = [[Fraction(0) for _ in range(order)] for _ in range(order)]
    for i, j in edges:
        adjacency[i][j] = adjacency[j][i] = Fraction(1)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(order)]
        for i in range(order)
    ]
    assert all(
        _determinant([row[:size] for row in comparison[:size]]) > 0
        for size in range(1, order + 1)
    )
    return Fraction(15) + 3 * threshold


def _pure_k4_sparse_cap(joined: int, active_vertices: int) -> Fraction:
    """Certify lambda_max<26 for the three tight pure K4 components."""

    assert (joined, active_vertices) in {(3, 5), (4, 7), (5, 9)}
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    if joined == 3:
        edges.extend([(0, 4), (1, 4), (2, 4)])
    elif joined == 4:
        edges.extend([(0, 4), (1, 4), (0, 5), (0, 6)])
    else:
        edges.extend((0, vertex) for vertex in range(4, 9))
    assert len(edges) == 6 + joined
    adjacency = [[Fraction(0) for _ in range(active_vertices)] for _ in range(active_vertices)]
    for i, j in edges:
        adjacency[i][j] = adjacency[j][i] = Fraction(1)
    threshold = Fraction(11, 3)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(active_vertices)]
        for i in range(active_vertices)
    ]
    assert all(
        _determinant([row[:size] for row in comparison[:size]]) > 0
        for size in range(1, active_vertices + 1)
    )
    return Fraction(26)


def _pure_k4_two_isolate_cap(joined: int, active_vertices: int) -> Fraction:
    """Sharpen the two endpoint-tight pure K4 rows with two isolates."""

    assert (joined, active_vertices) in {(6, 9), (7, 11)}
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    if joined == 6:
        # One added vertex has two clique neighbours; the other four edges
        # are moved to leaves at a Perron-maximal clique vertex.
        edges.extend([(0, 4), (1, 4)])
        edges.extend((0, vertex) for vertex in range(5, 9))
        threshold = Fraction(11, 3)
    else:
        # Seven forest edges are moved to leaves at one clique vertex.
        edges.extend((0, vertex) for vertex in range(4, 11))
        threshold = Fraction(37, 10)
    adjacency = [
        [Fraction(0) for _ in range(active_vertices)]
        for _ in range(active_vertices)
    ]
    for i, j in edges:
        adjacency[i][j] = adjacency[j][i] = Fraction(1)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(active_vertices)]
        for i in range(active_vertices)
    ]
    assert all(
        _determinant([row[:size] for row in comparison[:size]]) > 0
        for size in range(1, active_vertices + 1)
    )
    return Fraction(15) + 3 * threshold


def _pure_k4_three_isolate_cap(joined: int, active_vertices: int) -> Fraction:
    """Certify lambda_max<27.3 for the three sparse endpoint-tight rows."""

    assert (joined, active_vertices) in {(6, 8), (7, 10), (8, 12)}
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    if joined == 6:
        edges.extend([(0, 4), (1, 4), (2, 4)])
        edges.extend((0, vertex) for vertex in range(5, 8))
    elif joined == 7:
        edges.extend([(0, 4), (1, 4)])
        edges.extend((0, vertex) for vertex in range(5, 10))
    else:
        edges.extend((0, vertex) for vertex in range(4, 12))
    assert len(edges) == 6 + joined
    threshold = Fraction(41, 10)
    adjacency = [
        [Fraction(0) for _ in range(active_vertices)]
        for _ in range(active_vertices)
    ]
    for i, j in edges:
        adjacency[i][j] = adjacency[j][i] = Fraction(1)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(active_vertices)]
        for i in range(active_vertices)
    ]
    assert all(
        _determinant([row[:size] for row in comparison[:size]]) > 0
        for size in range(1, active_vertices + 1)
    )
    return Fraction(273, 10)


def _weighted_triangle_sparse_cap(active_vertices: int) -> Fraction:
    """Cap the 11- and 12-vertex triangle components with twelve edges."""

    assert active_vertices in (11, 12)
    assert Fraction(26, 15) ** 2 > 3
    edges = [(0, 1), (0, 2), (1, 2)]
    first_leaf = 3
    if active_vertices == 11:
        # The bicyclic extremal is a diamond with all remaining vertices
        # moved to leaves at an endpoint of its shared edge.
        edges.extend([(0, 3), (1, 3)])
        first_leaf = 4
        threshold = Fraction(47, 12)
        cap = Fraction(107, 4)
    else:
        # The unicyclic extremal is a triangle with all remaining vertices
        # moved to leaves at one endpoint of the heavy triangle edge.
        threshold = Fraction(23, 6)
        cap = Fraction(53, 2)
    edges.extend((0, vertex) for vertex in range(first_leaf, active_vertices))
    assert len(edges) == 12
    adjacency = [[Fraction(0) for _ in range(active_vertices)] for _ in range(active_vertices)]
    for i, j in edges:
        weight = Fraction(26, 15) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(active_vertices)]
        for i in range(active_vertices)
    ]
    assert all(
        _determinant([row[:size] for row in comparison[:size]]) > 0
        for size in range(1, active_vertices + 1)
    )
    return cap


def _weighted_triangle_weight4_unicyclic_cap(active_vertices: int) -> Fraction:
    """Cap a weight-4 unicyclic triangle at active orders eight through eleven."""

    assert active_vertices in (8, 9, 10, 11)
    # Perron edge-moving puts the rooted forest at one triangle vertex and
    # the weight-2 edge on an incident triangle edge.  The rational thresholds
    # below give Gram caps 15+3*rho of 25.3, 25.7, 26.1, and 26.5.
    edges = [(0, 1), (0, 2), (1, 2)]
    edges.extend((0, vertex) for vertex in range(3, active_vertices))
    adjacency = [
        [Fraction(0) for _ in range(active_vertices)]
        for _ in range(active_vertices)
    ]
    for i, j in edges:
        weight = Fraction(2) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    thresholds = {
        8: Fraction(103, 30),
        9: Fraction(107, 30),
        10: Fraction(37, 10),
        11: Fraction(23, 6),
    }
    threshold = thresholds[active_vertices]
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(active_vertices)]
        for i in range(active_vertices)
    ]
    assert all(
        _determinant([row[:size] for row in comparison[:size]]) > 0
        for size in range(1, active_vertices + 1)
    )
    return Fraction(15) + 3 * threshold


def _weighted_triangle_eleven_edge_cap(active_vertices: int) -> Fraction:
    """Certify lambda_max<26.25 for sparse 11-edge triangle components."""

    assert active_vertices in (10, 11)
    assert Fraction(26, 15) ** 2 > 3
    edges = [(0, 1), (0, 2), (1, 2)]
    first_leaf = 3
    if active_vertices == 10:
        # Bicyclic extremal: a diamond and six leaves at a shared-edge root.
        edges.extend([(0, 3), (1, 3)])
        first_leaf = 4
    # Unicyclic extremal: a triangle and eight leaves at a heavy-edge root.
    edges.extend((0, vertex) for vertex in range(first_leaf, active_vertices))
    assert len(edges) == 11
    adjacency = [[Fraction(0) for _ in range(active_vertices)] for _ in range(active_vertices)]
    for i, j in edges:
        weight = Fraction(26, 15) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    threshold = Fraction(15, 4)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(active_vertices)]
        for i in range(active_vertices)
    ]
    assert all(
        _determinant([row[:size] for row in comparison[:size]]) > 0
        for size in range(1, active_vertices + 1)
    )
    return Fraction(105, 4)


def _weighted_triangle_ten_edge_cap(active_vertices: int) -> Fraction:
    """Certify lambda_max<26 for sparse 10-edge triangle components."""

    assert active_vertices in (9, 10)
    assert Fraction(26, 15) ** 2 > 3
    edges = [(0, 1), (0, 2), (1, 2)]
    first_leaf = 3
    if active_vertices == 9:
        edges.extend([(0, 3), (1, 3)])
        first_leaf = 4
    edges.extend((0, vertex) for vertex in range(first_leaf, active_vertices))
    assert len(edges) == 10
    adjacency = [[Fraction(0) for _ in range(active_vertices)] for _ in range(active_vertices)]
    for i, j in edges:
        weight = Fraction(26, 15) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    threshold = Fraction(11, 3) if active_vertices == 9 else Fraction(53, 15)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(active_vertices)]
        for i in range(active_vertices)
    ]
    assert all(
        _determinant([row[:size] for row in comparison[:size]]) > 0
        for size in range(1, active_vertices + 1)
    )
    return Fraction(15) + 3 * threshold


def _weighted_triangle_seven_edge_cap() -> Fraction:
    """Certify lambda_max<25 for the tight weighted seven-edge component."""

    edges = [(0, 1), (0, 2), (1, 2)] + [(0, vertex) for vertex in range(3, 7)]
    adjacency = [[Fraction(0) for _ in range(7)] for _ in range(7)]
    for i, j in edges:
        weight = Fraction(2) if (i, j) == (0, 1) else Fraction(1)
        adjacency[i][j] = adjacency[j][i] = weight
    threshold = Fraction(10, 3)
    comparison = [
        [threshold * (i == j) - adjacency[i][j] for j in range(7)]
        for i in range(7)
    ]
    assert all(_determinant([row[:size] for row in comparison[:size]]) > 0 for size in range(1, 8))
    return Fraction(25)


def certificate() -> dict[str, object]:
    q144 = max(stationary_upper(288, multiplicity) for multiplicity in range(1, 15))

    # The Q=141 certificate uses bounds uniform in total energy for every
    # size-13 internal energy 9,...,63; its additional e=45,c=78 row covers
    # the only smaller-total allocation not literally present at Q=141.
    color = q141_certificate()
    assert color["size13_max_over_q144"] < 1
    assert color["size13_energy45_minimal_cross_over_q144"] < 1
    assert energy9_top_obstruction()["next_arithmetic_energy9_upper"] == 289220156718750000
    assert energy9_second_obstruction()["next_arithmetic_energy9_upper"] == ENERGY9_REFINED_UPPER

    # No-isolate (14,1) supports have internal energy 63,...,99.  For e=63
    # the support is a matching.  For e=72 it is P4+5K2 or 2P3+4K2, so the
    # adjacency spectral radius is <5 and lambda_max<20.  The later exact
    # internal layers use safe weighted clique caps.
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    e72_layer = q9_layer(connected_types(3, norm9), 8, required_isolates=0, order=14)
    assert e72_layer["maximum_determinant"] <= EXPECTED_MAXIMA[72]
    e81_unit = q9_layer(connected_types(5, norm9), 9, required_isolates=0, order=14)
    e90_weighted = weighted_layer(
        connected_types(8, norm9), (3,) + (1,) * 7,
        catalogue, required_isolates=0, order=14,
    )
    e90_unit = q9_layer(connected_types(7, norm9), 10, required_isolates=0, order=14)
    internal = {
        63: 216**7,
        72: e72_layer["maximum_determinant"],
        81: max(198 * 216**6, e81_unit["maximum_determinant"]),
        90: max(e90_weighted["maximum_determinant"], e90_unit["maximum_determinant"]),
        99: max(color["size14_internal99_maxima"].values()),
    }
    lambda_caps = {
        63: Fraction(18),
        72: Fraction(20),
        81: Fraction(261, 10),
        90: Fraction(27),
        99: Fraction(138, 5),
    }
    size14 = {
        energy: Fraction(value) * (Fraction(15) - Fraction(42) / lambda_caps[energy])
        for energy, value in internal.items()
    }
    assert all(value < q144 for value in size14.values())

    # At all-same Q=126, first suppose the support has no K5.  A component
    # without K4 has rho^2<=168 and hence lambda_max<28.  Otherwise choose a
    # spectral-maximizing component containing K4.  If the whole support has
    # m edges and r of the edges beyond the K4 remain in that component, then
    # its energy is at most 126-9(m-6-r), while the full support has at least
    # max(0,23+r-2m) isolates.  Weighted Motzkin--Straus gives rho^2<=3Q_c/2.
    supported_edges = sorted({len(weights) for weights in weight_partitions(14) if len(weights) >= 6})
    assert supported_edges == [6, 7, 8, 9, 10, 11, 12, 14]
    triangle_free = TRIANGLE_FREE_UPPER
    triangle_rows = []
    supported_all = sorted({len(weights) for weights in weight_partitions(14)})
    exact_triangle = q126_triangle_certificate()
    exact_triangle_active6 = q126_triangle_active6_certificate()
    exact_triangle_active5 = q126_triangle_active5_certificate()
    exact_triangle_active4 = q126_triangle_active4_certificate()
    for edges in (value for value in supported_all if value >= 3):
        for joined in range(edges - 2):
            outside_edges = edges - 3 - joined
            if edges == 14:
                component_edges = 3 + joined
                vertices = (3,) if joined == 0 else range(4, 4 + joined)
                for active_vertices in vertices:
                    if active_vertices + 2 * outside_edges > 15:
                        continue
                    parts = [active_vertices // 3 + (index < active_vertices % 3) for index in range(3)]
                    turan_edges = (active_vertices**2 - sum(part**2 for part in parts)) // 2
                    if component_edges > turan_edges:
                        continue
                    isolates = max(0, 15 - active_vertices - 2 * outside_edges)
                    cap = Fraction(15) + _sqrt_upper(
                        Fraction(9 * (2 * component_edges - active_vertices + 1))
                    )
                    # This branch has clique number at most three, so the
                    # weighted Motzkin--Straus inequality also gives
                    # rho(E)^2<=168 and hence lambda_max<28.
                    cap = min(cap, Fraction(28))
                    upper = NO_K4_PURE_UPPER
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "upper": upper,
                    })
                continue
            component_energy = 126 - 9 * outside_edges
            if edges == 9 and joined == 1:
                upper = Fraction(exact_triangle_active4["maxima"]["edges9_component4"])
                triangle_rows.append({
                    "edges": edges,
                    "joined_to_triangle": joined,
                    "active_vertices": 4,
                    "isolates": 1,
                    "lambda_cap": Fraction(127, 5),
                    "exact_determinant": True,
                    "upper": upper,
                })
                continue
            if edges == 10 and joined == 2:
                for active_vertices in (4, 5):
                    isolates = 5 - active_vertices
                    cap = Fraction(127, 5)
                    upper = (
                        Fraction(exact_triangle_active5["maxima"]["edges10_component5"])
                        if active_vertices == 5
                        else Fraction(exact_triangle_active4["maxima"]["edges10_component5"])
                    )
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "exact_determinant": True,
                        "upper": upper,
                    })
                continue
            if edges == 11 and joined == 3:
                # Active order four would be K4 and belongs to the other
                # branch, so the K4-free component has exactly five vertices.
                upper = Fraction(exact_triangle_active5["maxima"]["edges11_component6"])
                triangle_rows.append({
                    "edges": edges,
                    "joined_to_triangle": joined,
                    "active_vertices": 5,
                    "isolates": 0,
                    "lambda_cap": Fraction(127, 5),
                    "exact_determinant": True,
                    "upper": upper,
                })
                continue
            if edges == 12 and joined == 4:
                upper = Fraction(exact_triangle_active5["maxima"]["edges12_component7"])
                triangle_rows.append({
                    "edges": edges,
                    "joined_to_triangle": joined,
                    "active_vertices": 5,
                    "isolates": 0,
                    "lambda_cap": Fraction(127, 5),
                    "exact_determinant": True,
                    "upper": upper,
                })
                continue
            if edges == 10 and joined == 3:
                for active_vertices in range(4, 7):
                    isolates = 7 - active_vertices
                    cap = Fraction(649, 25)
                    upper = (
                        Fraction(exact_triangle_active6["maxima"]["edges10_component6"])
                        if active_vertices == 6
                        else _capped_upper(isolates, cap)
                    )
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "exact_determinant": active_vertices == 6,
                        "upper": upper,
                    })
                continue
            if edges == 11 and joined in (5, 6, 7, 8):
                # At maximal active order the component is a weighted
                # unicyclic triangle plus rooted forest.  Smaller active
                # orders have additional isolates and retain the trace cap.
                outside_edges = 8 - joined
                maximal_active = 3 + joined
                generic_cap = Fraction(15) + _sqrt_upper(
                    Fraction(4 * component_energy, 3)
                )
                for active_vertices in range(4, maximal_active + 1):
                    isolates = 15 - active_vertices - 2 * outside_edges
                    cap = (
                        _weighted_triangle_weight4_unicyclic_cap(active_vertices)
                        if active_vertices == maximal_active
                        else generic_cap
                    )
                    upper = _capped_upper(isolates, cap)
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "upper": upper,
                    })
                continue
            if edges == 11 and joined == 4:
                # Four outside K2 components remain.  Only active order seven
                # has no forced isolate; it is a weighted unicyclic triangle
                # with four rooted-tree edges.
                for active_vertices in range(4, 8):
                    isolates = 7 - active_vertices
                    cap = (
                        _weighted_triangle_seven_edge_cap()
                        if active_vertices == 7
                        else Fraction(649, 25)
                    )
                    upper = (
                        Fraction(exact_triangle_active6["maxima"]["edges11_component7"])
                        if active_vertices == 6
                        else _capped_upper(isolates, cap)
                    )
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "exact_determinant": active_vertices == 6,
                        "upper": upper,
                    })
                continue
            if edges == 12 and joined == 7:
                # Two outside K2 components remain.  The top two active
                # orders are unicyclic/bicyclic ten-edge components.
                for active_vertices in range(4, 11):
                    isolates = 11 - active_vertices
                    cap = (
                        _weighted_triangle_ten_edge_cap(active_vertices)
                        if active_vertices in (9, 10)
                        else Fraction(27)
                    )
                    upper = (
                        Fraction(exact_triangle["maxima"]["active10_unicyclic"])
                        if active_vertices == 10
                        else _capped_upper(isolates, cap)
                    )
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "exact_determinant": active_vertices == 10,
                        "upper": upper,
                    })
                continue
            if edges == 12 and joined == 8:
                # One outside K2 remains.  The maximal active orders are an
                # eleven-vertex unicyclic and ten-vertex bicyclic triangle
                # component; Perron edge-moving gives the exact caps above.
                for active_vertices in range(4, 12):
                    isolates = 13 - active_vertices
                    cap = (
                        _weighted_triangle_eleven_edge_cap(active_vertices)
                        if active_vertices in (10, 11)
                        else Fraction(2749, 100)
                    )
                    upper = _capped_upper(isolates, cap)
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "upper": upper,
                    })
                continue
            if edges == 12 and joined == 5:
                # Four outside edges are forced to be K2 components at the
                # maximal active order seven.  Enumerate every seven-vertex,
                # eight-edge triangle/no-K4 support and every gain phase.
                for active_vertices in range(4, 8):
                    isolates = 7 - active_vertices
                    cap = Fraction(649, 25)
                    upper = (
                        Fraction(exact_triangle["maxima"]["active7_bicyclic"])
                        if active_vertices == 7
                        else (
                            Fraction(exact_triangle_active6["maxima"]["edges12_component8"])
                            if active_vertices == 6
                            else _capped_upper(isolates, cap)
                        )
                    )
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "exact_determinant": active_vertices in (6, 7),
                        "upper": upper,
                    })
                continue
            if edges == 12 and joined == 6:
                # Three outside K2 components remain.  Active component
                # orders nine and eight are exactly the unicyclic/bicyclic
                # catalogues; smaller orders force at least two isolates.
                for active_vertices in range(4, 10):
                    isolates = 9 - active_vertices
                    cap = Fraction(2649, 100)
                    if active_vertices == 9:
                        upper = Fraction(exact_triangle["maxima"]["active9_unicyclic"])
                    elif active_vertices == 8:
                        upper = Fraction(exact_triangle["maxima"]["active8_bicyclic"])
                    else:
                        upper = _capped_upper(isolates, cap)
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "exact_determinant": active_vertices in (8, 9),
                        "upper": upper,
                    })
                continue
            if edges == 12 and joined == 9:
                # All twelve edges form one K4-free component.  Split by
                # active order.  Orders 11 and 12 are respectively bicyclic
                # and unicyclic and admit the exact Perron comparisons above;
                # smaller orders already win through their extra isolates.
                for active_vertices in range(4, 13):
                    isolates = 15 - active_vertices
                    cap = (
                        _weighted_triangle_sparse_cap(active_vertices)
                        if active_vertices in (11, 12)
                        else Fraction(2797, 100)
                    )
                    upper = _capped_upper(isolates, cap)
                    triangle_rows.append({
                        "edges": edges,
                        "joined_to_triangle": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": cap,
                        "upper": upper,
                    })
                continue
            isolates = max(0, 18 + joined - 2 * edges)
            cap = Fraction(15) + _sqrt_upper(Fraction(4 * component_energy, 3))
            upper = _capped_upper(isolates, cap)
            triangle_rows.append({
                "edges": edges,
                "joined_to_triangle": joined,
                "isolates": isolates,
                "lambda_cap": cap,
                "upper": upper,
            })
    no_k4 = max(triangle_free, *(row["upper"] for row in triangle_rows))
    assert no_k4 < q144
    k4_rows = []
    exact_k4_four_forest = q126_k4_four_forest_certificate()
    exact_isolated_k4 = q126_isolated_k4_certificate()
    exact_k4_three_forest = q126_k4_three_forest_certificate()
    exact_k4_leaf = q126_k4_leaf_certificate()
    exact_pure_k4_six_forest = q126_pure_k4_six_forest_certificate()
    for edges in supported_edges:
        for joined in range(edges - 5):
            if edges == 14:
                component_edges = 6 + joined
                outside_edges = 8 - joined
                vertices = (4,) if joined == 0 else range(5, 5 + joined)
                for active_vertices in vertices:
                    if active_vertices + 2 * outside_edges > 15:
                        continue
                    # A K5-free graph has at most the balanced four-partite
                    # Turan number of edges.  Skip impossible dense orders.
                    parts = [active_vertices // 4 + (index < active_vertices % 4) for index in range(4)]
                    turan_edges = (active_vertices**2 - sum(part**2 for part in parts)) // 2
                    if component_edges > turan_edges:
                        continue
                    isolates = max(0, 15 - active_vertices - 2 * outside_edges)
                    cap = Fraction(15) + _sqrt_upper(
                        Fraction(9 * (2 * component_edges - active_vertices + 1))
                    )
                    # K5 is absent in this branch.  Motzkin--Straus gives
                    # rho(E)^2<=189; 115/4 is a strict rational cap for the
                    # resulting Gram eigenvalue.
                    cap = min(cap, Fraction(115, 4))
                    tight_matching_cases = {(3, 5), (4, 7), (5, 9)}
                    one_isolate_matching_cases = {(4, 6), (5, 8)}
                    two_isolate_factor_case = (5, 7)
                    if (joined, active_vertices) in {(6, 9), (7, 11)}:
                        cap = min(cap, _pure_k4_two_isolate_cap(joined, active_vertices))
                    if (joined, active_vertices) == (7, 6):
                        cap = min(cap, _dense_pure_clique_cap("k4_order6"))
                    if (joined, active_vertices) == (8, 8):
                        cap = min(cap, _dense_pure_clique_cap("k4_order8"))
                    if (joined, active_vertices) in {(6, 8), (7, 10), (8, 12)}:
                        cap = min(cap, _pure_k4_three_isolate_cap(joined, active_vertices))
                    if (joined, active_vertices) in tight_matching_cases:
                        # Here active_vertices + 2*outside_edges = 15, so
                        # every edge outside the chosen connected component
                        # is forced to be an isolated K2.  Factor those K2
                        # determinants exactly and apply the trace envelope
                        # only to the active component.
                        upper = (
                            stationary_upper_dimension(active_vertices, 9 * component_edges)
                            * 216**outside_edges
                        )
                    elif (joined, active_vertices) in one_isolate_matching_cases:
                        # These endpoint-tight rows similarly force every
                        # outside edge to be a K2, now with one isolate.
                        upper = (
                            stationary_upper_dimension(active_vertices, 9 * component_edges)
                            * 216**outside_edges
                            * 15
                        )
                    elif (joined, active_vertices) == (6, 10):
                        upper = Fraction(exact_pure_k4_six_forest["maximum"])
                    elif (joined, active_vertices) == two_isolate_factor_case:
                        upper = (
                            stationary_upper_dimension(active_vertices, 9 * component_edges)
                            * 216**outside_edges
                            * 15**2
                        )
                    elif (joined, active_vertices) == (5, 6):
                        upper = (
                            stationary_upper_dimension(active_vertices, 9 * component_edges)
                            * 216**outside_edges
                            * 15**3
                        )
                    else:
                        upper = _capped_upper(isolates, cap)
                    assert upper < q144
                    k4_rows.append({
                        "edges": edges,
                        "joined_to_k4": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                        **(
                            {"exact_upper": [upper.numerator, upper.denominator]}
                            if (joined, active_vertices)
                            in tight_matching_cases
                            | one_isolate_matching_cases
                            | {(6, 10), two_isolate_factor_case, (5, 6)}
                            else {}
                        ),
                    })
                continue
            if edges == 12 and joined == 6:
                # All twelve support edges lie in the K4 component.  Orders
                # at most nine force at least six isolates, where the global
                # Motzkin cap suffices.  At order ten use the exact rooted-
                # forest comparison above.
                for active_vertices in range(5, 9):
                    cap = (
                        _seven_isolate_k4_cap("sqrt3_tricyclic")
                        if active_vertices == 8
                        else Fraction(115, 4)
                    )
                    isolates = 15 - active_vertices
                    upper = _capped_upper(isolates, cap)
                    assert upper < q144
                    k4_rows.append({
                        "edges": edges,
                        "joined_to_k4": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                    })
                cap = _weighted_k4_five_vertex_cap()
                upper = _capped_upper(6, cap)
                assert upper < q144
                k4_rows.append({
                    "edges": edges,
                    "joined_to_k4": joined,
                    "active_vertices": 9,
                    "isolates": 6,
                    "lambda_cap": [cap.numerator, cap.denominator],
                    "over_q144": float(upper / q144),
                })
                cap = _weighted_k4_six_leaf_cap()
                upper = _capped_upper(5, cap)
                assert upper < q144
                k4_rows.append({
                    "edges": edges,
                    "joined_to_k4": joined,
                    "active_vertices": 10,
                    "isolates": 5,
                    "lambda_cap": [cap.numerator, cap.denominator],
                    "over_q144": float(upper / q144),
                })
                continue
            if edges == 12 and joined == 4:
                # Two outside edges and maximal active order eight force a
                # four-edge rooted forest off K4, two K2 blocks, and three
                # isolates.  Its complete phase enumeration is exact.
                for active_vertices in range(5, 9):
                    isolates = 11 - active_vertices
                    cap = Fraction(2773, 100)
                    upper = (
                        Fraction(exact_k4_four_forest["maximum"])
                        if active_vertices == 8
                        else _capped_upper(isolates, cap)
                    )
                    k4_rows.append({
                        "edges": edges,
                        "joined_to_k4": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                        "exact_determinant": active_vertices == 8,
                        **(
                            {"exact_upper": [upper.numerator, upper.denominator]}
                            if active_vertices == 8
                            else {}
                        ),
                    })
                continue
            if edges == 12 and joined == 3:
                # At maximal active order seven, the three added edges form
                # a rooted forest; the three outside edges are K2 blocks and
                # two vertices are isolated.  Enumerate both heavy placements.
                for active_vertices in range(5, 8):
                    isolates = 9 - active_vertices
                    cap = Fraction(2719, 100)
                    upper = (
                        Fraction(exact_k4_three_forest["maximum"])
                        if active_vertices == 7
                        else _capped_upper(isolates, cap)
                    )
                    k4_rows.append({
                        "edges": edges,
                        "joined_to_k4": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                        "exact_determinant": active_vertices == 7,
                        **(
                            {"exact_upper": [upper.numerator, upper.denominator]}
                            if active_vertices == 7
                            else {}
                        ),
                    })
                continue
            if edges == 12 and joined == 2:
                # Four outside K2 blocks leave one isolate at maximal active
                # order six.  Split the heavy edge between the component and
                # outside blocks and apply the trace envelope componentwise.
                for active_vertices in (5, 6):
                    isolates = 7 - active_vertices
                    cap = Fraction(1331, 50)
                    if active_vertices == 6:
                        upper = max(
                            stationary_upper_dimension(6, 72) * 198 * 216**3 * 15,
                            stationary_upper_dimension(6, 90) * 216**4 * 15,
                        )
                    else:
                        upper = _capped_upper(isolates, cap)
                    k4_rows.append({
                        "edges": edges,
                        "joined_to_k4": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                        "component_trace_factorization": active_vertices == 6,
                        **(
                            {"exact_upper": [upper.numerator, upper.denominator]}
                            if active_vertices == 6
                            else {}
                        ),
                    })
                continue
            if edges == 12 and joined == 1:
                # Five outside K2 components remain.  Active order five is
                # exactly weighted K4-plus-leaf; active order four forces an
                # additional isolate and is already covered by the cap.
                for active_vertices in (5,):
                    isolates = 5 - active_vertices
                    cap = Fraction(2603, 100)
                    upper = (
                        Fraction(exact_k4_leaf["maximum"])
                        if active_vertices == 5
                        else _capped_upper(isolates, cap)
                    )
                    k4_rows.append({
                        "edges": edges,
                        "joined_to_k4": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                        "exact_determinant": active_vertices == 5,
                        **({"exact_upper": [upper.numerator, upper.denominator]} if active_vertices == 5 else {}),
                    })
                continue
            if edges == 11 and joined == 0:
                upper = Fraction(exact_isolated_k4["maxima"]["edges11_weight4"])
                k4_rows.append({
                    "edges": edges,
                    "joined_to_k4": joined,
                    "active_vertices": 4,
                    "isolates": 1,
                    "lambda_cap": [2603, 100],
                    "over_q144": float(upper / q144),
                    "exact_determinant": True,
                    "exact_upper": [upper.numerator, upper.denominator],
                })
                continue
            if edges == 12 and joined == 0:
                upper = Fraction(exact_isolated_k4["maxima"]["edges12_weight3"])
                k4_rows.append({
                    "edges": edges,
                    "joined_to_k4": joined,
                    "active_vertices": 4,
                    "isolates": 0,
                    "lambda_cap": [127, 5],
                    "over_q144": float(upper / q144),
                    "exact_determinant": True,
                    "exact_upper": [upper.numerator, upper.denominator],
                })
                continue
            if edges == 11 and joined == 4:
                # One outside K2 and maximal active order eight force the
                # same four-edge rooted-forest catalogue, now with a single
                # weight-4 edge and five isolates.
                for active_vertices in range(5, 9):
                    isolates = 13 - active_vertices
                    cap = Fraction(113, 4)
                    upper = (
                        Fraction(exact_k4_four_forest["weight4_maximum"])
                        if active_vertices == 8
                        else _capped_upper(isolates, cap)
                    )
                    k4_rows.append({
                        "edges": edges,
                        "joined_to_k4": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                        "exact_determinant": active_vertices == 8,
                        **(
                            {"exact_upper": [upper.numerator, upper.denominator]}
                            if active_vertices == 8
                            else {}
                        ),
                    })
                continue
            if edges == 11 and joined == 5:
                # The sole eleven-edge partition has one weight-4 edge.  At
                # active order nine the five edges beyond K4 form a rooted
                # forest adding five vertices, so Perron edge-moving reduces
                # to five leaves and moves the weight-2 scaled edge onto an
                # incident clique edge.  Smaller active orders retain the
                # global cap but force additional isolates.
                for active_vertices in range(5, 9):
                    cap = (
                        _seven_isolate_k4_cap("weight4_bicyclic")
                        if active_vertices == 8
                        else Fraction(115, 4)
                    )
                    isolates = 15 - active_vertices
                    upper = _capped_upper(isolates, cap)
                    assert upper < q144
                    k4_rows.append({
                        "edges": edges, "joined_to_k4": joined,
                        "active_vertices": active_vertices, "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                    })
                cap = _weighted_k4_five_leaf_cap()
                upper = _capped_upper(6, cap)
                assert upper < q144
                k4_rows.append({
                    "edges": edges, "joined_to_k4": joined,
                    "active_vertices": 9, "isolates": 6,
                    "lambda_cap": [cap.numerator, cap.denominator],
                    "over_q144": float(upper / q144),
                })
                continue
            if edges == 12 and joined == 5:
                # One edge lies outside the K4 component.  Split the active
                # component order.  At order nine the five added edges are a
                # rooted forest; at order eight there is one extra cycle and
                # the exact Perron comparison above applies.
                for active_vertices in range(5, 10):
                    isolates = 13 - active_vertices
                    if active_vertices == 9:
                        cap = _weighted_k4_five_leaf_cap()
                    elif active_vertices == 8:
                        cap = _weighted_k4_four_vertex_cap()
                    else:
                        cap = Fraction(113, 4)
                    upper = _capped_upper(isolates, cap)
                    assert upper < q144
                    k4_rows.append({
                        "edges": edges,
                        "joined_to_k4": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                    })
                continue
            if edges == 10 and joined == 4:
                # Both norm-27 edges and all ten support edges lie in the K4
                # component.  At maximal active order the four added edges
                # form a rooted forest; Perron shifting gives the displayed
                # two-adjacent-heavy-edge model.
                for active_vertices in range(5, 9):
                    isolates = 15 - active_vertices
                    cap = (
                        _seven_isolate_k4_cap("two_sqrt3_forest")
                        if active_vertices == 8
                        else Fraction(115, 4)
                    )
                    upper = _capped_upper(isolates, cap)
                    k4_rows.append({
                        "edges": edges,
                        "joined_to_k4": joined,
                        "active_vertices": active_vertices,
                        "isolates": isolates,
                        "lambda_cap": [cap.numerator, cap.denominator],
                        "over_q144": float(upper / q144),
                    })
                continue
            component_energy = 126 - 9 * (edges - 6 - joined)
            rho_square = Fraction(3 * component_energy, 2)
            cap = Fraction(15) + _sqrt_upper(rho_square)
            assert (cap - 15) ** 2 >= rho_square
            isolates = max(0, 23 + joined - 2 * edges)
            upper = _capped_upper(isolates, cap)
            assert upper < q144
            k4_rows.append({
                "edges": edges,
                "joined_to_k4": joined,
                "isolates": isolates,
                "lambda_cap": [cap.numerator, cap.denominator],
                "over_q144": float(upper / q144),
            })

    # If Q=126 contains K5, the prior component/isolate split is already
    # strong enough for Q=144.  Weighted supports have 10,11,12 edges.  The
    # pure 14-edge support has r=0,...,4 extra edges in the K5 component.
    weighted_k5_rows = {
        10: [fixed_isolate_capped_upper(10, 252, _weighted_k5_two_heavy_edges_cap())],
        11: [
            fixed_isolate_capped_upper(8, 252, Fraction(287, 10)),
            fixed_isolate_capped_upper(9, 252, _weighted_k5_one_leaf_cap()),
        ],
        12: [
            fixed_isolate_capped_upper(6, 252, Fraction(563, 20)),
            fixed_isolate_capped_upper(7, 252, Fraction(287, 10)),
            fixed_isolate_capped_upper(8, 252, _weighted_k5_two_leaf_cap()),
        ],
    }
    weighted_k5 = {edges: max(values) for edges, values in weighted_k5_rows.items()}

    # For the pure fourteen-edge support, split by r extra edges in the K5
    # component and by its active vertex count v.  Hong's connected-graph
    # bound rho(A)^2<=2m-v+1 is paired with the resulting isolate count.  The
    # only loose small-vertex instance is K6 minus one edge, whose Perron root
    # is (3+sqrt(41))/2 and is handled by the exact rational cap 29.11.
    assert Fraction(961, 150) ** 2 > 41
    dense_k6 = q126_k6_minus_edge_certificate()
    exact_pure_k5 = q126_pure_k5_certificate()
    pure_k5_rows = []
    for joined in range(5):
        vertices = (5,) if joined == 0 else range(6, 6 + joined)
        outside_edges = 4 - joined
        for active_vertices in vertices:
            isolates = max(0, 15 - active_vertices - 2 * outside_edges)
            cap = Fraction(15) + _sqrt_upper(
                Fraction(9 * (2 * (10 + joined) - active_vertices + 1))
            )
            if (joined, active_vertices) == (3, 6):
                cap = min(cap, _dense_pure_clique_cap("k5_order6"))
            if (joined, active_vertices) == (4, 8):
                cap = min(cap, _dense_pure_clique_cap("k5_order8"))
            if joined == 0:
                upper = Fraction(exact_pure_k5["maxima"]["joined0"])
            elif joined == 1:
                upper = Fraction(exact_pure_k5["maxima"]["joined1_leaf"])
            elif joined == 4 and active_vertices == 6:
                cap = Fraction(2911, 100)
                upper = Fraction(dense_k6["full_maximum"])
            else:
                if joined == 4 and active_vertices == 7:
                    cap = _pure_k5_two_vertex_cap()
                upper = fixed_isolate_capped_upper(isolates, 252, cap)
            pure_k5_rows.append((joined, active_vertices, isolates, cap, upper))
    pure_k5 = {
        joined: max(row[4] for row in pure_k5_rows if row[0] == joined)
        for joined in range(5)
    }
    assert all(value < q144 for value in weighted_k5.values())
    assert all(value < q144 for value in pure_k5.values())

    # The only later all-same slice is Q=135.  Fifteen support edges imply
    # clique number at most six, rho(E)^2<=225, lambda_max<=30.
    # Refine Q=135 by its spectral component.  Without K5, clique number at
    # most four gives lambda_max<29.25.  With K5 but no K6, component energy
    # and isolates handle weighted supports; the pure fifteen-edge support is
    # split by active order and Hong/Turan.  A K6 consumes all fifteen unit
    # edges and leaves nine isolates.
    # In the no-K5 branch weighted Motzkin--Straus with clique number four
    # gives rho(E)^2 <= 3Q/2 = 405/2.
    # Thus 29.24 is a strict rational Gram-eigenvalue cap.
    q135_no_k5_cap = Fraction(731, 25)
    assert (q135_no_k5_cap - 15) ** 2 > Fraction(405, 2)
    q135_rows = [_capped_upper(0, q135_no_k5_cap, 270)]
    support135 = sorted({len(weights) for weights in weight_partitions(15)})
    for edges in (value for value in support135 if value >= 10):
        for joined in range(edges - 9):
            outside_edges = edges - 10 - joined
            if edges == 15:
                component_edges = 10 + joined
                vertices = (5,) if joined == 0 else range(6, 6 + joined)
                for active_vertices in vertices:
                    parts = [active_vertices // 5 + (index < active_vertices % 5) for index in range(5)]
                    turan_edges = (active_vertices**2 - sum(part**2 for part in parts)) // 2
                    if component_edges > turan_edges:
                        continue
                    isolates = max(0, 15 - active_vertices - 2 * outside_edges)
                    cap = Fraction(15) + _sqrt_upper(
                        Fraction(9 * (2 * component_edges - active_vertices + 1))
                    )
                    q135_rows.append(_capped_upper(isolates, cap, 270))
                continue
            component_energy = 135 - 9 * outside_edges
            isolates = max(0, 30 + joined - 2 * edges)
            cap = Fraction(15) + _sqrt_upper(Fraction(8 * component_energy, 5))
            q135_rows.append(_capped_upper(isolates, cap, 270))
    q135_rows.append(_capped_upper(9, Fraction(30), 270))
    q135_all_same = max(q135_rows)
    assert q135_all_same < q144

    color_upper = max(
        Fraction(ENERGY9_REFINED_UPPER),
        Fraction(
            color["certified_excluding_energy9_numerator"],
            color["certified_excluding_energy9_denominator"],
        ),
    )
    energy45_upper = Fraction(
        color["size13_energy45_minimal_cross_numerator"],
        color["size13_energy45_minimal_cross_denominator"],
    )
    certified_upper = max(
        color_upper,
        energy45_upper,
        *size14.values(),
        no_k4,
        *weighted_k5.values(),
        *pure_k5.values(),
        q135_all_same,
    )
    # k4_rows store human-readable ratios; recompute their exact values for
    # the aggregate rather than rounding those ratios back to rationals.
    k4_exact = []
    for row in k4_rows:
        if "exact_upper" in row:
            k4_exact.append(Fraction(*row["exact_upper"]))
        else:
            cap = Fraction(*row["lambda_cap"])
            k4_exact.append(_capped_upper(row["isolates"], cap))
    certified_upper = max(certified_upper, *k4_exact)

    return {
        "q144_envelope_numerator": q144.numerator,
        "q144_envelope_denominator": q144.denominator,
        "size14_internal_maxima": internal,
        "size14_over_q144": {str(key): float(value / q144) for key, value in size14.items()},
        "q126_no_k4_over_q144": float(no_k4 / q144),
        "q126_triangle_free_over_q144": float(triangle_free / q144),
        "q126_triangle_rows": [
            {
                **{key: value for key, value in row.items() if key not in {"lambda_cap", "upper"}},
                "lambda_cap": [row["lambda_cap"].numerator, row["lambda_cap"].denominator],
                "over_q144": float(row["upper"] / q144),
            }
            for row in triangle_rows
        ],
        "q126_k4_max_over_q144": max(row["over_q144"] for row in k4_rows),
        "q126_k4_rows": k4_rows,
        "q126_weighted_k5_over_q144": {
            str(key): float(value / q144) for key, value in weighted_k5.items()
        },
        "q126_pure_k5_over_q144": {
            str(key): float(value / q144) for key, value in pure_k5.items()
        },
        "q135_all_same_over_q144": float(q135_all_same / q144),
        "certified_upper_numerator": certified_upper.numerator,
        "certified_upper_denominator": certified_upper.denominator,
        "theorem": "Every strict counterexample with Q<144 lies below the Q=144 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
