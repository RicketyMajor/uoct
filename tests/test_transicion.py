import itertools

import networkx as nx
from uoct.cobertura import cobertura
from uoct.transicion import hungaro, por_mover, voraz_intercambio


def linea():
    G = nx.MultiDiGraph(crs="EPSG:32719")
    for i in range(5):
        G.add_node(i, x=100.0 * i, y=0.0)
    for u in range(4):
        G.add_edge(u, u + 1, length=100.0)
        G.add_edge(u + 1, u, length=100.0)
    return G


def test_por_mover_incluye_las_repetidas_y_las_que_salen():
    assert sorted(por_mover([0, 0, 1, 4], X={0, 2, 3})) == [0, 1, 4]


def test_hungaro_iguala_a_la_fuerza_bruta():
    G, orig, dest = linea(), [0, 4, 1], [2, 3]
    plan = hungaro(G, orig, dest)
    costo = sum(m for *_, m in plan)
    bruto = min(sum(abs(o - d) * 100 for o, d in zip(sel, dest))
                for sel in itertools.permutations(orig, len(dest)))
    assert costo == bruto and len(plan) == 2


def test_voraz_intercambio_nunca_baja_y_mueve_de_a_una():
    F = list(range(6))
    N = {0: {0}, 1: {1}, 2: {2}, 3: {3, 4}, 4: {4}, 5: {5}}
    w = {0: 1, 1: 1, 2: 5, 3: 1, 4: 5, 5: 1}
    pasos = voraz_intercambio(N, w, F, [0, 0, 1])     # dos cámaras en el mismo cruce
    valores = [f for *_, f in pasos]
    assert valores[0] == cobertura([0, 1], N, w)
    assert all(a < b for a, b in zip(valores, valores[1:]))
    assert pasos[1][1] == 0                             # primero sale la repetida (aporta 0)
