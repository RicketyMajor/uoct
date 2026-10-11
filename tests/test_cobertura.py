import networkx as nx
from uoct.cobertura import cobertura, distancias, pesos_jerarquia, vecinos


def linea(doble=True):
    G = nx.MultiDiGraph(crs="EPSG:32719")
    for i in range(4):
        G.add_node(i, x=100.0 * i, y=0.0)
    for u, v in [(0, 1), (1, 2), (2, 3)]:
        G.add_edge(u, v, length=100.0, highway="primary", name="Alameda")
        if doble:
            G.add_edge(v, u, length=100.0, highway="primary", name="Alameda")
    return G


def test_distancias_con_corte_y_vecinos():
    d = distancias(linea(), [0, 1, 2, 3], rmax=150)
    assert d[0] == {0: 0.0, 1: 100.0}
    N = vecinos(d, 100)
    assert N[1] == {0, 1, 2}
    w = {i: 1 for i in range(4)}
    assert cobertura({1}, N, w) == 3
    assert cobertura({1}, vecinos(d, 0), w) == 1          # R = 0: solo el cruce de la cámara


def test_distancia_dirigida_mira_aguas_arriba():
    d = distancias(linea(doble=False), [0, 1, 2, 3], rmax=250, dirigida=True)
    assert d[2] == {2: 0.0, 1: 100.0, 0: 200.0}           # se llega a 2 desde 1 y 0, no desde 3


def test_peso_es_la_mitad_del_flujo_de_las_ramas():
    G = linea()
    G.add_node(9, x=100.0, y=100.0)
    G.add_edge(1, 9, length=100.0, highway="residential", name="Calle")
    G.add_edge(9, 1, length=100.0, highway="residential_link", name="Calle")
    w = pesos_jerarquia(G, [0, 1])
    assert w[1] == (13022 + 13022 + 1392) / 2   # dos ramas de primaria y una residencial
    assert w[0] == 13022 / 2                    # extremo: una sola rama
