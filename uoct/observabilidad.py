"""Observabilidad de los flujos de tramo con las cámaras (Li et al., 2023; Zhuo et al., 2024).

Supuestos: una cámara cuenta sin error el flujo de todos los tramos de su cruce; el flujo se
conserva en los nodos de `conserva`, y todo lo demás (el exterior de la comuna y, en la variante,
los centroides) se contrae a un nodo virtual que también conserva. Así, un tramo no medido se
deduce si y solo si es un puente de la red sin los tramos medidos, y los flujos independientes
que se conocen son el rango cográfico |M| − comp(G − M) + comp(G)."""
import heapq
from collections import Counter
from dataclasses import dataclass

import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .cobertura import flujo_clase

VIRTUAL = 0


@dataclass
class Red:
    u: list            # extremos de cada tramo en la red de OSM
    v: list
    clave: list
    a: np.ndarray      # extremos en la red de conservación (0 = nodo virtual)
    b: np.ndarray
    n: int             # nodos de la red de conservación, con el virtual
    largo: np.ndarray  # metros
    flujo: np.ndarray  # tránsito de clase (veh/día): el mismo proxy que w_i
    incidentes: dict   # nodo de OSM -> tramos que lo tocan


def red(G, dentro, conserva):
    """Tramos: arcos dirigidos de G con algún extremo en `dentro` (la comuna).
    conserva ⊆ dentro; los demás nodos se contraen al nodo virtual."""
    # Un nodo fuera de `dentro` tiene solo parte de sus arcos en la red: conservar ahí inventa restricciones.
    assert set(conserva) <= set(dentro), "conserva debe estar dentro de la comuna"
    idx = {x: i + 1 for i, x in enumerate(sorted(conserva))}
    u, v, clave, largo, flujo, inc = [], [], [], [], [], {}
    for x, y, k, d in G.edges(keys=True, data=True):
        if x in dentro or y in dentro:
            t = len(u)
            u.append(x)
            v.append(y)
            clave.append(k)
            largo.append(d.get("length", 0.0))
            flujo.append(flujo_clase(d.get("highway")))
            for z in {x, y}:
                inc.setdefault(z, []).append(t)
    a = np.array([idx.get(x, VIRTUAL) for x in u], dtype=int)
    b = np.array([idx.get(y, VIRTUAL) for y in v], dtype=int)
    return Red(u, v, clave, a, b, len(conserva) + 1, np.array(largo, float), np.array(flujo, float), inc)


def medidos(r, X):
    """Tramos que cuentan las cámaras de los nodos X: todos los de su cruce."""
    return {t for j in X for t in r.incidentes.get(j, [])}


def _componentes(r, libre):
    m = coo_matrix((np.ones(int(libre.sum())), (r.a[libre], r.b[libre])), shape=(r.n, r.n))
    return connected_components(m, directed=False)[0]


def _todos(r):
    return np.ones(len(r.u), bool)


def rango(r, M):
    """Flujos independientes que se conocen midiendo M (más allá de la conservación)."""
    libre = _todos(r)
    libre[list(M)] = False
    return len(M) - _componentes(r, libre) + _componentes(r, _todos(r))


def rango_completo(r):
    """Lo que hay que conocer para deducirlo todo: la dimensión del espacio de ciclos."""
    return len(r.u) - r.n + _componentes(r, _todos(r))


def _par(r, t):
    return (min(r.a[t], r.b[t]), max(r.a[t], r.b[t]))


def identificables(r, M):
    """M y los tramos no medidos que no están en ningún ciclo de la red sin M (puentes)."""
    libres = [t for t in range(len(r.u)) if t not in M and r.a[t] != r.b[t]]
    pares = Counter(_par(r, t) for t in libres)
    H = nx.Graph()
    H.add_nodes_from(range(r.n))
    H.add_edges_from(pares)
    puentes = {(min(x, y), max(x, y)) for x, y in nx.bridges(H)}
    return set(M) | {t for t in libres if pares[_par(r, t)] == 1 and _par(r, t) in puentes}


def voraz(r, candidatos, B):
    """Voraz perezoso sobre el rango, que es submodular: garantía 1 − 1/e (Li et al., 2023).
    Devuelve [(nodo, rango acumulado)]; para cuando ninguna cámara más aporta."""
    c0, libre, M, actual, orden = _componentes(r, _todos(r)), _todos(r), set(), 0, []

    def ganancia(j):
        nuevos = [t for t in r.incidentes.get(j, []) if t not in M]
        if not nuevos:
            return 0, nuevos
        prueba = libre.copy()
        prueba[nuevos] = False
        return len(M) + len(nuevos) - _componentes(r, prueba) + c0 - actual, nuevos

    pila = [(-ganancia(j)[0], j) for j in candidatos]
    heapq.heapify(pila)
    while pila and len(orden) < B:
        _, j = heapq.heappop(pila)
        g, nuevos = ganancia(j)  # la ganancia vieja es cota superior: solo puede haber bajado
        if pila and g < -pila[0][0]:
            heapq.heappush(pila, (-g, j))
            continue
        if g <= 0:
            break
        M.update(nuevos)
        libre[nuevos] = False
        actual += g
        orden.append((j, actual))
    return orden


def resumen(r, X):
    """Cuánto se deduce con las cámaras de X: tramos, km, tránsito (proxy) y rango."""
    M = medidos(r, X)
    ident = list(identificables(r, M))
    return {"cámaras": len(set(X)), "tramos medidos": len(M), "tramos deducidos": len(ident) - len(M),
            "tramos identificables %": 100 * len(ident) / len(r.u),
            "km identificables %": 100 * r.largo[ident].sum() / r.largo.sum(),
            "tránsito identificable %": 100 * r.flujo[ident].sum() / r.flujo.sum(),
            "rango %": 100 * rango(r, M) / rango_completo(r)}
