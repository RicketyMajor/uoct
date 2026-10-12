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
