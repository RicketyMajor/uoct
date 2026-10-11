"""La transición: qué cámara va a dónde (húngaro) y el voraz de a una cámara."""
from collections import Counter

import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment

from .cobertura import cobertura, cubre_desde


def por_mover(nodos_camaras, X):
    """Cámaras que pueden moverse (por su nodo actual): las que quedan fuera de X y,
    si un cruce de X tiene varias, todas menos una."""
    vistos, salen = set(), []
    for n in nodos_camaras:
        if n in X and n not in vistos:
            vistos.add(n)
        else:
            salen.append(n)
    return salen


def hungaro(G, origenes, destinos):
    """Asignación de costo mínimo (Kuhn, 1955) entre cámaras por mover y ubicaciones nuevas.
    Con un costo fijo por traslado, el número de traslados no depende de la asignación:
    la distancia por la red solo decide qué cámara va a cada lugar."""
    H = G.to_undirected(as_view=True)
    filas = []
    for o in origenes:
        lon = nx.single_source_dijkstra_path_length(H, o, weight="length")
        filas.append([lon.get(d, 1e9) for d in destinos])
    c = np.array(filas)
    r, k = linear_sum_assignment(c)                    # rectangular: elige las cámaras más cercanas
    return [(origenes[a], destinos[b], float(c[a, b])) for a, b in zip(r, k)]


def voraz_intercambio(N, w, F, nodos_camaras, delta=0.0):
    """Coaching ▶ 18:34–20:19: sacar la cámara que menos aporta y ponerla donde más aporta,
    si la ganancia neta es ≥ delta y positiva. Cada paso es un traslado."""
    cubre = cubre_desde(N)
    X = list(nodos_camaras)                            # una entrada por cámara
    cnt = Counter(i for j in X for i in cubre.get(j, ()))
    f = cobertura(X, N, w)
    pasos = [(0, None, None, f)]
    while len(pasos) <= len(X):
        perdida = {a: sum(w[i] for i in cubre.get(X[a], ()) if cnt[i] == 1) for a in range(len(X))}
        a = min(perdida, key=perdida.get)
        sale = X[a]
        for i in cubre.get(sale, ()):
            cnt[i] -= 1
        ocupados = set(X[:a] + X[a + 1:])
        ganancia = {l: sum(w[i] for i in cubre.get(l, ()) if cnt[i] == 0) for l in F if l not in ocupados}
        entra = max(ganancia, key=ganancia.get)
        neto = ganancia[entra] - perdida[a]
        if neto <= 0 or neto < delta:
            for i in cubre.get(sale, ()):
                cnt[i] += 1
            break
        for i in cubre.get(entra, ()):
            cnt[i] += 1
        X[a] = entra
        f += neto
        pasos.append((len(pasos), sale, entra, f))
    return pasos
