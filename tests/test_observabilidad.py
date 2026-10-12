import random

import networkx as nx
import numpy as np

from uoct import observabilidad as obs


def red_de(arcos, dentro, conserva):
    G = nx.MultiDiGraph()
    for x, y in arcos:
        G.add_edge(x, y, length=100.0, highway="primary")
    return obs.red(G, set(dentro), set(conserva))


def matriz_B(r):
    """Conservación en los nodos que conservan (el virtual, 0, queda fuera: su fila es redundante)."""
    B = np.zeros((max(r.n - 1, 1), len(r.u)))
    for t in range(len(r.u)):
        if r.a[t]:
            B[r.a[t] - 1, t] -= 1
        if r.b[t]:
            B[r.b[t] - 1, t] += 1
    return B


def algebra(r, M):
    """Referencia: rango de [B; I_M] − rango de B, y los tramos t con e_t en el espacio fila."""
    l, B = len(r.u), matriz_B(r)
    A = np.vstack([B] + [np.eye(l)[t] for t in sorted(M)])
    rk = np.linalg.matrix_rank(A)
    ident = {t for t in range(l) if np.linalg.matrix_rank(np.vstack([A, np.eye(l)[t]])) == rk}
    return rk - np.linalg.matrix_rank(B), ident


def test_en_un_ciclo_basta_medir_un_tramo():
    r = red_de([(1, 2), (2, 3), (3, 1)], dentro=[1, 2, 3], conserva=[1, 2, 3])
    assert obs.rango_completo(r) == 1
    assert obs.identificables(r, {0}) == {0, 1, 2}
    assert obs.rango(r, {0}) == 1


def test_la_camara_mide_los_tramos_de_su_cruce():
    r = red_de([(1, 2), (2, 3), (3, 1), (2, 2)], dentro=[1, 2, 3], conserva=[1, 2, 3])
    M = obs.medidos(r, [2])
    assert len(M) == 3                                   # (1,2), (2,3) y el lazo (2,2); no (3,1)
    assert all(2 in (r.u[t], r.v[t]) for t in M)


def test_sin_conservacion_solo_se_sabe_lo_medido():
    r = red_de([(1, 2), (2, 3), (3, 1)], dentro=[1, 2, 3], conserva=[])
    assert obs.identificables(r, {1}) == {1}


def test_coincide_con_el_algebra_en_redes_al_azar():
    azar = random.Random(13)
    for _ in range(300):
        nodos = list(range(1, 8))
        arcos = [(azar.choice(nodos), azar.choice(nodos)) for _ in range(azar.randint(3, 12))]
        dentro = azar.sample(nodos, azar.randint(1, 7))
        conserva = azar.sample(dentro, azar.randint(0, len(dentro)))
        r = red_de(arcos, dentro, conserva)
        if not r.u:
            continue
        M = set(azar.sample(range(len(r.u)), azar.randint(0, len(r.u))))
        rk, ident = algebra(r, M)
        assert obs.rango(r, M) == rk
        assert obs.identificables(r, M) == ident
        assert obs.rango_completo(r) == len(r.u) - np.linalg.matrix_rank(matriz_B(r))


def test_el_voraz_perezoso_elige_siempre_la_mayor_ganancia():
    azar = random.Random(7)
    nodos = list(range(1, 11))
    arcos = [(azar.choice(nodos), azar.choice(nodos)) for _ in range(30)]
    r = red_de(arcos, nodos, nodos)
    orden = obs.voraz(r, nodos, len(nodos))
    X, previo = [], 0
    for j, rk in orden:
        mejor = max(obs.rango(r, obs.medidos(r, X + [c])) for c in nodos if c not in X)
        assert rk == mejor and rk > previo
        X.append(j)
        previo = rk
    assert previo == obs.rango_completo(r)      # sin tope, llega a observarlo todo


def test_resumen_de_una_configuracion():
    r = red_de([(1, 2), (2, 3), (3, 1)], dentro=[1, 2, 3], conserva=[1, 2, 3])
    s = obs.resumen(r, [1])
    assert s["tramos medidos"] == 2 and s["tramos deducidos"] == 1
    assert s["tramos identificables %"] == 100.0 and s["rango %"] == 100.0
