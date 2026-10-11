import itertools
import math
import random

from uoct.cobertura import cobertura
from uoct.modelo import curva_B, exacto, lexicografico, voraz


def instancia(semilla, n=9):
    rnd = random.Random(semilla)
    F = list(range(n))
    N = {i: {i} | {j for j in F if rnd.random() < 0.25} for i in F}
    w = {i: rnd.randint(1, 8) for i in F}
    return N, w, F


def optimo(N, w, F, p):
    return max(cobertura(X, N, w) for X in itertools.combinations(F, p))


def test_exacto_iguala_a_la_fuerza_bruta():
    for s in range(5):
        N, w, F = instancia(s)
        X, f = exacto(N, w, F, 3)
        assert len(X) <= 3 and f == optimo(N, w, F, 3)


def test_voraz_cumple_la_garantia():
    for s in range(5):
        N, w, F = instancia(s)
        assert cobertura(voraz(N, w, F, 3), N, w) >= (1 - 1 / math.e) * optimo(N, w, F, 3)


def test_lexicografico_maxima_cobertura_y_luego_lo_que_se_queda():
    for s in range(5):
        N, w, F = instancia(s)
        C = {0, 1, 2}
        X, f, quedan = lexicografico(N, w, F, 3, C)
        fstar = optimo(N, w, F, 3)
        assert f == fstar
        mejor = max(len(set(Y) & C) for Y in itertools.combinations(F, 3) if cobertura(Y, N, w) == fstar)
        assert quedan == mejor


def test_curva_B_crece_y_llega_al_optimo():
    N, w, F = instancia(3)
    C = {0, 1, 2}
    curva = curva_B(N, w, F, 3, C, Bmax=3)
    valores = [f for _, f in curva]
    assert valores[0] == cobertura(C, N, w)                          # B = 0: quedarse como está
    assert all(a <= b for a, b in zip(valores, valores[1:]))
    assert valores[-1] == optimo(N, w, F, 3)                         # con 3 traslados se llega a f*
