import itertools
import math
import random

from uoct.cobertura import cobertura
from uoct.modelo import exacto, voraz


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
