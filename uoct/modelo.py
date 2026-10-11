"""MCLP con p cámaras (Church y ReVelle, 1974): voraz, exacto, lexicográfico y curva f*(B)."""
import pyomo.environ as pyo

from .cobertura import cobertura, cubre_desde

SOLVERS = {"highs": "appsi_highs", "gurobi": "gurobi"}
BRECHA = {"highs": "mip_rel_gap", "gurobi": "MIPGap"}


def voraz(N, w, F, p):
    """Agrega de a una la ubicación que más cobertura nueva aporta.
    Garantía: al menos 1 − 1/e ≈ 63 % del óptimo (Nemhauser, Wolsey y Fisher, 1978)."""
    cubre = cubre_desde(N)
    X, cubiertos = [], set()
    for _ in range(min(p, len(F))):
        j = max((j for j in F if j not in X), key=lambda j: sum(w[i] for i in cubre.get(j, set()) - cubiertos))
        X.append(j)
        cubiertos |= cubre.get(j, set())
    return X


def _modelo(N, w, F, p):
    m = pyo.ConcreteModel()
    m.x = pyo.Var(F, domain=pyo.Binary)                  # x_j = 1: cámara en j
    m.y = pyo.Var(list(N), bounds=(0, 1))                # y_i = 1: cruce i cubierto
    m.cob = pyo.Constraint(list(N), rule=lambda m, i: m.y[i] <= sum(m.x[j] for j in N[i]))
    m.p = pyo.Constraint(expr=sum(m.x[j] for j in F) <= p)
    m.f = pyo.Expression(expr=sum(w[i] * m.y[i] for i in N))
    return m


def _resolver(m, solver):
    s = pyo.SolverFactory(SOLVERS[solver])
    s.options[BRECHA[solver]] = 0                        # óptimo exacto, no a un 0,01 %
    s.solve(m)
    return {j for j in m.x if pyo.value(m.x[j]) > 0.5}


def exacto(N, w, F, p, C=(), B=None, solver="highs"):
    """Óptimo del MCLP. Con B, a lo más B ubicaciones fuera de C (los traslados).
    Σ x_j ≤ p y no = p: la cobertura no baja al agregar cámaras, y con B = 0 el modelo
    sigue factible aunque C tenga dos cámaras en un mismo cruce."""
    m = _modelo(N, w, F, p)
    if B is not None:
        fuera = [j for j in F if j not in set(C)]
        m.B = pyo.Constraint(expr=sum(m.x[j] for j in fuera) <= B)
    m.obj = pyo.Objective(expr=m.f, sense=pyo.maximize)
    X = _resolver(m, solver)
    return X, cobertura(X, N, w)


def lexicografico(N, w, F, p, C, solver="highs"):
    """Primero la máxima cobertura f*; entre las soluciones con f*, la que conserva más
    ubicaciones de C. Es el orden lexicográfico del coaching (▶ 15:26–17:31), resuelto en
    dos fases en vez de ponderar por una cota: el mismo orden, sin errores numéricos con
    pesos reales."""
    _, fstar = exacto(N, w, F, p, solver=solver)
    m = _modelo(N, w, F, p)
    m.piso = pyo.Constraint(expr=m.f >= fstar - 1e-6 * max(1.0, fstar))
    m.obj = pyo.Objective(expr=sum(m.x[j] for j in F if j in set(C)), sense=pyo.maximize)
    X = _resolver(m, solver)
    return X, cobertura(X, N, w), len(X & set(C))


def curva_B(N, w, F, p, C, Bmax, solver="highs"):
    """f*(B) para B = 0..Bmax: la mejor cobertura con a lo más B ubicaciones nuevas."""
    return [(B, exacto(N, w, F, p, C, B, solver)[1]) for B in range(Bmax + 1)]
