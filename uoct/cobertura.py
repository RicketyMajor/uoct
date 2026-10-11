"""Distancias por la red entre cruces, vecindades N_i, pesos w_i y la cobertura f(X)."""
import networkx as nx

# Tránsito medio diario (veh/día) por clase de vía de OSM, medido en 10.698 puntos de conteo
# europeos: Shen et al. (2024), Atmospheric Environment 335, Tabla 3. El resto vale «residential».
FLUJO_CLASE = {"motorway": 46632, "trunk": 46632, "primary": 13022, "secondary": 5852, "tertiary": 5852}
FLUJO_RESTO = 1392


def distancias(G, F, rmax, dirigida=False):
    """d[j][i] = metros de i a j por la red, solo si ≤ rmax.
    No dirigida: la red sin sentidos de tránsito. Dirigida: el camino de i hacia j,
    es decir, i aguas arriba de j (la congestión en j crece hacia los cruces que la alimentan)."""
    H = G.reverse(copy=False) if dirigida else G.to_undirected(as_view=True)
    Fs = set(F)
    d = {}
    for j in F:
        lon = nx.single_source_dijkstra_path_length(H, j, cutoff=rmax, weight="length")
        d[j] = {i: m for i, m in lon.items() if i in Fs}
    return d


def vecinos(d, R):
    """N[i] = ubicaciones j desde donde una cámara cubre el cruce i (d ≤ R)."""
    N = {}
    for j, fila in d.items():
        for i, m in fila.items():
            if m <= R:
                N.setdefault(i, set()).add(j)
    return N


def cubre_desde(N):
    """Lo inverso de N: qué cruces cubre una cámara en j."""
    cubre = {}
    for i, js in N.items():
        for j in js:
            cubre.setdefault(j, set()).add(i)
    return cubre


def flujo_clase(highway):
    h = highway[0] if isinstance(highway, list) else highway
    return FLUJO_CLASE.get(str(h).removesuffix("_link"), FLUJO_RESTO)


def pesos_jerarquia(G, F):
    """w_i provisional (hasta tener el flujo de la UOCT): la mitad de la suma del tránsito
    de clase de cada rama que llega al cruce. Cada vehículo entra por una rama y sale por
    otra, así que la mitad aproxima el volumen que entra al cruce (veh/día)."""
    Gu = G.to_undirected(as_view=True)
    w = {}
    for i in F:
        ramas = {}
        for _, v, dat in Gu.edges(i, data=True):
            if v != i:
                ramas[v] = max(ramas.get(v, 0), flujo_clase(dat.get("highway")))
        w[i] = sum(ramas.values()) / 2
    return w


def cobertura(X, N, w):
    """f(X): suma de w_i de los cruces con al menos una cámara de X a distancia ≤ R."""
    X = set(X)
    return sum(w[i] for i, js in N.items() if js & X)
