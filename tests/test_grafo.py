import geopandas as gpd
import networkx as nx
from shapely.geometry import Point, box
from uoct.grafo import asignar_camaras, cruces_F


def grafo_juguete():
    """Cuatro cruces en línea, cada 100 m, calle de doble sentido."""
    G = nx.MultiDiGraph(crs="EPSG:32719")
    for i, x in enumerate([0, 100, 200, 300]):
        G.add_node(i, x=float(x), y=0.0)
    for u, v in [(0, 1), (1, 2), (2, 3)]:
        G.add_edge(u, v, length=100.0, highway="primary", name="Alameda")
        G.add_edge(v, u, length=100.0, highway="primary", name="Alameda")
    return G


def puntos(xy):
    return gpd.GeoDataFrame(geometry=[Point(x, y) for x, y in xy], crs=32719)


def test_cruces_F_agrupa_semaforos_por_interseccion():
    sem = puntos([(5, 3), (95, 0), (100, 8), (260, 0), (150, 40), (9000, 0)])
    # (150, 40) queda a 64 m de todo cruce; (9000, 0) está fuera de la comuna
    F = cruces_F(grafo_juguete(), sem, box(-50, -50, 350, 50))
    assert F == [0, 1, 3]


def test_asignar_camaras_al_cruce_mas_cercano_y_descarta_las_de_fuera():
    cams = puntos([(10, 0), (190, 0), (290, 5), (5000, 0)])
    a = asignar_camaras(grafo_juguete(), [0, 1, 3], cams, box(-50, -50, 350, 50))
    assert list(a["nodo"]) == [0, 1, 3]                 # la de (5000, 0) se descarta
    assert round(a["metros"].iloc[1]) == 90
