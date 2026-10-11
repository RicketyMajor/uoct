"""Red vial de OSM y los conjuntos F (cruces semaforizados) y C (cámaras sobre F)."""
from pathlib import Path

import geopandas as gpd
import numpy as np
import osmnx as ox
from scipy.spatial import cKDTree

from .datos import CRS_UTM

LUGAR = "Santiago, Región Metropolitana de Santiago, Chile"


def descargar(carpeta, lugar=LUGAR, buffer_m=500):
    """Baja de OSM la comuna, su red para autos y sus semáforos, y los deja en `carpeta`.
    La red lleva un margen de buffer_m alrededor de la comuna para que las distancias
    de los cruces del borde no se corten en el límite comunal."""
    carpeta = Path(carpeta)
    comuna = ox.geocode_to_gdf(lugar).to_crs(CRS_UTM)
    comuna[["geometry"]].to_file(carpeta / "comuna.gpkg")
    zona = comuna.buffer(buffer_m).to_crs(4326).iloc[0]
    ox.save_graphml(ox.graph_from_polygon(zona, network_type="drive"), carpeta / "red.graphml")
    sem = ox.features_from_polygon(comuna.to_crs(4326).geometry.iloc[0], tags={"highway": "traffic_signals"})
    sem[sem.geom_type == "Point"][["geometry"]].reset_index().to_file(carpeta / "semaforos.gpkg")


def cargar(carpeta):
    carpeta = Path(carpeta)
    G = ox.project_graph(ox.load_graphml(carpeta / "red.graphml"), to_crs=CRS_UTM)
    comuna = gpd.read_file(carpeta / "comuna.gpkg").to_crs(CRS_UTM)
    sem = gpd.read_file(carpeta / "semaforos.gpkg").to_crs(CRS_UTM)
    return G, comuna.geometry.iloc[0], sem


def cruces_F(G, semaforos, poligono, max_m=50):
    """F: intersecciones a ≤ max_m de algún semáforo de la comuna. OSM marca un semáforo
    por cada acceso al cruce; asignar cada uno a su intersección más cercana los agrupa."""
    sem = semaforos[semaforos.within(poligono)]
    n, d = ox.distance.nearest_nodes(G, sem.geometry.x.to_numpy(), sem.geometry.y.to_numpy(), return_dist=True)
    n, d = np.asarray(n), np.asarray(d)
    return sorted(set(n[d <= max_m].tolist()))


def asignar_camaras(G, F, camaras, poligono, fuera_m=300):
    """Cada cámara va al cruce de F más cercano. Se descartan las que caen a más de
    fuera_m del límite comunal (comuna mal etiquetada en el inventario)."""
    dentro = camaras[camaras.distance(poligono) <= fuera_m].copy()
    xy = np.array([(G.nodes[j]["x"], G.nodes[j]["y"]) for j in F])
    d, k = cKDTree(xy).query(np.c_[dentro.geometry.x, dentro.geometry.y])
    dentro["nodo"] = np.asarray(F)[k]
    dentro["metros"] = d
    return dentro.reset_index(drop=True)
