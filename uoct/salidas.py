"""Tablas (csv y LaTeX) y figuras del informe."""
from pathlib import Path

import contextily as cx
import folium
import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import osmnx as ox  # noqa: E402
from shapely.geometry import Point  # noqa: E402


def tabla(df, carpeta, nombre):
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    df.to_csv(carpeta / f"{nombre}.csv", index=False)
    # Todo a texto con coma decimal: el `decimal` de pandas falla con columnas de enteros.
    tex = df.apply(lambda col: col.map(lambda v: f"{v:.1f}".replace(".", ",") if isinstance(v, float) else str(v)))
    (carpeta / f"{nombre}.tex").write_text(tex.to_latex(index=False, escape=True), encoding="utf-8")


def puntos(G, nodos):
    return gpd.GeoDataFrame(geometry=[Point(G.nodes[n]["x"], G.nodes[n]["y"]) for n in nodos], crs=G.graph["crs"])


def mapa(G, poligono, capas, ruta, titulo, flechas=()):
    """capas: [(nodos, estilo de GeoDataFrame.plot con label)]; flechas: [(origen, destino)]."""
    fig, ax = plt.subplots(figsize=(8, 8))
    ox.plot_graph(G, ax=ax, node_size=0, edge_color="#9a9a9a", edge_linewidth=0.4,
                  bgcolor="white", show=False, close=False)
    gpd.GeoSeries([poligono], crs=G.graph["crs"]).boundary.plot(ax=ax, color="black", linewidth=0.8)
    for nodos, estilo in capas:
        puntos(G, nodos).plot(ax=ax, **estilo)
    for o, d in flechas:
        ax.annotate("", xy=(G.nodes[d]["x"], G.nodes[d]["y"]), xytext=(G.nodes[o]["x"], G.nodes[o]["y"]),
                    arrowprops=dict(arrowstyle="->", color="tab:red", lw=1))
    x0, y0, x1, y1 = poligono.bounds
    ax.set_xlim(x0 - 150, x1 + 150)
    ax.set_ylim(y0 - 150, y1 + 150)
    cx.add_basemap(ax, crs=G.graph["crs"], source=cx.providers.Esri.WorldGrayCanvas, attribution_size=6)  # OSM y CARTO bloquean a los scripts
    ax.set_axis_off()
    ax.set_title(titulo)
    ax.legend(loc="lower left", fontsize=8)
    fig.savefig(ruta, dpi=200, bbox_inches="tight")
    plt.close(fig)


def figura_curva(exacta, intercambio, total, ruta):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot([b for b, _ in exacta], [100 * f / total for _, f in exacta], "o-", label="Óptimo con a lo más B traslados")
    ax.plot([k for k, *_ in intercambio], [100 * p[-1] / total for p in intercambio], "s--", label="Voraz de a una cámara")
    ax.set_xlabel("Traslados")
    ax.set_ylabel("Cobertura ponderada (%)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.savefig(ruta, dpi=200, bbox_inches="tight")
    plt.close(fig)


def figura_sensibilidad(df, ruta):
    fig, ax = plt.subplots(figsize=(6, 4))
    for (peso, g), color in zip(df[df["distancia"] == "no dirigida"].groupby("peso"), ["tab:blue", "tab:orange"]):
        ax.plot(g["R"], g["f(C) %"], "o--", color=color, label=f"Actual (peso {peso})")
        ax.plot(g["R"], g["exacto %"], "o-", color=color, label=f"Óptimo (peso {peso})")
    ax.set_xlabel("Radio de cobertura R (m, por la red)")
    ax.set_ylabel("Cobertura (%)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.savefig(ruta, dpi=200, bbox_inches="tight")
    plt.close(fig)


def mapa_html(G, capas, ruta, flechas=()):
    """capas: [(nombre, nodos, color)]; flechas: [(origen, destino)]."""
    m = folium.Map(location=[-33.45, -70.65], zoom_start=14, tiles="OpenStreetMap")
    for nombre, nodos, color in capas:
        grupo = folium.FeatureGroup(name=nombre)
        for p in puntos(G, nodos).to_crs(4326).geometry:
            folium.CircleMarker([p.y, p.x], radius=4, color=color, fill=True, fill_opacity=0.8).add_to(grupo)
        grupo.add_to(m)
    if flechas:
        grupo = folium.FeatureGroup(name="Traslados")
        for o, d in flechas:
            a, b = puntos(G, [o, d]).to_crs(4326).geometry
            folium.PolyLine([[a.y, a.x], [b.y, b.x]], color="red", weight=2).add_to(grupo)
        grupo.add_to(m)
    folium.LayerControl().add_to(m)
    m.save(str(ruta))


def figura_curva_obs(curva, rango_completo, p, ruta):
    """Rango alcanzado por el voraz de observabilidad contra el número de cámaras."""
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(range(1, len(curva) + 1), [100 * rk / rango_completo for _, rk in curva], "-", color="tab:blue")
    ax.axvline(p, color="gray", ls="--", lw=1)
    ax.text(p + 5, 3, f"{p} cámaras", color="gray", fontsize=8)
    ax.set_xlabel("Cámaras (voraz de observabilidad, sobre los cruces semaforizados)")
    ax.set_ylabel("Flujos independientes conocidos (% del rango completo)")
    ax.grid(alpha=0.3)
    fig.savefig(ruta, dpi=200, bbox_inches="tight")
    plt.close(fig)


def mapa_tramos(G, poligono, r, X, ruta, titulo):
    """Tramos medidos por las cámaras de X, deducidos por conservación y ciegos."""
    from matplotlib.lines import Line2D

    from .observabilidad import identificables, medidos
    M = medidos(r, X)
    deducidos = identificables(r, M) - M
    aristas = ox.graph_to_gdfs(G, nodes=False).loc[list(zip(r.u, r.v, r.clave))].reset_index(drop=True)
    estado = ["medido" if t in M else "deducido" if t in deducidos else "ciego" for t in range(len(r.u))]
    estilos = {"ciego": ("#c8c8c8", 0.5), "deducido": ("tab:blue", 1.8), "medido": ("tab:green", 1.8)}
    fig, ax = plt.subplots(figsize=(8, 8))
    for e, (color, ancho) in estilos.items():
        aristas[[s == e for s in estado]].plot(ax=ax, color=color, linewidth=ancho)
    gpd.GeoSeries([poligono], crs=G.graph["crs"]).boundary.plot(ax=ax, color="black", linewidth=0.8)
    puntos(G, set(X)).plot(ax=ax, color="black", markersize=10)
    x0, y0, x1, y1 = poligono.bounds
    ax.set_xlim(x0 - 150, x1 + 150)
    ax.set_ylim(y0 - 150, y1 + 150)
    cx.add_basemap(ax, crs=G.graph["crs"], source=cx.providers.Esri.WorldGrayCanvas, attribution_size=6)
    ax.set_axis_off()
    ax.set_title(titulo)
    n = {e: estado.count(e) for e in estilos}
    ax.legend(handles=[Line2D([], [], color="tab:green", lw=2, label=f"Medido por una cámara ({n['medido']})"),
                       Line2D([], [], color="tab:blue", lw=2, label=f"Deducido por conservación ({n['deducido']})"),
                       Line2D([], [], color="#c8c8c8", lw=2, label=f"Sin información ({n['ciego']})"),
                       Line2D([], [], color="black", marker="o", lw=0, label="Cámara")],
              loc="lower left", fontsize=8)
    fig.savefig(ruta, dpi=200, bbox_inches="tight")
    plt.close(fig)
