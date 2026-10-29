"""Todo el prototipo con un comando:  python run.py  (o  python run.py --descargar  para volver a bajar OSM).
Escribe out/tablas, out/figuras y out/mapa.html, y contrasta con el Estado del proyecto."""
import argparse
import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd
from scipy.spatial import cKDTree
from shapely.geometry import Point

from uoct import cobertura as cob
from uoct import datos, grafo, modelo, salidas, transicion
from uoct import observabilidad as obs

RAIZ = Path(__file__).parent
DATA, OUT = RAIZ / "data", RAIZ / "out"
INVENTARIO = DATA / "Inventario_CCTV_RM_Actualizado_31-08-2026.xlsx"
RADIOS = [0, 50, 100, 150, 200, 250, 300, 400]   # m, por la red (50–100: el alcance visual de Boldt et al., 2026)
R_BASE = 300                            # m: el que usan Dell'Olmo et al. (2014) y Wheeler (2026); ver el informe
DELTAS = [0, 0.25, 0.5, 1.0]            # ganancia neta mínima para trasladar, en % de la cobertura total
B_PLAN = 10                             # traslados del plan que se dibuja aparte (el codo de la curva)
SOLVER = "highs"
ESTADO = {"F (cruces semaforizados)": 527, "C (cámaras)": 56, "C en línea": 43,
          "Cruces con cámara": 54, "Cámaras a < 30 m de un cruce (%)": 93}


def fila(R, peso, distancia, N, w, F, C, C_ef):
    total = sum(w.values())
    Xv = modelo.voraz(N, w, F, len(C))
    _, fe = modelo.exacto(N, w, F, len(C), solver=SOLVER)
    Xl, _, _ = modelo.lexicografico(N, w, F, len(C), set(C), solver=SOLVER)
    fv = cob.cobertura(Xv, N, w)
    return {"R": R, "peso": peso, "distancia": distancia,
            "f(C) %": 100 * cob.cobertura(C, N, w) / total, "f(C en línea) %": 100 * cob.cobertura(C_ef, N, w) / total,
            "voraz %": 100 * fv / total, "exacto %": 100 * fe / total,
            "brecha voraz %": 100 * (fe - fv) / fe if fe else 0.0,
            "traslados lexicográfica": len(Xl - set(C))}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--descargar", action="store_true", help="vuelve a bajar la red y los semáforos de OSM")
    if ap.parse_args().descargar:
        grafo.descargar(DATA)

    G, pol, sem = grafo.cargar(DATA)
    F = grafo.cruces_F(G, sem, pol)
    cams = grafo.asignar_camaras(G, F, datos.camaras_comuna(datos.leer_inventario(INVENTARIO)), pol)
    C, C_ef = list(cams["nodo"]), list(cams.loc[cams["en_linea"], "nodo"])
    tablas, figuras = OUT / "tablas", OUT / "figuras"
    figuras.mkdir(parents=True, exist_ok=True)

    medido = {"F (cruces semaforizados)": len(F), "C (cámaras)": len(cams), "C en línea": len(C_ef),
              "Cruces con cámara": len(set(C)),
              "Cámaras a < 30 m de un cruce (%)": round(100 * (cams["metros"] < 30).mean(), 1)}
    control = pd.DataFrame({"indicador": list(ESTADO), "Estado del proyecto": list(ESTADO.values()),
                            "este prototipo": pd.Series([medido[k] for k in ESTADO], dtype=object)})
    salidas.tabla(control, tablas, "control")
    tramos = [d["length"] for u, v, d in G.to_undirected().edges(data=True)
              if pol.contains(Point(G.nodes[u]["x"], G.nodes[u]["y"])) and pol.contains(Point(G.nodes[v]["x"], G.nodes[v]["y"]))]
    inv = datos.leer_inventario(INVENTARIO).dropna(subset=["lat", "lon"])
    dentro_pol = gpd.GeoSeries(gpd.points_from_xy(inv["lon"], inv["lat"]), crs=4326).to_crs(G.graph["crs"]).within(pol)
    ajenas = int((dentro_pol.to_numpy() & (inv["comuna"].fillna("") != "Santiago").to_numpy()).sum())
    pares_F = len(cKDTree([(G.nodes[j]["x"], G.nodes[j]["y"]) for j in F]).query_pairs(30))
    red = pd.DataFrame({"indicador": ["Nodos (con 500 m de margen)", "Arcos dirigidos (con margen)",
                                      "Calles dentro de la comuna (sin sentido, los dos extremos dentro)",
                                      "Largo mediano de esas calles (m)",
                                      "Pares de nodos de F a menos de 30 m entre sí",
                                      "Cámaras a más de 30 m de su cruce", "Distancia máxima cámara-cruce (m)",
                                      "Cámaras dentro de la comuna sin rótulo «Santiago» (no se usan)"],
                        "valor": pd.Series([G.number_of_nodes(), G.number_of_edges(), len(tramos),
                                            round(float(pd.Series(tramos).median()), 1), pares_F,
                                            int((cams["metros"] > 30).sum()), round(float(cams["metros"].max()), 1),
                                            ajenas], dtype=object)})
    salidas.tabla(red, tablas, "red")

    pesos = {"jerarquía": cob.pesos_jerarquia(G, F), "uniforme": {i: 1 for i in F}}
    d = cob.distancias(G, F, max(RADIOS))
    filas = [fila(R, nombre, "no dirigida", cob.vecinos(d, R), w, F, C, C_ef)
             for R in RADIOS for nombre, w in pesos.items()]
    d_dir = cob.distancias(G, F, R_BASE, dirigida=True)
    filas.append(fila(R_BASE, "jerarquía", "dirigida", cob.vecinos(d_dir, R_BASE), pesos["jerarquía"], F, C, C_ef))
    sens = pd.DataFrame(filas)
    salidas.tabla(sens, tablas, "sensibilidad")

    # Caso base: R_BASE, peso por jerarquía, distancia no dirigida
    N, w = cob.vecinos(d, R_BASE), pesos["jerarquía"]
    total = sum(w.values())
    Xl, fl, quedan = modelo.lexicografico(N, w, F, len(C), set(C), solver=SOLVER)
    curva = modelo.curva_B(N, w, F, len(C), set(C), len(Xl - set(C)), solver=SOLVER)
    salidas.tabla(pd.DataFrame([(B, 100 * f / total) for B, f in curva], columns=["B", "f*(B) %"]), tablas, "curva")
    plan = transicion.hungaro(G, transicion.por_mover(C, Xl), sorted(Xl - set(C)))
    salidas.tabla(pd.DataFrame(plan, columns=["desde (nodo)", "hasta (nodo)", "metros"]), tablas, "traslados")
    inter = {pct: transicion.voraz_intercambio(N, w, F, C, pct / 100 * total) for pct in DELTAS}
    salidas.tabla(pd.DataFrame([{"delta %": f"{k:g}".replace(".", ","), "traslados": len(v) - 1, "f %": 100 * v[-1][-1] / total}
                                for k, v in inter.items()]), tablas, "intercambio")

    fuera = set(cams.loc[~cams["en_linea"], "nodo"])
    salidas.mapa(G, pol, [(F, dict(color="#bbbbbb", markersize=6, label=f"Cruces semaforizados (F = {len(F)})")),
                          (set(C_ef), dict(color="tab:green", markersize=22, label="Cámara en línea")),
                          (fuera, dict(color="tab:red", markersize=22, marker="x", label="Cámara fuera de línea"))],
                 figuras / "estado.png", "Situación actual")
    salidas.mapa(G, pol, [(set(C) - Xl, dict(color="tab:red", markersize=22, label="Cámara que se mueve")),
                          (set(C) & Xl, dict(color="tab:green", markersize=22, label="Cámara que se queda")),
                          (Xl - set(C), dict(color="tab:blue", markersize=22, marker="^", label="Ubicación nueva"))],
                 figuras / "transicion.png", f"Transición lexicográfica, R = {R_BASE} m",
                 flechas=[(o, dd) for o, dd, _ in plan])
    Xb, fb = modelo.exacto(N, w, F, len(C), set(C), B_PLAN, solver=SOLVER)
    plan_b = transicion.hungaro(G, transicion.por_mover(C, Xb), sorted(Xb - set(C)))
    salidas.tabla(pd.DataFrame(plan_b, columns=["desde (nodo)", "hasta (nodo)", "metros"]), tablas, f"traslados_B{B_PLAN}")
    salidas.mapa(G, pol, [(set(C) - Xb, dict(color="tab:red", markersize=22, label="Cámara que se mueve")),
                          (set(C) & Xb, dict(color="tab:green", markersize=22, label="Cámara que se queda")),
                          (Xb - set(C), dict(color="tab:blue", markersize=22, marker="^", label="Ubicación nueva"))],
                 figuras / f"transicion_B{B_PLAN}.png",
                 f"Plan con {B_PLAN} traslados, R = {R_BASE} m: {100 * fb / total:.1f} % de cobertura",
                 flechas=[(o, dd) for o, dd, _ in plan_b])
    salidas.figura_curva(curva, inter[0], total, figuras / "curva.png")
    salidas.figura_sensibilidad(sens, figuras / "sensibilidad.png")
    salidas.mapa_html(G, [("Cruces semaforizados", F, "gray"), ("Cámaras hoy", set(C), "green"),
                          ("Configuración objetivo", Xl, "blue")], OUT / "mapa.html",
                      flechas=[(o, dd) for o, dd, _ in plan])

    # Observabilidad: cuántos flujos de tramo se deducen con las cámaras (supuestos en uoct/observabilidad.py)
    dentro = {x for x, dat in G.nodes(data=True) if pol.contains(Point(dat["x"], dat["y"]))}
    redes = {"toda la comuna": obs.red(G, dentro, dentro), "solo F (centroides)": obs.red(G, dentro, set(F) & dentro)}
    base = redes["toda la comuna"]
    Xe, _ = modelo.exacto(N, w, F, len(C), solver=SOLVER)
    curva_o = obs.voraz(base, F, len(F))
    Xo = [j for j, _ in curva_o[:len(C)]]
    configs = {"Actual": C, "Actual en línea": C_ef, "MCLP exacto": Xe, "MCLP lexicográfica": Xl,
               "Voraz de observabilidad": Xo, f"Todos los cruces F ({len(F)})": F}
    observ = pd.DataFrame([{"conservación": nr, "configuración": nc, **obs.resumen(r, X),
                            "cobertura %": 100 * cob.cobertura(X, N, w) / total}
                           for nr, r in redes.items() for nc, X in configs.items()])
    salidas.tabla(pd.DataFrame([{"conservación": nr, "tramos (arcos con un extremo en la comuna)": len(r.u),
                                 "nodos (con el virtual)": r.n, "rango completo": obs.rango_completo(r)}
                                for nr, r in redes.items()]), tablas, "red_observabilidad")
    salidas.tabla(observ.rename(columns={"cámaras": "cruces", "tramos medidos": "medidos",
                                         "tramos deducidos": "deducidos", "tramos identificables %": "tramos %",
                                         "km identificables %": "km %", "tránsito identificable %": "tránsito %"}),
                  tablas, "observabilidad")   # encabezados cortos: la tabla tiene que caber en el informe
    rc = obs.rango_completo(base)
    marcas = sorted(set(range(10, len(curva_o) + 1, 10)) | {len(C), len(curva_o)})
    salidas.tabla(pd.DataFrame([(b, 100 * curva_o[b - 1][1] / rc) for b in marcas], columns=["cámaras", "rango %"]),
                  tablas, "curva_observabilidad")
    salidas.figura_curva_obs(curva_o, rc, len(C), figuras / "curva_observabilidad.png")
    salidas.mapa_tramos(G, pol, base, C, figuras / "observabilidad_actual.png", "Flujos de tramo con las cámaras actuales")
    salidas.mapa_tramos(G, pol, base, Xo, figuras / "observabilidad_voraz.png",
                        f"Flujos de tramo con {len(Xo)} cámaras del voraz de observabilidad")

    print(control.to_string(index=False))
    print(sens.round(1).to_string(index=False))
    print(f"Caso base R = {R_BASE} m: f(C) = {100 * cob.cobertura(C, N, w) / total:.1f} %, "
          f"f* = {100 * fl / total:.1f} %, {len(plan)} traslados, {sum(m for *_, m in plan) / 1000:.1f} km; "
          f"se quedan {quedan} de {len(set(C))} ubicaciones")
    print(observ.round(1).to_string(index=False))
    print(f"Voraz de observabilidad: {len(curva_o)} cámaras llegan al {100 * curva_o[-1][1] / rc:.1f} % del rango "
          f"({rc} flujos independientes en {len(base.u)} tramos)")


if __name__ == "__main__":
    main()
