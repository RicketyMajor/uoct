# UOCT — Redistribución de cámaras CCTV (Grafos y Algoritmos, UDP)

Prototipo del proyecto semestral: ¿dónde deberían estar las cámaras de la UOCT en la comuna de Santiago y cómo se llega ahí desde donde están hoy?

## Correr

    python -m venv .venv
    .venv/Scripts/python -m pip install -r requirements.txt
    .venv/Scripts/python -m pytest
    .venv/Scripts/python run.py

`run.py` escribe en `out/` (no va al repo):

- `tablas/` — cada tabla en `.csv` y en `.tex` (para `\input` en el informe): `control` (contraste con el Estado del proyecto), `sensibilidad` (radio, peso y sentido), `curva` (cobertura óptima con a lo más B traslados), `traslados` (plan de la transición lexicográfica), `traslados_B10` (plan con 10 traslados), `intercambio` (voraz de a una cámara según el umbral δ), `observabilidad` (cuánto se deduce con cada configuración, con y sin centroides), `curva_observabilidad` (rango del voraz de observabilidad según las cámaras).
- `figuras/` — `estado.png`, `transicion.png`, `transicion_B10.png`, `curva.png`, `sensibilidad.png`, `curva_observabilidad.png`, `observabilidad_actual.png` y `observabilidad_voraz.png` (tramos medidos, deducidos y sin información).
- `mapa.html` — mapa interactivo (cruces, cámaras de hoy, configuración objetivo y traslados).

`python run.py --descargar` vuelve a bajar la red y los semáforos de OpenStreetMap. Sin esa opción se usa la copia de `data/`, para que todo el grupo corra sobre la misma red.

## Qué hace cada módulo

| Módulo | Qué hace |
|---|---|
| `uoct/datos.py` | Lee el inventario CCTV (hoja «RM») y filtra las cámaras de la comuna, en UTM 19S (EPSG:32719) |
| `uoct/grafo.py` | Red para autos de OSM con 500 m de margen; F = intersecciones a ≤ 50 m de un semáforo; cada cámara a su cruce de F más cercano |
| `uoct/cobertura.py` | Distancias por la red (no dirigida, o dirigida aguas arriba), vecindades N_i para un radio R, pesos w_i y la cobertura f(X) |
| `uoct/modelo.py` | MCLP con p cámaras: voraz, exacto (Pyomo + HiGHS), lexicográfico en dos fases y curva f*(B) |
| `uoct/transicion.py` | Asignación húngara (qué cámara va a dónde) y voraz de a una cámara |
| `uoct/observabilidad.py` | Cuántos flujos de tramo se deducen con las cámaras: red de conservación, tramos identificables (puentes), rango y voraz perezoso |
| `uoct/salidas.py` | Tablas, mapas, gráficos y mapa HTML |

## Parámetros (arriba de `run.py`)

- `R_BASE = 300` m y `RADIOS` para la sensibilidad. 300 m es el radio de Dell'Olmo et al. (2014) y de Wheeler (2026), los dos MCLP de cámaras más cercanos a este problema.
- `DELTAS`: ganancia neta mínima (en % de la cobertura total) para que el voraz de a una cámara haga un traslado.
- `B_PLAN = 10`: traslados del plan que se dibuja aparte.
- `SOLVER = "highs"`; con licencia de Gurobi, `"gurobi"`.

## Datos y supuestos

- **Inventario CCTV de la UOCT** al 31-08-2026 (`data/`): 1.037 cámaras en la RM, 497 con coordenadas, 58 rotuladas «Santiago». Dos caen a ~1 km fuera de la comuna y se descartan: quedan 56, 43 en línea («Operativa» cuenta como en línea).
- **OpenStreetMap** (descarga del 2026-10-11): red `drive` de la comuna con 500 m de margen y los nodos `highway=traffic_signals` de la comuna.
- **Peso w_i provisional** hasta tener el flujo de la UOCT: la mitad de la suma del tránsito diario medio de la clase de cada calle que llega al cruce. Los valores por clase (autopista o troncal 46.632; primaria 13.022; secundaria o terciaria 5.852; el resto 1.392 veh/día) son los medidos por Shen et al. (2024, *Atmospheric Environment* 335, 120719, Tabla 3) en 10.698 puntos de conteo europeos. La clase sola explica alrededor del 44 % de la varianza del flujo: es un proxy, no una medición.
- Las cámaras fuera de línea se pueden reubicar (supuesto S5 del Informe 1); la línea base se mide con todas y con las en línea.
- **Observabilidad** (Li et al., 2023; Zhuo et al., 2024), con supuestos que la UOCT no ha confirmado:
  - una cámara cuenta sin error el flujo de todos los tramos de su cruce (cada sentido de una calle es un tramo);
  - el flujo se conserva en todos los nodos de la comuna, y el exterior se contrae a un nodo virtual. En la variante «solo F», los nodos que no son cruces semaforizados se tratan como centroides (generan o atraen viajes) y también se contraen al nodo virtual;
  - un tramo sin cámara se deduce si y solo si es un puente de la red sin los tramos medidos. Los tests lo comprueban contra el rango de la matriz de conservación en 300 redes al azar.
