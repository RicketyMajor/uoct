"""Inventario CCTV de la UOCT → cámaras de una comuna, en UTM 19S."""
import geopandas as gpd
import pandas as pd

CRS_UTM = 32719  # UTM 19S, el de Santiago
COLUMNAS = ["nombre", "comuna", "estado", "marca", "modelo", "cod_marca", "cod_modelo", "lat", "lon"]


def leer_inventario(ruta):
    """Hoja «RM». Tiene dos filas de encabezado que se repiten a mitad de hoja y coordenadas
    a veces como texto (con punto o coma decimal). «Operativa» se cuenta como en línea."""
    df = pd.read_excel(ruta, sheet_name="RM", header=None, skiprows=2, names=COLUMNAS)
    df = df[df["nombre"].notna() & (df["nombre"] != "NOMBRE / UBICACIÓN")]
    for c in ("lat", "lon"):
        df[c] = pd.to_numeric(df[c].astype(str).str.strip().str.replace(",", "."), errors="coerce")
    df = df.dropna(subset=["lat", "lon"]).copy()
    df["en_linea"] = df["estado"].astype(str).str.strip().isin(["En línea", "Operativa"])
    return df.reset_index(drop=True)


def camaras_comuna(df, comuna="Santiago"):
    sel = df[df["comuna"].astype(str).str.strip().str.casefold() == comuna.casefold()]
    g = gpd.GeoDataFrame(sel, geometry=gpd.points_from_xy(sel["lon"], sel["lat"]), crs=4326)
    return g.to_crs(CRS_UTM).reset_index(drop=True)
