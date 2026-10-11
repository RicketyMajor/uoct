import pandas as pd
from uoct.datos import camaras_comuna, leer_inventario

ENC = ["NOMBRE / UBICACIÓN", "COMUNA", "ESTADO", "MARCA", "MODELO", "MARCA", "MODELO", "LATITUD", "LONGITUD"]
SUPER = [None, None, None, "DATOS CÁMARA", None, "DATOS CODIFICADOR DE VIDEO", None, None, None]


def inventario_falso(tmp_path):
    filas = [SUPER, ENC,
             ["A - B", "Santiago", "En línea", "X", "Y", None, None, "-33.45", "-70.65"],
             SUPER, ENC,                                    # encabezado repetido a mitad de hoja
             ["C - D", "Santiago ", "Operativa", "X", "Y", None, None, -33.44, "-70,66"],
             ["E - F", "Buin", "Fuera de línea", "X", "Y", None, None, -33.70, -70.70],
             ["G - H", "Santiago", "Fuera de línea", "X", "Y", None, None, None, None]]
    ruta = tmp_path / "inv.xlsx"
    pd.DataFrame(filas).to_excel(ruta, sheet_name="RM", header=False, index=False)
    return ruta


def test_leer_inventario_limpia_encabezados_y_coordenadas(tmp_path):
    df = leer_inventario(inventario_falso(tmp_path))
    assert list(df["nombre"]) == ["A - B", "C - D", "E - F"]          # sin encabezados ni filas sin coordenadas
    assert df.loc[1, "lon"] == -70.66                                  # coma decimal
    assert list(df["en_linea"]) == [True, True, False]                 # «Operativa» cuenta como en línea


def test_camaras_comuna_filtra_y_proyecta(tmp_path):
    g = camaras_comuna(leer_inventario(inventario_falso(tmp_path)))
    assert len(g) == 2 and g.crs.to_epsg() == 32719
    assert 340_000 < g.geometry.x.iloc[0] < 360_000                    # Santiago en UTM 19S
