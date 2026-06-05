// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// MÓDULO 00: GEOMETRÍAS DEL ÁREA DE ESTUDIO
// ----------------------------------------------------------------------------
// Este script NO produce salidas visibles. Define las 13 geometrías del área
// de estudio del trabajo y las exporta para que otros scripts las importen
// mediante require(). Cualquier ajuste futuro de coordenadas se hace solo en
// este archivo y se propaga automáticamente al resto.
//
// CONTENIDO:
//   1. Polígono de contexto general (AN)
//   2. Polígono de exclusión (MED)
//   3. Polígonos de análisis regional (GSE, LAB, NOR, EUR)
//   4. Transecto oceánico (T53)
//   5. Ventanas costeras terrestres (L40w, L40e, L53w, L53e, L67w, L67e)
//   6. Geometría operativa AN_limpio = AN.difference(MED)
//
// USO DESDE OTRO SCRIPT:
//   var geom = require('users/[TU_USUARIO_GEE]/TFM_SST:00_geometrias_module');
//   var poligono = geom.GSE;             // accede al polígono Gulf Stream
//   var costera  = geom.L40e;            // accede a una ventana costera
//   var todas    = geom.lista_completa;  // accede a la lista completa
//
// VERSIÓN: v5 (modularizada, basada en v4 validada visualmente)
// FECHA:   2026
// ============================================================================


// ----------------------------------------------------------------------------
// 1. POLÍGONO DE CONTEXTO GENERAL
// ----------------------------------------------------------------------------
exports.AN = ee.Geometry.Rectangle([-80, 20, 20, 75]);
//                                  lon_oeste, lat_sur, lon_este, lat_norte
//   Atlántico Norte. Cubre desde el Caribe norte hasta el borde ártico.
//   Se utiliza siempre con la exclusión MED aplicada (AN_limpio).


// ----------------------------------------------------------------------------
// 2. POLÍGONO DE EXCLUSIÓN
// ----------------------------------------------------------------------------
//exports.MED = ee.Geometry.Rectangle([-6, 30, 20, 45]);
exports.MED = ee.Geometry.Polygon([
  [
    [-6, 30],   // SW
    [20, 30],   // SE
    [20, 45.9],   // NE
    [3, 45.9],    // Costa sur de Francia
    [0, 44],
    [-2, 43],
    [-4, 42],
    [-6, 40],   // Baja por el Atlántico portugués
    [-6, 30]    // Cierre
  ]
]);

// ----------------------------------------------------------------------------
// 3. POLÍGONOS DE ANÁLISIS REGIONAL
// ----------------------------------------------------------------------------
exports.GSE = ee.Geometry.Rectangle([-75, 36, -50, 44]);
//   Gulf Stream Extension. Zona de calentamiento intensificado.
//   Referencia: Wu et al. (2012), Todd & Ren (2023).

exports.LAB = ee.Geometry.Rectangle([-60, 53, -45, 60]);
//   Mar de Labrador. Convección profunda y propagación de la GSA reciente.
//   Referencia: Biló et al. (2022), Asbjørnsen et al. (2024).

exports.NOR = ee.Geometry.Rectangle([-20, 65, 20, 75]);
//   Mares Nórdicos. Amplificación ártica y formación de aguas densas.
//   Referencia: Rantanen et al. (2022), Isaksen et al. (2022).

exports.EUR = ee.Geometry.Rectangle([-12, 36, -6, 48]);
//   Costa atlántica europea. Fachada oceánica de la influencia climática
//   de la Corriente del Golfo sobre Europa occidental.
//   Borde oriental ajustado a -6° para no solapar con MED.
//   Referencia: Biguino et al. (2023), Piedracoba et al. (2024).

exports.CBL = ee.Geometry.Rectangle([-45, 45, -15, 60]);
//   Cold Blob / warming hole del Atlántico subpolar. Captura específicamente
//   la zona de enfriamiento relativo identificada por Rahmstorf et al. (2015)
//   y Caesar et al. (2018) como huella térmica del debilitamiento de la AMOC.
//   Solapa parcialmente con LAB (Mar de Labrador) por diseño, ya que ambas
//   regiones representan aspectos complementarios del mismo fenómeno.


// ----------------------------------------------------------------------------
// 4. TRANSECTO OCEÁNICO (doble uso: promedio regional + perfil longitudinal)
// ----------------------------------------------------------------------------
exports.T53 = ee.Geometry.Rectangle([-55, 52.5, -10, 53.5]);
//   Franja oceánica a 53°N entre Terranova e Irlanda. Tramo central de la
//   Corriente Noratlántica. Se reduce como región (media) o como perfil
//   longitudinal (gradiente este-oeste).


// ----------------------------------------------------------------------------
// 5. VENTANAS COSTERAS TERRESTRES
// ----------------------------------------------------------------------------
//   Pares de ventanas en tres latitudes (40°, 53°, 67°N), una a cada lado
//   del Atlántico. Tamaño: 2°lat × 3°lon ≈ 220 × 240 km a estas latitudes.
//   Función: extraer temperatura del aire (ERA5-Land) sobre tierra.

// Costa americana                        // Costa europea
exports.L40w = ee.Geometry.Rectangle([-76, 39, -73, 41]);  // New Jersey
exports.L40e = ee.Geometry.Rectangle([ -9, 39,  -6, 41]);  // Portugal central
exports.L53w = ee.Geometry.Rectangle([-58, 52, -55, 54]);  // Terranova
exports.L53e = ee.Geometry.Rectangle([-10, 52,  -7, 54]);  // Irlanda
exports.L67w = ee.Geometry.Rectangle([-68, 66, -64, 68]);  // Baffin
exports.L67e = ee.Geometry.Rectangle([ 12, 66,  15, 68]);  // Bodø-Lofoten


// ----------------------------------------------------------------------------
// 6. GEOMETRÍA OPERATIVA: AN sin Mediterráneo
// ----------------------------------------------------------------------------
exports.AN_limpio = exports.AN.difference(exports.MED);
//   Dominio efectivo para cálculos oceánicos. Excluye el Mediterráneo
//   pero conserva Mar del Norte, Báltico y mares Nórdicos (sí oceánicos).


// ----------------------------------------------------------------------------
// 7. LISTAS AGRUPADAS (opcionales, para bucles sobre geometrías)
// ----------------------------------------------------------------------------

// Polígonos de análisis regional (5: incluye AN_limpio + 4 regionales)
exports.lista_regionales = [
  {nombre: 'AN_limpio', geometria: exports.AN_limpio},
  {nombre: 'GSE', geometria: exports.GSE},
  {nombre: 'LAB', geometria: exports.LAB},
  {nombre: 'NOR', geometria: exports.NOR},
  {nombre: 'EUR', geometria: exports.EUR},
  {nombre: 'CBL', geometria: exports.CBL}
];

// Ventanas costeras (6)
exports.lista_costeras = [
  {nombre: 'L40w', geometria: exports.L40w},
  {nombre: 'L40e', geometria: exports.L40e},
  {nombre: 'L53w', geometria: exports.L53w},
  {nombre: 'L53e', geometria: exports.L53e},
  {nombre: 'L67w', geometria: exports.L67w},
  {nombre: 'L67e', geometria: exports.L67e}
];

// Lista completa: las 12 geometrías sobre las que se calculan series
// (AN_limpio + 4 regionales + T53 + 6 costeras)
exports.lista_completa = exports.lista_regionales.concat([
  {nombre: 'T53', geometria: exports.T53}
]).concat(exports.lista_costeras);


// ----------------------------------------------------------------------------
// 8. INFORMACIÓN DIAGNÓSTICA (solo si se ejecuta este script directamente)
// ----------------------------------------------------------------------------
// Comentar las siguientes líneas en producción si no quieres salida visible
// al importar este módulo desde otros scripts. En GEE, los exports se
// procesan sin ejecutar los prints si el módulo se importa con require().

print('==================================================================');
print('MÓDULO 00 — GEOMETRÍAS DEL TFM');
print('==================================================================');
print('Geometrías exportadas: 14 (1 contexto + 1 exclusión + 5 regionales');
print('                        + 1 transecto + 6 costeras)');
print('Geometría operativa principal: AN_limpio = AN.difference(MED)');
print('==================================================================');
