// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 04c: NOAA OISST v2.1 — EXPORTACIÓN DE MAPAS DE TENDENCIA A DRIVE
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Exportar a Google Drive como archivos GeoTIFF los siguientes mapas:
//     · Tendencia decadal periodo completo (1982-2024).
//     · p-valor Mann-Kendall periodo completo.
//     · Tendencia decadal subperiodo 1 (1982-2000).
//     · Tendencia decadal subperiodo 2 (2001-2024).
//
//   Una vez exportados, los archivos se descargan desde Drive y se abren
//   en QGIS o ArcGIS para componer las figuras.
//
// FLUJO DE USO:
//   1. Ejecutar el script con Run.
//   2. Ir al panel Tasks (pestaña a la derecha de Console).
//   3. Para cada tarea, pulsar el botón Run individual.
//   4. Esperar (5-30 min). Las tareas terminan en estado Completed.
//   5. Los archivos aparecen en tu Google Drive, carpeta TFM_SST_Exports.
// ============================================================================


// ----------------------------------------------------------------------------
// 1. IMPORTAR MÓDULO DE GEOMETRÍAS
// ----------------------------------------------------------------------------
var geom = require('users/jacintod/TFM1:00_geometrias_module');


// ----------------------------------------------------------------------------
// 2. PARÁMETROS
// ----------------------------------------------------------------------------
var FECHA_INICIO = '1981-09-01';
var FECHA_FIN    = '2025-01-01';
var SST_MIN = -3;
var SST_MAX = 35;

var ANIO_INICIO_SERIE = 1982;
var ANIO_FIN_SERIE    = 2024;
var ANIO_CORTE = 2000;

var CARPETA_DRIVE = 'TFM_SST_Exports';
var ESCALA_EXPORT = 25000;  // 25 km, resolución nativa OISST


// ----------------------------------------------------------------------------
// 3. RECONSTRUIR COLECCIÓN MENSUAL Y ANUAL (idéntico a script 04)
// ----------------------------------------------------------------------------
var oisst_raw = ee.ImageCollection('NOAA/CDR/OISST/V2_1')
  .filterDate(FECHA_INICIO, FECHA_FIN)
  .filterBounds(geom.AN_limpio)
  .select('sst');

var oisst_diario = oisst_raw.map(function(img) {
  var sst = img.select('sst').multiply(0.01).rename('sst');
  var mascara = sst.gte(SST_MIN).and(sst.lte(SST_MAX));
  return sst.updateMask(mascara)
            .clip(geom.AN_limpio)
            .copyProperties(img, ['system:time_start']);
});

var anio_inicio_calc = ee.Number(ee.Date(FECHA_INICIO).get('year'));
var anio_fin_calc    = ee.Number(ee.Date(FECHA_FIN).get('year'));
var meses_total = anio_fin_calc.subtract(anio_inicio_calc).multiply(12).toInt();

var oisst_mensual = ee.ImageCollection.fromImages(
  ee.List.sequence(0, meses_total.subtract(1)).map(function(offset) {
    var fecha_mes = ee.Date(FECHA_INICIO).advance(ee.Number(offset), 'month');
    var fecha_fin_mes = fecha_mes.advance(1, 'month');
    var imgs_del_mes = oisst_diario.filterDate(fecha_mes, fecha_fin_mes);
    var num_imagenes = imgs_del_mes.size();
    return imgs_del_mes.median()
      .set('system:time_start', fecha_mes.millis())
      .set('anio', fecha_mes.get('year'))
      .set('mes', fecha_mes.get('month'))
      .set('n_imagenes', num_imagenes);
  })
)
.filter(ee.Filter.gt('n_imagenes', 0));

// Serie anual completa
var anios_serie = ee.List.sequence(ANIO_INICIO_SERIE, ANIO_FIN_SERIE);
var oisst_anual = ee.ImageCollection.fromImages(
  anios_serie.map(function(a) {
    var meses_del_anio = oisst_mensual.filter(ee.Filter.eq('anio', a));
    return meses_del_anio.mean()
      .set('anio', a)
      .addBands(ee.Image.constant(a).toFloat().rename('year'))
      .set('system:time_start', ee.Date.fromYMD(a, 7, 1).millis());
  })
);


// ----------------------------------------------------------------------------
// 4. FUNCIÓN PARA CALCULAR MAPA DE TENDENCIA EN UN PERIODO
// ----------------------------------------------------------------------------
var calcular_mapa_tendencia = function(coleccion_anual, anio_min, anio_max) {
  var subset = coleccion_anual
    .filter(ee.Filter.gte('anio', anio_min))
    .filter(ee.Filter.lte('anio', anio_max))
    .select(['year', 'sst']);

  var tendencia = subset
    .reduce(ee.Reducer.sensSlope())
    .select('slope')
    .multiply(10)
    .rename('tendencia_C_por_decada');

  return tendencia;
};

// Mapa de p-valor píxel a píxel (solo lo calculamos para el periodo completo,
// que es el más informativo)
var calcular_mapa_pvalor = function(coleccion_anual) {
  // Reducer kendallsCorrelation(1) sobre la banda SST ordenada cronológicamente
  // devuelve tau y p-value para cada píxel.
  var resultado = coleccion_anual.select('sst')
    .reduce(ee.Reducer.kendallsCorrelation(1));
  // Bandas devueltas: 'sst_tau' y 'sst_p-value'
  return resultado.select('sst_p-value').rename('p_valor');
};


// ----------------------------------------------------------------------------
// 5. CALCULAR LOS CUATRO MAPAS
// ----------------------------------------------------------------------------
var mapa_tendencia_completo = calcular_mapa_tendencia(oisst_anual, 
                                                       ANIO_INICIO_SERIE, 
                                                       ANIO_FIN_SERIE);

var mapa_pvalor_completo = calcular_mapa_pvalor(oisst_anual);

var mapa_tendencia_sub1 = calcular_mapa_tendencia(oisst_anual, 
                                                   ANIO_INICIO_SERIE, 
                                                   ANIO_CORTE);

var mapa_tendencia_sub2 = calcular_mapa_tendencia(oisst_anual, 
                                                   ANIO_CORTE + 1, 
                                                   ANIO_FIN_SERIE);


// ----------------------------------------------------------------------------
// 6. EXPORTACIONES A GOOGLE DRIVE
// ----------------------------------------------------------------------------
// Cada Export.image.toDrive crea una tarea en el panel Tasks.
// Hay que activarla manualmente desde el panel.

Export.image.toDrive({
  image: mapa_tendencia_completo,
  description: 'oisst_tendencia_AN_1982-2024',
  folder: CARPETA_DRIVE,
  fileNamePrefix: 'oisst_tendencia_AN_1982-2024',
  region: geom.AN_limpio,
  scale: ESCALA_EXPORT,
  crs: 'EPSG:4326',
  maxPixels: 1e9
});

Export.image.toDrive({
  image: mapa_pvalor_completo,
  description: 'oisst_pvalor_AN_1982-2024',
  folder: CARPETA_DRIVE,
  fileNamePrefix: 'oisst_pvalor_AN_1982-2024',
  region: geom.AN_limpio,
  scale: ESCALA_EXPORT,
  crs: 'EPSG:4326',
  maxPixels: 1e9
});

Export.image.toDrive({
  image: mapa_tendencia_sub1,
  description: 'oisst_tendencia_AN_1982-2000',
  folder: CARPETA_DRIVE,
  fileNamePrefix: 'oisst_tendencia_AN_1982-2000',
  region: geom.AN_limpio,
  scale: ESCALA_EXPORT,
  crs: 'EPSG:4326',
  maxPixels: 1e9
});

Export.image.toDrive({
  image: mapa_tendencia_sub2,
  description: 'oisst_tendencia_AN_2001-2024',
  folder: CARPETA_DRIVE,
  fileNamePrefix: 'oisst_tendencia_AN_2001-2024',
  region: geom.AN_limpio,
  scale: ESCALA_EXPORT,
  crs: 'EPSG:4326',
  maxPixels: 1e9
});


// ----------------------------------------------------------------------------
// 7. INFORMACIÓN EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 04c — EXPORTACIONES A DRIVE PREPARADAS');
print('==================================================================');
print('Se han creado 4 tareas en el panel Tasks:');
print('  1. oisst_tendencia_AN_1982-2024 (tendencia periodo completo)');
print('  2. oisst_pvalor_AN_1982-2024 (significancia Mann-Kendall)');
print('  3. oisst_tendencia_AN_1982-2000 (tendencia subperiodo 1)');
print('  4. oisst_tendencia_AN_2001-2024 (tendencia subperiodo 2)');
print('');
print('SIGUIENTE PASO:');
print('  - Ir al panel Tasks (pestaña junto a Console).');
print('  - Para cada tarea, pulsar el botón Run individual.');
print('  - Esperar a que termine cada exportación (puede tardar varios');
print('    minutos por tarea).');
print('  - Los archivos GeoTIFF aparecerán en Drive, carpeta '+CARPETA_DRIVE);
