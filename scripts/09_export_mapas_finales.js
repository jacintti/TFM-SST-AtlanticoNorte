// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 09: EXPORTACIÓN DE MAPAS FINALES A GOOGLE DRIVE
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Exportar a Google Drive como archivos GeoTIFF y Shapefile los datos
//   necesarios para componer las cuatro figuras cartográficas centrales
//   del trabajo en QGIS:
//
//     · Figura 1 (Capítulo 3): Mapa de contexto y regiones de análisis.
//     · Figura 3 (Capítulo 5): Mapa de climatología OISST 1991-2020.
//     · Figura 5 (Capítulo 5): Mapa de salinidad media HYCOM 1993-2023.
//     · Figura 7 (Capítulo 5): Mapa de asimetría costera ERA5-Land
//                              (ventanas como polígonos vectoriales).
//
//   No se regenera el mapa de tendencia OISST píxel a píxel: ya está
//   exportado por el script 04c.
//
// USO:
//   1. Ejecutar el script con Run.
//   2. Ir al panel Tasks (junto a Console).
//   3. Para cada tarea pendiente, pulsar el botón Run individual.
//   4. Aceptar el diálogo de cada tarea.
//   5. Esperar a que terminen (varios minutos por tarea).
//   6. Los archivos aparecerán en Drive, carpeta TFM_SST_Mapas.
// ============================================================================


// ----------------------------------------------------------------------------
// 1. IMPORTAR MÓDULO DE GEOMETRÍAS
// ----------------------------------------------------------------------------
var geom = require('users/jacintod/TFM1:00_geometrias_module');


// ----------------------------------------------------------------------------
// 2. PARÁMETROS GENERALES
// ----------------------------------------------------------------------------
var CARPETA_DRIVE = 'TFM_SST_Mapas';

var SST_MIN = -3;
var SST_MAX = 35;


// ============================================================================
// EXPORTACIÓN 1 — GEOMETRÍAS COMO SHAPEFILE (mapa de contexto)
// ============================================================================
// Generamos una FeatureCollection con todas las geometrías del módulo 00,
// cada una etiquetada con su nombre y función. QGIS la abre directamente
// como capa vectorial y permite aplicar estilo personalizado.

var geometrias_lista = [
  {nombre: 'AN_limpio',  geom: geom.AN_limpio,  tipo: 'Dominio operativo'},
  {nombre: 'GSE',        geom: geom.GSE,        tipo: 'Regional'},
  {nombre: 'LAB',        geom: geom.LAB,        tipo: 'Regional'},
  {nombre: 'NOR',        geom: geom.NOR,        tipo: 'Regional'},
  {nombre: 'EUR',        geom: geom.EUR,        tipo: 'Regional'},
  {nombre: 'CBL',        geom: geom.CBL,        tipo: 'Regional (cold blob)'},
  {nombre: 'T53',        geom: geom.T53,        tipo: 'Transecto'},
  {nombre: 'L40w',       geom: geom.L40w,       tipo: 'Ventana costera'},
  {nombre: 'L40e',       geom: geom.L40e,       tipo: 'Ventana costera'},
  {nombre: 'L53w',       geom: geom.L53w,       tipo: 'Ventana costera'},
  {nombre: 'L53e',       geom: geom.L53e,       tipo: 'Ventana costera'},
  {nombre: 'L67w',       geom: geom.L67w,       tipo: 'Ventana costera'},
  {nombre: 'L67e',       geom: geom.L67e,       tipo: 'Ventana costera'}
];

var geometrias_fc = ee.FeatureCollection(geometrias_lista.map(function(g) {
  return ee.Feature(g.geom, {nombre: g.nombre, tipo: g.tipo});
}));

Export.table.toDrive({
  collection:    geometrias_fc,
  description:   'geometrias_TFM_contexto',
  folder:        CARPETA_DRIVE,
  fileNamePrefix:'geometrias_TFM_contexto',
  fileFormat:    'SHP'
});


// ============================================================================
// EXPORTACIÓN 2 — CLIMATOLOGÍA OISST 1991-2020
// ============================================================================
// SST media multianual de OISST sobre el periodo de referencia estándar
// definido por la OMM (1991-2020). Servirá como mapa base para situar
// gradientes térmicos.

var oisst_raw = ee.ImageCollection('NOAA/CDR/OISST/V2_1')
  .filterDate('1991-01-01', '2021-01-01')
  .filterBounds(geom.AN_limpio)
  .select('sst');

var oisst_diario = oisst_raw.map(function(img) {
  var sst = img.select('sst').multiply(0.01).rename('sst');
  var mascara = sst.gte(SST_MIN).and(sst.lte(SST_MAX));
  return sst.updateMask(mascara).clip(geom.AN_limpio);
});

var climatologia_oisst = oisst_diario.mean().rename('SST_media_C');

Export.image.toDrive({
  image:          climatologia_oisst,
  description:    'oisst_climatologia_1991_2020',
  folder:         CARPETA_DRIVE,
  fileNamePrefix: 'oisst_climatologia_1991_2020',
  region:         geom.AN_limpio,
  scale:          25000,
  crs:            'EPSG:4326',
  maxPixels:      1e9
});


// ============================================================================
// EXPORTACIÓN 3 — SALINIDAD HYCOM 1993-2023 (superficie)
// ============================================================================
// Salinidad media multianual de HYCOM en superficie (0 m), aplicando
// el desescalado correspondiente (× 0.001 + 20).

var hycom_raw = ee.ImageCollection('HYCOM/sea_temp_salinity')
  .filterDate('1993-01-01', '2024-01-01')
  .filterBounds(geom.AN_limpio)
  .select(['water_temp_0', 'salinity_0']);

var hycom_diario = hycom_raw.map(function(img) {
  var temp_sup = img.select('water_temp_0').multiply(0.001).add(20)
                    .rename('water_temp_0');
  var sal = img.select('salinity_0').multiply(0.001).add(20).rename('Sal_PSU');
  var mascara = temp_sup.gte(SST_MIN).and(temp_sup.lte(SST_MAX));
  return sal.updateMask(mascara).clip(geom.AN_limpio);
});

var salinidad_hycom = hycom_diario.mean().rename('Sal_PSU');

Export.image.toDrive({
  image:          salinidad_hycom,
  description:    'hycom_salinidad_superficie_1993_2023',
  folder:         CARPETA_DRIVE,
  fileNamePrefix: 'hycom_salinidad_superficie_1993_2023',
  region:         geom.AN_limpio,
  scale:          9000,
  crs:            'EPSG:4326',
  maxPixels:      1e9
});


// ============================================================================
// EXPORTACIÓN 4 — VENTANAS COSTERAS CON TENDENCIA ERA5-LAND
// ============================================================================
// Las seis ventanas costeras como polígonos vectoriales, cada uno con
// sus tendencias anual, invierno y verano como atributos. QGIS los
// pintará coloreados según el atributo elegido.

var resultados_era5 = ee.FeatureCollection([
  ee.Feature(geom.L40w, {
    ventana: 'L40w', costa: 'Oeste (America)', latitud: 40,
    tend_anual: 0.34, tend_invierno: 0.46, tend_verano: 0.29
  }),
  ee.Feature(geom.L40e, {
    ventana: 'L40e', costa: 'Este (Europa)',   latitud: 40,
    tend_anual: 0.30, tend_invierno: 0.20, tend_verano: 0.39
  }),
  ee.Feature(geom.L53w, {
    ventana: 'L53w', costa: 'Oeste (America)', latitud: 53,
    tend_anual: 0.58, tend_invierno: 0.62, tend_verano: 0.55
  }),
  ee.Feature(geom.L53e, {
    ventana: 'L53e', costa: 'Este (Europa)',   latitud: 53,
    tend_anual: 0.17, tend_invierno: 0.24, tend_verano: 0.12
  }),
  ee.Feature(geom.L67w, {
    ventana: 'L67w', costa: 'Oeste (America)', latitud: 67,
    tend_anual: 0.70, tend_invierno: 0.70, tend_verano: 0.38
  }),
  ee.Feature(geom.L67e, {
    ventana: 'L67e', costa: 'Este (Europa)',   latitud: 67,
    tend_anual: 0.38, tend_invierno: 0.35, tend_verano: 0.53
  })
]);

Export.table.toDrive({
  collection:     resultados_era5,
  description:    'era5land_ventanas_tendencias',
  folder:         CARPETA_DRIVE,
  fileNamePrefix: 'era5land_ventanas_tendencias',
  fileFormat:     'SHP'
});


// ----------------------------------------------------------------------------
// INFORMACIÓN EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 09 — EXPORTACIONES DE MAPAS FINALES PREPARADAS');
print('==================================================================');
print('Carpeta destino en Drive:', CARPETA_DRIVE);
print('');
print('Se han creado 4 tareas en el panel Tasks:');
print('  1. geometrias_TFM_contexto (SHP)');
print('  2. oisst_climatologia_1991_2020 (GeoTIFF)');
print('  3. hycom_salinidad_superficie_1993_2023 (GeoTIFF)');
print('  4. era5land_ventanas_tendencias (SHP)');
print('');
print('SIGUIENTE PASO:');
print('  - Ir al panel Tasks (junto a Console).');
print('  - Para cada tarea, pulsar el botón Run.');
print('  - Aceptar el diálogo (todos los parámetros están prerellenados).');
print('  - Esperar entre 5 y 15 minutos por tarea.');
print('  - Descargar los archivos desde Drive cuando estén Completed.');