// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 07a: HYCOM — CARGA, PREPROCESAMIENTO Y VERIFICACIÓN
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Cargar el producto HYCOM (Hybrid Coordinate Ocean Model) del catálogo
//   de GEE. Extraer temperatura y salinidad en cuatro niveles de
//   profundidad: 0, 100, 500 y 1000 metros. Aplicar el desescalado
//   apropiado, máscara térmica, y agregar a escala mensual.
//
//   HYCOM es un modelo OCEÁNICO ASIMILATIVO global. A diferencia de OISST
//   (interpolación de observaciones in situ + satélite) y MODIS (medida
//   directa por satélite), HYCOM aporta TRES novedades fundamentales:
//     · Datos a múltiples profundidades (no solo superficie).
//     · Salinidad además de temperatura.
//     · Cobertura sin huecos (al ser modelo, no depende de nubosidad).
//
// COBERTURA TEMPORAL: 1992-10-02 a 2024-09-05.
// AÑOS COMPLETOS ANALIZABLES: 1993-2023 (31 años).
// RESOLUCIÓN ESPACIAL: 0,08° (≈9 km).
//
// DESESCALADO:
//   · Temperatura: valor_real = valor_bruto × 0,001 + 20 (°C)
//   · Salinidad:   valor_real = valor_bruto × 0,001 + 20 (PSU)
//
// VERIFICACIÓN ESPERADA:
//   · Temperatura HYCOM 0 m, junio 2010, AN_limpio: ~12,2-12,5 °C
//     (coherente con OISST y MODIS para el mismo mes/región).
//   · Salinidad HYCOM 0 m media: ~34,5-35 PSU sobre Atlántico abierto.
//
// REFERENCIA DEL PRODUCTO:
//   Chassignet, E.P., Hurlburt, H.E., Smedstad, O.M., et al. (2007).
//   The HYCOM (HYbrid Coordinate Ocean Model) data assimilative system.
//   Journal of Marine Systems, 65(1-4), 60-83.
//   https://doi.org/10.1016/j.jmarsys.2005.09.016
// ============================================================================


// ----------------------------------------------------------------------------
// 1. IMPORTAR MÓDULO DE GEOMETRÍAS
// ----------------------------------------------------------------------------
var geom = require('users/jacintod/TFM1:00_geometrias_module');


// ----------------------------------------------------------------------------
// 2. PARÁMETROS DEL SCRIPT
// ----------------------------------------------------------------------------
var FECHA_INICIO = '1992-10-02';   // primer día disponible HYCOM
var FECHA_FIN    = '2024-09-06';   // un día después del último

var SST_MIN = -3;    // máscara térmica para temperatura
var SST_MAX = 35;

// Bandas que nos interesan
var BANDAS_TEMP = ['water_temp_0', 'water_temp_100', 'water_temp_500', 'water_temp_1000'];
var BANDAS_SAL  = ['salinity_0',   'salinity_100',   'salinity_500',   'salinity_1000'];
var TODAS_BANDAS = BANDAS_TEMP.concat(BANDAS_SAL);

// Paletas
var PALETA_SST = [
  '040274','040281','0502a3','0502b8','0502ce','0502e6','0602ff',
  '235cb1','307ef3','269db1','30c8e2','32d3ef','3be285','3ff38f',
  '86e26f','3ae237','b5e22e','d6e21f','fff705','ffd611','ffb613',
  '17f80e','ff8b13','ff6e08','ff500d','ff0000','de0101','c21301','a71001','911003'
];
// Paleta divergente para salinidad: marrón (dulce) → blanco (35 PSU) → azul (salado)
var PALETA_SAL = [
  '8c510a','bf812d','dfc27d','f6e8c3','f5f5f5',
  'c7eae5','80cdc1','35978f','01665e'
];


// ----------------------------------------------------------------------------
// 3. CARGAR PRODUCTO Y FILTRAR
// ----------------------------------------------------------------------------
var hycom_raw = ee.ImageCollection('HYCOM/sea_temp_salinity')
  .filterDate(FECHA_INICIO, FECHA_FIN)
  .filterBounds(geom.AN_limpio)
  .select(TODAS_BANDAS);


// ----------------------------------------------------------------------------
// 4. PREPROCESAMIENTO: DESESCALADO + MÁSCARA TÉRMICA + RECORTE
// ----------------------------------------------------------------------------
var hycom_diario = hycom_raw.map(function(img) {
  // Desescalar todas las bandas (factor 0,001 y offset 20)
  var temp = img.select(BANDAS_TEMP).multiply(0.001).add(20);
  var sal  = img.select(BANDAS_SAL).multiply(0.001).add(20);

  // Renombrar para mantener nombres claros
  temp = temp.rename(BANDAS_TEMP);
  sal  = sal.rename(BANDAS_SAL);

  // Máscara térmica solo a la temperatura de superficie (water_temp_0)
  // (las temperaturas profundas pueden estar fuera del rango 0-35 °C
  // y son físicamente válidas)
  var temp_sup = temp.select('water_temp_0');
  var mascara = temp_sup.gte(SST_MIN).and(temp_sup.lte(SST_MAX));

  return temp.addBands(sal)
             .updateMask(mascara)
             .clip(geom.AN_limpio)
             .copyProperties(img, ['system:time_start']);
});


// ----------------------------------------------------------------------------
// 5. AGREGACIÓN MENSUAL (MEDIANA)
// ----------------------------------------------------------------------------
var anio_inicio = ee.Number(ee.Date(FECHA_INICIO).get('year'));
var anio_fin    = ee.Number(ee.Date(FECHA_FIN).get('year'));
var meses_total = anio_fin.subtract(anio_inicio).multiply(12).toInt();

var hycom_mensual = ee.ImageCollection.fromImages(
  ee.List.sequence(0, meses_total.subtract(1)).map(function(offset) {
    var fecha_mes = ee.Date(FECHA_INICIO).advance(ee.Number(offset), 'month');
    var anio = fecha_mes.get('year');
    var mes  = fecha_mes.get('month');
    var fi = ee.Date.fromYMD(anio, mes, 1);
    var ff = fi.advance(1, 'month');
    var imgs = hycom_diario.filterDate(fi, ff);
    var n = imgs.size();
    return imgs.median()
      .set('system:time_start', fi.millis())
      .set('anio', anio)
      .set('mes', mes)
      .set('n_imagenes', n);
  })
).filter(ee.Filter.gt('n_imagenes', 0));


// ----------------------------------------------------------------------------
// 6. INFORMACIÓN DIAGNÓSTICA EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 07a — HYCOM CARGA Y PREPROCESAMIENTO');
print('==================================================================');
print('Periodo:', FECHA_INICIO, '→', FECHA_FIN);
print('Bandas extraídas (temperatura):', BANDAS_TEMP);
print('Bandas extraídas (salinidad):', BANDAS_SAL);
print('');
print('Imágenes diarias cargadas:', hycom_diario.size());
print('Imágenes mensuales generadas:', hycom_mensual.size());
print('Primer mes con dato:', ee.Date(
  hycom_mensual.aggregate_min('system:time_start')).format('YYYY-MM'));
print('Último mes con dato:', ee.Date(
  hycom_mensual.aggregate_max('system:time_start')).format('YYYY-MM'));

// Verificación 1: temperatura HYCOM 0 m, junio 2010
var junio_2010 = hycom_mensual.filterDate('2010-06-01', '2010-07-01').first();
var temp_jun2010 = junio_2010.select('water_temp_0').reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: geom.AN_limpio,
  scale: 25000,
  maxPixels: 1e9,
  bestEffort: true
});
print('Temp. HYCOM 0m junio 2010 sobre AN_limpio (esperado ~12-13°C):',
  temp_jun2010);

// Verificación 2: salinidad HYCOM 0 m, junio 2010
var sal_jun2010 = junio_2010.select('salinity_0').reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: geom.AN_limpio,
  scale: 25000,
  maxPixels: 1e9,
  bestEffort: true
});
print('Salinidad HYCOM 0m junio 2010 sobre AN_limpio (esperado ~34,5-35 PSU):',
  sal_jun2010);


// ----------------------------------------------------------------------------
// 7. VISUALIZACIÓN: MEDIAS CLIMATOLÓGICAS
// ----------------------------------------------------------------------------
Map.centerObject(geom.AN_limpio, 3);
Map.setOptions('HYBRID');

// Mapa 1: temperatura superficie media del periodo
var temp0_media = hycom_mensual.select('water_temp_0').mean();
Map.addLayer(
  temp0_media,
  {min: 0, max: 28, palette: PALETA_SST},
  '01. HYCOM — Temperatura 0m media 1992-2024 (°C)',
  true
);

// Mapa 2: salinidad superficie media del periodo (la novedad)
var sal0_media = hycom_mensual.select('salinity_0').mean();
Map.addLayer(
  sal0_media,
  {min: 33, max: 37, palette: PALETA_SAL},
  '02. HYCOM — Salinidad 0m media 1992-2024 (PSU)',
  false
);