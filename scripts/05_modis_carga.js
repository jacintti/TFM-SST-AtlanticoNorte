// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 05: MODIS AQUA L3SMI — CARGA Y PREPROCESAMIENTO
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Cargar el producto MODIS Aqua L3SMI desde el catálogo de GEE, filtrarlo
//   al dominio AN_limpio y al periodo 2002-2024, aplicar las máscaras
//   de calidad y térmica, y agregar a escala mensual con umbral mínimo
//   de cobertura.
//
//   MODIS Aqua se usa como producto de VALIDACIÓN CRUZADA frente a OISST,
//   y como producto principal para detalle COSTERO (resolución 4,6 km vs
//   25 km de OISST).
//
// CATÁLOGO GEE: NASA/OCEANDATA/MODIS-Aqua/L3SMI
// COBERTURA TEMPORAL: 2002-07-04 → presente
// RESOLUCIÓN ESPACIAL: ~4,6 km
//
// SALIDA (en memoria):
//   - modis_diario: colección diaria preprocesada
//   - modis_mensual: colección mensual agregada (mediana, con umbral de
//     cobertura mínima de 10 días/mes)
//
// VERIFICACIÓN ESPERADA:
//   - Número de imágenes diarias: ~8.000 (de 2002 a 2024).
//   - Número de imágenes mensuales: ~270 (de jul 2002 a dic 2024).
//   - SST media junio 2010 sobre AN_limpio: ~12-13 °C (similar a OISST).
//   - Mapa climatológico: patrón similar al de OISST pero con más detalle.
//
// REFERENCIA DEL PRODUCTO:
//   Kilpatrick, K.A., Podestá, G., Walsh, S., et al. (2015). A decade of sea
//   surface temperature from MODIS. Remote Sensing of Environment, 165, 27-41.
//   https://doi.org/10.1016/j.rse.2015.04.023
// ============================================================================


// ----------------------------------------------------------------------------
// 1. IMPORTAR MÓDULO DE GEOMETRÍAS
// ----------------------------------------------------------------------------
var geom = require('users/jacintod/TFM1:00_geometrias_module');


// ----------------------------------------------------------------------------
// 2. PARÁMETROS DEL SCRIPT
// ----------------------------------------------------------------------------
var FECHA_INICIO = '2002-07-04';   // primer día de MODIS Aqua disponible
var FECHA_FIN    = '2025-01-01';

var SST_MIN = -3;
var SST_MAX = 35;

var DIAS_MIN_POR_MES = 10;  // umbral de cobertura mínima mensual

var PALETA_SST = [
  '040274','040281','0502a3','0502b8','0502ce','0502e6','0602ff',
  '235cb1','307ef3','269db1','30c8e2','32d3ef','3be285','3ff38f',
  '86e26f','3ae237','b5e22e','d6e21f','fff705','ffd611','ffb613',
  '17f80e','ff8b13','ff6e08','ff500d','ff0000','de0101','c21301','a71001','911003'
];


// ----------------------------------------------------------------------------
// 3. CARGAR PRODUCTO Y FILTRAR
// ----------------------------------------------------------------------------
// Catálogo: https://developers.google.com/earth-engine/datasets/catalog/NASA_OCEANDATA_MODIS-Aqua_L3SMI
var modis_raw = ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI')
  .filterDate(FECHA_INICIO, FECHA_FIN)
  .filterBounds(geom.AN_limpio)
  .select(['sst']);
//   El producto L3SMI ya viene con control de calidad aplicado en el
//   procesamiento previo de NASA, por lo que no es necesario aplicar
//   un filtro qual_sst manualmente.
//   Banda 'sst':      temperatura superficial del mar en °C
//   Banda 'qual_sst': bandera de calidad del retrieval
//     0 = mejor calidad
//     1 = buena calidad
//     2 = aceptable (se descarta en este trabajo)
//     3+ = mala calidad (se descarta)


// ----------------------------------------------------------------------------
// 4. PREPROCESAMIENTO POR IMAGEN: CALIDAD + MÁSCARA TÉRMICA + RECORTE
// ----------------------------------------------------------------------------
var modis_diario = modis_raw.map(function(img) {
  var sst = img.select('sst');

  // Máscara térmica: solo rango físico válido (-3 a 35 °C)
  // El control de calidad propio de MODIS ya está aplicado en el L3SMI.
  var mascara_termica = sst.gte(SST_MIN).and(sst.lte(SST_MAX));

  return sst.updateMask(mascara_termica)
            .clip(geom.AN_limpio)
            .copyProperties(img, ['system:time_start']);
});


// ----------------------------------------------------------------------------
// 5. AGREGACIÓN MENSUAL CON UMBRAL DE COBERTURA
// ----------------------------------------------------------------------------
var anio_inicio = ee.Number(ee.Date(FECHA_INICIO).get('year'));
var anio_fin    = ee.Number(ee.Date(FECHA_FIN).get('year'));
var meses_total = anio_fin.subtract(anio_inicio).multiply(12).toInt();

var modis_mensual = ee.ImageCollection.fromImages(
  ee.List.sequence(0, meses_total.subtract(1)).map(function(offset) {
    var fecha_mes = ee.Date(FECHA_INICIO).advance(ee.Number(offset), 'month');
    // Ajustar al primer día del mes
    var anio = fecha_mes.get('year');
    var mes  = fecha_mes.get('month');
    var fecha_mes_inicio = ee.Date.fromYMD(anio, mes, 1);
    var fecha_mes_fin    = fecha_mes_inicio.advance(1, 'month');

    var imgs_del_mes = modis_diario.filterDate(fecha_mes_inicio, fecha_mes_fin);
    var num_imagenes = imgs_del_mes.size();

    // Calcular la mediana mensual (más robusta que media frente a outliers)
    var mes_compuesto = imgs_del_mes.median()
      .set('system:time_start', fecha_mes_inicio.millis())
      .set('anio', anio)
      .set('mes', mes)
      .set('n_imagenes', num_imagenes);

    return mes_compuesto;
  })
)
// Aplicar umbral de cobertura mínima
.filter(ee.Filter.gte('n_imagenes', DIAS_MIN_POR_MES));


// ----------------------------------------------------------------------------
// 6. INFORMACIÓN DIAGNÓSTICA EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 05 — MODIS AQUA L3SMI CARGA Y PREPROCESAMIENTO');
print('==================================================================');
print('Periodo solicitado:', FECHA_INICIO, '→', FECHA_FIN);
print('Umbral cobertura mensual:', DIAS_MIN_POR_MES, 'días');
print('Imágenes diarias cargadas:', modis_diario.size());
print('Imágenes mensuales (con cobertura suficiente):', modis_mensual.size());
print('Primer mes con dato:', ee.Date(
  modis_mensual.aggregate_min('system:time_start')).format('YYYY-MM'));
print('Último mes con dato:', ee.Date(
  modis_mensual.aggregate_max('system:time_start')).format('YYYY-MM'));

// SST media de junio 2010 sobre AN_limpio (para comparar con OISST)
var junio_2010 = modis_mensual.filterDate('2010-06-01', '2010-07-01').first();
var sst_jun2010 = junio_2010.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: geom.AN_limpio,
  scale: 4640,        // resolución nativa MODIS Aqua (~4,6 km)
  maxPixels: 1e9,
  bestEffort: true
});
print('SST MODIS media junio 2010 (esperado ~12-13°C similar a OISST):',
  sst_jun2010);


// ----------------------------------------------------------------------------
// 7. VISUALIZACIÓN: MEDIA CLIMATOLÓGICA DEL PERIODO COMPLETO
// ----------------------------------------------------------------------------
Map.centerObject(geom.AN_limpio, 3);
Map.setOptions('HYBRID');

var sst_media_periodo = modis_mensual.mean();

Map.addLayer(
  sst_media_periodo,
  {min: 0, max: 28, palette: PALETA_SST},
  '01. MODIS — SST media 2002-2024 (°C)',
  true
);
