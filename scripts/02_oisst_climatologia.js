// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 02: NOAA OISST v2.1 — CLIMATOLOGÍAS Y ANOMALÍAS
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   A partir de la colección mensual preprocesada (replicada aquí desde el
//   script 01), calcular:
//     · La climatología 1991-2020: 12 imágenes, una por mes calendario.
//     · La colección de anomalías mensuales: serie completa menos climatología.
//
// SALIDA (en memoria):
//   · climatologia: 12 imágenes mensuales (referencia 1991-2020).
//   · oisst_anomalia: colección mensual de anomalías.
//
// VERIFICACIÓN ESPERADA:
//   · climatologia.size() = 12.
//   · oisst_anomalia.size() ≈ 520.
//   · Anomalía media sobre AN_limpio en 1991-2020 ≈ 0 (por construcción).
//   · Anomalía media sobre AN_limpio en 2010-2019 ≈ +0.3 a +0.5 °C (positiva).
//
// VISUALIZACIÓN:
//   · Climatología enero (mes frío de referencia).
//   · Climatología julio (mes cálido de referencia).
//   · Anomalía media 2010-2019 — mapa clave del trabajo.
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

var REF_INICIO = '1991-01-01';    // climatología base 1991-2020 (OMM)
var REF_FIN    = '2021-01-01';

var DECADA_INICIO = '2010-01-01';  // década reciente para verificar anomalía
var DECADA_FIN    = '2020-01-01';

var SST_MIN = -3;
var SST_MAX = 35;

var PALETA_SST = [
  '040274','040281','0502a3','0502b8','0502ce','0502e6','0602ff',
  '235cb1','307ef3','269db1','30c8e2','32d3ef','3be285','3ff38f',
  '86e26f','3ae237','b5e22e','d6e21f','fff705','ffd611','ffb613',
  '17f80e','ff8b13','ff6e08','ff500d','ff0000','de0101','c21301','a71001','911003'
];

// Paleta divergente azul-blanco-rojo para anomalías
var PALETA_ANOMALIA = [
  '053061','2166ac','4393c3','92c5de','d1e5f0',
  'ffffff',
  'fddbc7','f4a582','d6604d','b2182b','67001f'
];


// ----------------------------------------------------------------------------
// 3. RECARGAR Y PREPROCESAR COLECCIÓN OISST (replica del script 01)
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

// Agregación mensual (mediana de los días del mes)
var anio_inicio = ee.Number(ee.Date(FECHA_INICIO).get('year'));
var anio_fin    = ee.Number(ee.Date(FECHA_FIN).get('year'));
var meses_total = anio_fin.subtract(anio_inicio).multiply(12).toInt();

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


// ----------------------------------------------------------------------------
// 4. CLIMATOLOGÍA 1991-2020: UNA IMAGEN POR MES CALENDARIO
// ----------------------------------------------------------------------------
// Filtrar la serie al periodo de referencia
var oisst_referencia = oisst_mensual.filterDate(REF_INICIO, REF_FIN);

// Lista de meses calendario (1 a 12)
var meses_calendario = ee.List.sequence(1, 12);

// Para cada mes calendario, calcular la media sobre los 30 años de referencia
var climatologia = ee.ImageCollection.fromImages(
  meses_calendario.map(function(m) {
    var imgs_este_mes = oisst_referencia.filter(ee.Filter.eq('mes', m));
    return imgs_este_mes.mean()
      .set('mes_calendario', m)
      .rename('sst_climatologia');
  })
);


// ----------------------------------------------------------------------------
// 5. ANOMALÍAS MENSUALES: SERIE COMPLETA − CLIMATOLOGÍA DEL MISMO MES
// ----------------------------------------------------------------------------
var oisst_anomalia = oisst_mensual.map(function(img) {
  var mes_de_img = ee.Number(img.get('mes'));
  // Encontrar la imagen climatológica del mes correspondiente
  var clim_mes = climatologia
    .filter(ee.Filter.eq('mes_calendario', mes_de_img))
    .first();
  // Renombramos a 'sst' la banda climatológica antes de restar
  // para mantener consistencia de nombres de banda
  return img.subtract(clim_mes.rename('sst'))
            .rename('anomalia')
            .copyProperties(img, ['system:time_start', 'anio', 'mes']);
});


// ----------------------------------------------------------------------------
// 6. INFORMACIÓN DIAGNÓSTICA EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 02 — CLIMATOLOGÍAS Y ANOMALÍAS OISST');
print('==================================================================');
print('Periodo de referencia climatológico:', REF_INICIO, '→', REF_FIN);
print('Climatología — número de imágenes (esperado 12):', climatologia.size());
print('Anomalías — número de imágenes:', oisst_anomalia.size());

// Verificación 1: anomalía media en periodo de referencia ≈ 0
var anom_referencia = oisst_anomalia
  .filterDate(REF_INICIO, REF_FIN)
  .mean()
  .reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geom.AN_limpio,
    scale: 25000,
    maxPixels: 1e9,
    bestEffort: true
  });
print('Anomalía media 1991-2020 sobre AN_limpio (esperado ~0):',
  anom_referencia);

// Verificación 2: anomalía media década reciente
var anom_decada = oisst_anomalia
  .filterDate(DECADA_INICIO, DECADA_FIN)
  .mean()
  .reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geom.AN_limpio,
    scale: 25000,
    maxPixels: 1e9,
    bestEffort: true
  });
print('Anomalía media 2010-2019 sobre AN_limpio (esperado +0.3 a +0.5 °C):',
  anom_decada);


// ----------------------------------------------------------------------------
// 7. VISUALIZACIÓN
// ----------------------------------------------------------------------------
Map.centerObject(geom.AN_limpio, 3);
Map.setOptions('HYBRID');

// Climatología enero
var clim_enero = climatologia.filter(ee.Filter.eq('mes_calendario', 1)).first();
Map.addLayer(
  clim_enero,
  {min: 0, max: 28, palette: PALETA_SST},
  '01. Climatología ENERO 1991-2020 (°C)',
  false
);

// Climatología julio
var clim_julio = climatologia.filter(ee.Filter.eq('mes_calendario', 7)).first();
Map.addLayer(
  clim_julio,
  {min: 0, max: 28, palette: PALETA_SST},
  '02. Climatología JULIO 1991-2020 (°C)',
  false
);

// Anomalía media 2010-2019 — el mapa clave
var anomalia_2010s = oisst_anomalia
  .filterDate(DECADA_INICIO, DECADA_FIN)
  .mean();
Map.addLayer(
  anomalia_2010s,
  {min: -1, max: 1, palette: PALETA_ANOMALIA},
  '03. Anomalía media 2010-2019 (°C)',
  true   // capa visible por defecto
);


// ----------------------------------------------------------------------------
// 8. NOTA OPERATIVA
// ----------------------------------------------------------------------------
// El recálculo desde cero (scripts 01 → 02 → 03 ...) hace que cada ejecución
// dure 30-90 segundos. En el siguiente script 03 (reducción regional) podemos
// importar este módulo con require() para evitar el recálculo.

