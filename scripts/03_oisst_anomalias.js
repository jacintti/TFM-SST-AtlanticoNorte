// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 03: NOAA OISST v2.1 — REDUCCIÓN REGIONAL POR GEOMETRÍAS
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   A partir de las colecciones mensual y de anomalías OISST (replicadas
//   aquí desde el script 02), reducir cada imagen a una cifra por geometría
//   mediante reduceRegion(media). Generar series temporales univariantes
//   por geometría y visualizarlas.
//
// SALIDA (en memoria):
//   · serie_sst[geom]: serie temporal mensual de SST absoluta (°C)
//   · serie_anomalia[geom]: serie temporal mensual de anomalía (°C)
//
// VERIFICACIÓN ESPERADA:
//   · GSE jun 2010: ~22-24 °C (Gulf Stream Extension, verano)
//   · NOR ene 2015: ~2-5 °C (Mares Nórdicos, invierno)
//   · LAB y giro subpolar: tendencia a anomalías negativas en 2010s
//
// GRÁFICOS:
//   · Evolución anual SST 1982-2024 para GSE, LAB, NOR, EUR.
//   · Heatmap de anomalías anuales por polígono (preludio de G.7).
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
var REF_INICIO   = '1991-01-01';
var REF_FIN      = '2021-01-01';

var SST_MIN = -3;
var SST_MAX = 35;
var ESCALA_OISST = 25000;   // 25 km, resolución nativa OISST

// Colores para los gráficos
var COLORES = {
  'AN_limpio': '#333333',
  'GSE':       '#d62728',   // rojo - Gulf Stream
  'LAB':       '#1f77b4',   // azul - Labrador
  'NOR':       '#9467bd',   // morado - Mares Nórdicos
  'EUR':       '#ff7f0e',   // naranja - Costa europea
  'T53':       '#2ca02c'    // verde - transecto
};


// ----------------------------------------------------------------------------
// 3. RECONSTRUIR COLECCIONES OISST (replica de scripts 01 y 02)
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

// Climatología 1991-2020
var oisst_referencia = oisst_mensual.filterDate(REF_INICIO, REF_FIN);
var climatologia = ee.ImageCollection.fromImages(
  ee.List.sequence(1, 12).map(function(m) {
    var imgs_este_mes = oisst_referencia.filter(ee.Filter.eq('mes', m));
    return imgs_este_mes.mean().set('mes_calendario', m).rename('sst_climatologia');
  })
);

// Anomalías
var oisst_anomalia = oisst_mensual.map(function(img) {
  var mes_de_img = ee.Number(img.get('mes'));
  var clim_mes = climatologia
    .filter(ee.Filter.eq('mes_calendario', mes_de_img))
    .first();
  return img.subtract(clim_mes.rename('sst'))
            .rename('anomalia')
            .copyProperties(img, ['system:time_start', 'anio', 'mes']);
});


// ----------------------------------------------------------------------------
// 4. FUNCIÓN PARA REDUCIR UNA COLECCIÓN SOBRE UNA GEOMETRÍA
// ----------------------------------------------------------------------------
// Toma una colección y una geometría, devuelve una FeatureCollection de
// Features sin geometría con propiedades: fecha, anio, mes, valor.

var reducir_a_serie = function(coleccion, geometria, banda) {
  return coleccion.map(function(img) {
    var reduccion = img.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geometria,
      scale: ESCALA_OISST,
      maxPixels: 1e9,
      bestEffort: true
    });
    return ee.Feature(null, {
      'fecha': img.date().format('YYYY-MM'),
      'anio':  img.get('anio'),
      'mes':   img.get('mes'),
      'valor': reduccion.get(banda)
    });
  })
  .filter(ee.Filter.notNull(['valor']));
};


// ----------------------------------------------------------------------------
// 5. CALCULAR SERIES DE SST ABSOLUTA PARA LAS GEOMETRÍAS REGIONALES + T53
// ----------------------------------------------------------------------------
var serie_sst_AN  = reducir_a_serie(oisst_mensual, geom.AN_limpio, 'sst');
var serie_sst_GSE = reducir_a_serie(oisst_mensual, geom.GSE,       'sst');
var serie_sst_LAB = reducir_a_serie(oisst_mensual, geom.LAB,       'sst');
var serie_sst_NOR = reducir_a_serie(oisst_mensual, geom.NOR,       'sst');
var serie_sst_EUR = reducir_a_serie(oisst_mensual, geom.EUR,       'sst');
var serie_sst_T53 = reducir_a_serie(oisst_mensual, geom.T53,       'sst');


// ----------------------------------------------------------------------------
// 6. CALCULAR SERIES DE ANOMALÍA PARA LAS MISMAS GEOMETRÍAS
// ----------------------------------------------------------------------------
var serie_anom_AN  = reducir_a_serie(oisst_anomalia, geom.AN_limpio, 'anomalia');
var serie_anom_GSE = reducir_a_serie(oisst_anomalia, geom.GSE,       'anomalia');
var serie_anom_LAB = reducir_a_serie(oisst_anomalia, geom.LAB,       'anomalia');
var serie_anom_NOR = reducir_a_serie(oisst_anomalia, geom.NOR,       'anomalia');
var serie_anom_EUR = reducir_a_serie(oisst_anomalia, geom.EUR,       'anomalia');
var serie_anom_T53 = reducir_a_serie(oisst_anomalia, geom.T53,       'anomalia');


// ----------------------------------------------------------------------------
// 7. SERIES DE LAS VENTANAS COSTERAS (OISST sobre franja oceánica costera)
// ----------------------------------------------------------------------------
var serie_sst_L40w = reducir_a_serie(oisst_mensual, geom.L40w, 'sst');
var serie_sst_L40e = reducir_a_serie(oisst_mensual, geom.L40e, 'sst');
var serie_sst_L53w = reducir_a_serie(oisst_mensual, geom.L53w, 'sst');
var serie_sst_L53e = reducir_a_serie(oisst_mensual, geom.L53e, 'sst');
var serie_sst_L67w = reducir_a_serie(oisst_mensual, geom.L67w, 'sst');
var serie_sst_L67e = reducir_a_serie(oisst_mensual, geom.L67e, 'sst');


// ----------------------------------------------------------------------------
// 8. INFORMACIÓN DIAGNÓSTICA EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 03 — REDUCCIÓN REGIONAL OISST');
print('==================================================================');
print('Número de meses en serie SST AN_limpio:', serie_sst_AN.size());
print('Número de meses en serie SST GSE:',        serie_sst_GSE.size());
print('Número de meses en serie SST LAB:',        serie_sst_LAB.size());

// Verificación 1: SST GSE junio 2010
var gse_jun2010 = serie_sst_GSE
  .filter(ee.Filter.and(ee.Filter.eq('anio', 2010), ee.Filter.eq('mes', 6)))
  .first();
print('SST GSE junio 2010 (esperado ~22-24 °C):',
  gse_jun2010.get('valor'));

// Verificación 2: SST NOR enero 2015
var nor_ene2015 = serie_sst_NOR
  .filter(ee.Filter.and(ee.Filter.eq('anio', 2015), ee.Filter.eq('mes', 1)))
  .first();
print('SST NOR enero 2015 (esperado ~2-5 °C):',
  nor_ene2015.get('valor'));


// ----------------------------------------------------------------------------
// 9. GRÁFICO: EVOLUCIÓN SST ANUAL DE LOS POLÍGONOS REGIONALES
// ----------------------------------------------------------------------------
// Para graficar, conviene calcular la media anual de cada serie.
// Construimos una función auxiliar que agrupa por año.

var agregar_anual = function(serie_mensual) {
  var anios = ee.List.sequence(1982, 2024);
  return ee.FeatureCollection(
    anios.map(function(a) {
      var del_anio = serie_mensual.filter(ee.Filter.eq('anio', a));
      var media = del_anio.aggregate_mean('valor');
      return ee.Feature(null, {'anio': a, 'valor': media});
    })
  ).filter(ee.Filter.notNull(['valor']));
};

var anual_GSE = agregar_anual(serie_sst_GSE);
var anual_LAB = agregar_anual(serie_sst_LAB);
var anual_NOR = agregar_anual(serie_sst_NOR);
var anual_EUR = agregar_anual(serie_sst_EUR);

// Para graficar varias series en el mismo plot, las combinamos en una sola
// FeatureCollection con una propiedad por polígono.
var combinado_anual = anual_GSE.map(function(f) {
  var anio = f.get('anio');
  var gse = f.get('valor');
  var lab = anual_LAB.filter(ee.Filter.eq('anio', anio)).first().get('valor');
  var nor = anual_NOR.filter(ee.Filter.eq('anio', anio)).first().get('valor');
  var eur = anual_EUR.filter(ee.Filter.eq('anio', anio)).first().get('valor');
  return ee.Feature(null, {
    'anio': anio, 'GSE': gse, 'LAB': lab, 'NOR': nor, 'EUR': eur
  });
});

var grafico_evolucion = ui.Chart.feature.byFeature(
    combinado_anual, 'anio', ['GSE', 'LAB', 'NOR', 'EUR'])
  .setChartType('LineChart')
  .setOptions({
    title: 'OISST — Evolución anual SST 1982-2024 (regional)',
    hAxis: {title: 'Año', format: '####'},
    vAxis: {title: 'SST media anual (°C)'},
    lineWidth: 2,
    pointSize: 3,
    colors: [COLORES.GSE, COLORES.LAB, COLORES.NOR, COLORES.EUR],
    legend: {position: 'top'}
  });
print(grafico_evolucion);


// ----------------------------------------------------------------------------
// 10. GRÁFICO: EVOLUCIÓN ANOMALÍA ANUAL POR POLÍGONO (preludio G.7)
// ----------------------------------------------------------------------------
var anomanual_GSE = agregar_anual(serie_anom_GSE);
var anomanual_LAB = agregar_anual(serie_anom_LAB);
var anomanual_NOR = agregar_anual(serie_anom_NOR);
var anomanual_EUR = agregar_anual(serie_anom_EUR);
var anomanual_AN  = agregar_anual(serie_anom_AN);

var combinado_anom = anomanual_AN.map(function(f) {
  var anio = f.get('anio');
  var an  = f.get('valor');
  var gse = anomanual_GSE.filter(ee.Filter.eq('anio', anio)).first().get('valor');
  var lab = anomanual_LAB.filter(ee.Filter.eq('anio', anio)).first().get('valor');
  var nor = anomanual_NOR.filter(ee.Filter.eq('anio', anio)).first().get('valor');
  var eur = anomanual_EUR.filter(ee.Filter.eq('anio', anio)).first().get('valor');
  return ee.Feature(null, {
    'anio': anio, 'AN_limpio': an, 'GSE': gse, 'LAB': lab, 'NOR': nor, 'EUR': eur
  });
});

var grafico_anomalia = ui.Chart.feature.byFeature(
    combinado_anom, 'anio', ['AN_limpio', 'GSE', 'LAB', 'NOR', 'EUR'])
  .setChartType('LineChart')
  .setOptions({
    title: 'OISST — Anomalía anual respecto a 1991-2020 (regional)',
    hAxis: {title: 'Año', format: '####'},
    vAxis: {title: 'Anomalía (°C)'},
    lineWidth: 2,
    pointSize: 3,
    colors: [COLORES.AN_limpio, COLORES.GSE, COLORES.LAB, COLORES.NOR, COLORES.EUR],
    legend: {position: 'top'}
  });
print(grafico_anomalia);


// ----------------------------------------------------------------------------
// 11. NOTA OPERATIVA
// ----------------------------------------------------------------------------
// Este script construye en memoria 12 FeatureCollection (6 SST + 6 anomalía
// para los 6 polígonos regionales + transecto + 6 ventanas costeras).
