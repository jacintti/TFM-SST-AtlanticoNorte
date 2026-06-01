// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 06c: NOAA OISST — TENDENCIAS POR REGIÓN EN PERIODO 2003-2021
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Calcular tendencias OISST en exactamente el mismo periodo que MODIS
//   (2003-2021, 19 años completos). Sirve como BENCHMARK para la
//   validación cruzada con MODIS (script 06a).
//
// COMPARATIVA:
//   Si MODIS y OISST dan tendencias similares en una región sobre el mismo
//   periodo → ambos productos se validan mutuamente para esa región.
//   Si discrepan → hay que investigar la causa (cobertura MODIS, sesgos
//   metodológicos, frentes oceánicos móviles, etc.).
//
// REFERENCIA DEL PRODUCTO:
//   Huang, B., Liu, C., Banzon, V., et al. (2021). Improvements of the Daily
//   Optimum Interpolation Sea Surface Temperature (DOISST) Version 2.1.
//   Journal of Climate, 34(8), 2923-2939.
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

// Periodo de análisis (mismo que MODIS 06a)
var TEND_ANIO_INICIO = 2003;
var TEND_ANIO_FIN    = 2021;


// ----------------------------------------------------------------------------
// 3. CARGAR Y PREPROCESAR OISST
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
    var imgs = oisst_diario.filterDate(fecha_mes, fecha_mes.advance(1, 'month'));
    var n = imgs.size();
    return imgs.median()
      .set('system:time_start', fecha_mes.millis())
      .set('anio', fecha_mes.get('year'))
      .set('mes', fecha_mes.get('month'))
      .set('n_imagenes', n);
  })
).filter(ee.Filter.gt('n_imagenes', 0));


// ----------------------------------------------------------------------------
// 4. AGREGACIÓN A SERIE ANUAL 2003-2021
// ----------------------------------------------------------------------------
var anios = ee.List.sequence(TEND_ANIO_INICIO, TEND_ANIO_FIN);

var oisst_anual = ee.ImageCollection.fromImages(
  anios.map(function(a) {
    var meses = oisst_mensual.filter(ee.Filter.eq('anio', a));
    return meses.mean()
      .set('anio', a)
      .set('system:time_start', ee.Date.fromYMD(a, 7, 1).millis());
  })
);


// ----------------------------------------------------------------------------
// 5. FUNCIÓN DE TENDENCIA POR REGIÓN
// ----------------------------------------------------------------------------
var calcular_tendencia = function(nombre, geometria) {
  var serie = oisst_anual.map(function(img) {
    var v = img.select('sst').reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geometria,
      scale: 25000,
      maxPixels: 1e9,
      bestEffort: true
    });
    return ee.Feature(null, {
      'anio': img.get('anio'),
      'sst':  v.get('sst')
    });
  }).filter(ee.Filter.notNull(['sst']))
    .sort('anio');

  var n = serie.size();

  var sen = serie.reduceColumns({
    reducer: ee.Reducer.sensSlope(),
    selectors: ['anio', 'sst']
  });
  var mk = serie.reduceColumns({
    reducer: ee.Reducer.kendallsCorrelation(2),
    selectors: ['anio', 'sst']
  });

  var tau = ee.Number(mk.get('tau'));
  var n_num = ee.Number(n);
  var z_factor = n_num.multiply(9).multiply(n_num.subtract(1))
    .divide(ee.Number(2).multiply(n_num.multiply(2).add(5)));
  var z = tau.multiply(z_factor.sqrt());
  var z_abs = z.abs();
  var t = ee.Number(1).divide(ee.Number(1).add(z_abs.multiply(0.2316419)));
  var d = ee.Number(0.3989423)
    .multiply(z_abs.multiply(-1).multiply(z_abs).divide(2).exp());
  var phi = ee.Number(1).subtract(d.multiply(
    t.multiply(0.3193815)
      .add(t.pow(2).multiply(-0.3565638))
      .add(t.pow(3).multiply(1.781478))
      .add(t.pow(4).multiply(-1.821256))
      .add(t.pow(5).multiply(1.330274))
  ));
  var p_valor = ee.Number(2).multiply(ee.Number(1).subtract(phi));

  return ee.Feature(null, {
    'producto':               'OISST',
    'region':                 nombre,
    'n_anios':                n,
    'tendencia_C_por_decada': ee.Number(sen.get('slope')).multiply(10),
    'tau_kendall':            tau,
    'p_valor':                p_valor,
    'significativo_p005':     p_valor.lt(0.05)
  });
};


// ----------------------------------------------------------------------------
// 6. APLICAR A LAS SIETE REGIONES
// ----------------------------------------------------------------------------
var regiones_lista = [
  {nombre: 'AN_limpio', geometria: geom.AN_limpio},
  {nombre: 'GSE',       geometria: geom.GSE},
  {nombre: 'LAB',       geometria: geom.LAB},
  {nombre: 'NOR',       geometria: geom.NOR},
  {nombre: 'EUR',       geometria: geom.EUR},
  {nombre: 'CBL',       geometria: geom.CBL},
  {nombre: 'T53',       geometria: geom.T53}
];

var resultados = ee.FeatureCollection(regiones_lista.map(function(r) {
  return calcular_tendencia(r.nombre, r.geometria);
}));


// ----------------------------------------------------------------------------
// 7. INFORMACIÓN EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 06c — TENDENCIAS OISST POR REGIÓN (2003-2021)');
print('==================================================================');
print('Periodo:', TEND_ANIO_INICIO, 'a', TEND_ANIO_FIN, '(19 años)');
print('Resolución de reducción regional:', 25000, 'm');
print('Este resultado es el BENCHMARK para validación cruzada con MODIS 06a.');
print('');
print('TABLA DE TENDENCIAS OISST 2003-2021:');
print(resultados);