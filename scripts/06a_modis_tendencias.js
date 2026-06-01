// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 06a: MODIS AQUA — TENDENCIAS POR REGIÓN (2003-2021)
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Calcular tendencias Mann-Kendall + Sen's slope para las siete regiones
//   del área de estudio usando MODIS Aqua sobre el periodo 2003-2021
//   (19 años completos).
//
// SALIDA:
//   · Tabla en consola con las tendencias por región (°C/década, tau, p-valor).
//
// VALIDACIÓN CRUZADA:
//   Esta tabla se comparará con la equivalente OISST 2003-2021 que produce
//   el script 06b. Si los números son consistentes, OISST queda validado.
//
// REFERENCIA DEL PRODUCTO:
//   Kilpatrick, K.A., Podestá, G., Walsh, S., et al. (2015). A decade of sea
//   surface temperature from MODIS. Remote Sensing of Environment, 165, 27-41.
// ============================================================================


// ----------------------------------------------------------------------------
// 1. IMPORTAR MÓDULO DE GEOMETRÍAS
// ----------------------------------------------------------------------------
var geom = require('users/jacintod/TFM1:00_geometrias_module');


// ----------------------------------------------------------------------------
// 2. PARÁMETROS
// ----------------------------------------------------------------------------
var FECHA_INICIO = '2002-07-04';
var FECHA_FIN    = '2022-03-01';
var SST_MIN = -3;
var SST_MAX = 35;
var DIAS_MIN_POR_MES = 10;

var TEND_ANIO_INICIO = 2003;
var TEND_ANIO_FIN    = 2021;


// ----------------------------------------------------------------------------
// 3. CARGAR Y PREPROCESAR MODIS (réplica del script 05)
// ----------------------------------------------------------------------------
var modis_raw = ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI')
  .filterDate(FECHA_INICIO, FECHA_FIN)
  .filterBounds(geom.AN_limpio)
  .select(['sst']);

var modis_diario = modis_raw.map(function(img) {
  var sst = img.select('sst');
  var mascara = sst.gte(SST_MIN).and(sst.lte(SST_MAX));
  return sst.updateMask(mascara)
            .clip(geom.AN_limpio)
            .copyProperties(img, ['system:time_start']);
});

var anio_inicio = ee.Number(ee.Date(FECHA_INICIO).get('year'));
var anio_fin    = ee.Number(ee.Date(FECHA_FIN).get('year'));
var meses_total = anio_fin.subtract(anio_inicio).multiply(12).toInt();

var modis_mensual = ee.ImageCollection.fromImages(
  ee.List.sequence(0, meses_total.subtract(1)).map(function(offset) {
    var fecha_mes = ee.Date(FECHA_INICIO).advance(ee.Number(offset), 'month');
    var anio = fecha_mes.get('year');
    var mes  = fecha_mes.get('month');
    var fi = ee.Date.fromYMD(anio, mes, 1);
    var ff = fi.advance(1, 'month');
    var imgs = modis_diario.filterDate(fi, ff);
    var n = imgs.size();
    return imgs.median()
      .set('system:time_start', fi.millis())
      .set('anio', anio)
      .set('mes', mes)
      .set('n_imagenes', n);
  })
).filter(ee.Filter.gte('n_imagenes', DIAS_MIN_POR_MES));


// ----------------------------------------------------------------------------
// 4. AGREGACIÓN A SERIE ANUAL 2003-2021
// ----------------------------------------------------------------------------
var anios = ee.List.sequence(TEND_ANIO_INICIO, TEND_ANIO_FIN);

var modis_anual = ee.ImageCollection.fromImages(
  anios.map(function(a) {
    var meses = modis_mensual.filter(ee.Filter.eq('anio', a));
    return meses.mean()
      .set('anio', a)
      .set('system:time_start', ee.Date.fromYMD(a, 7, 1).millis());
  })
);


// ----------------------------------------------------------------------------
// 5. FUNCIÓN DE TENDENCIA POR REGIÓN
// ----------------------------------------------------------------------------
var calcular_tendencia = function(nombre, geometria) {
  var serie = modis_anual.map(function(img) {
    var v = img.select('sst').reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: geometria,
  scale: 9280,
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
    'producto':               'MODIS',
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
print('SCRIPT 06a — TENDENCIAS MODIS POR REGIÓN');
print('==================================================================');
print('Periodo:', TEND_ANIO_INICIO, 'a', TEND_ANIO_FIN, '(19 años)');
print('Resolución de reducción regional:', 9280, 'm');
print('');
print('TABLA DE TENDENCIAS MODIS:');
print(resultados);
