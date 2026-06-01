// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 06b: MODIS AQUA — DIAGNÓSTICO DE COBERTURA POR REGIÓN
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Contar, para cada combinación de región y año (2003-2021), cuántos
//   meses tienen un dato válido tras aplicar nuestro filtro de cobertura
//   mensual (mínimo 10 días con dato por mes).
//
//   Sirve para diagnosticar dónde MODIS es fiable y dónde no. Las regiones
//   en latitudes altas (LAB, NOR) suelen tener huecos importantes por
//   nubosidad e hielo marino. Las regiones más al sur (GSE, EUR) deberían
//   tener cobertura prácticamente completa.
//
// SALIDA:
//   · Tabla en consola con número de meses con dato por región y año.
//   · Resumen por región: media y mínimo de meses por año.
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
// 3. CARGAR Y PREPROCESAR MODIS (réplica del 06a)
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
// 4. FUNCIÓN PARA CONTAR MESES CON DATO POR REGIÓN Y AÑO
// ----------------------------------------------------------------------------
// Para cada combinación (región, año), contamos cuántos de los 12 meses
// del año dan un valor no nulo al reducir sobre la región.

var contar_meses_validos = function(nombre, geometria) {
  var anios = ee.List.sequence(TEND_ANIO_INICIO, TEND_ANIO_FIN);

  return anios.map(function(a) {
    var meses_del_anio = modis_mensual.filter(ee.Filter.eq('anio', a));

    // Para cada mes del año, comprobamos si la reducción sobre la región
    // devuelve un valor no nulo.
    var cuenta = meses_del_anio.map(function(img) {
      var v = img.select('sst').reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: geometria,
        scale: 25000,
        maxPixels: 1e9,
        bestEffort: true
      });
      var valor = v.get('sst');
      return ee.Feature(null, {
        'tiene_dato': ee.Algorithms.If(valor, 1, 0)
      });
    }).aggregate_sum('tiene_dato');

    return ee.Feature(null, {
      'region':  nombre,
      'anio':    a,
      'n_meses_validos': cuenta
    });
  });
};


// ----------------------------------------------------------------------------
// 5. APLICAR A LAS SIETE REGIONES
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

// Construir una FeatureCollection con todos los pares (región, año)
var todos_resultados = ee.FeatureCollection(
  regiones_lista.reduce(function(acumulado, r) {
    return acumulado.concat(contar_meses_validos(r.nombre, r.geometria));
  }, [])
);

print('==================================================================');
print('SCRIPT 06b — DIAGNÓSTICO DE COBERTURA MODIS POR REGIÓN Y AÑO');
print('==================================================================');
print('Para cada región y año, n_meses_validos indica cuántos meses');
print('del año tienen un valor calculable en la reducción regional.');
print('Valor máximo posible: 12 (todos los meses con dato).');
print('');
print('Tabla detallada:');
print(todos_resultados);


// ----------------------------------------------------------------------------
// 6. RESUMEN POR REGIÓN
// ----------------------------------------------------------------------------
// Para cada región, calculamos media y mínimo de meses por año.

var resumen_por_region = regiones_lista.map(function(r) {
  var filtrada = todos_resultados.filter(ee.Filter.eq('region', r.nombre));
  var media = filtrada.aggregate_mean('n_meses_validos');
  var minimo = filtrada.aggregate_min('n_meses_validos');
  return ee.Feature(null, {
    'region': r.nombre,
    'media_meses_por_anio': media,
    'minimo_meses_por_anio': minimo
  });
});

print('');
print('Resumen por región:');
print(ee.FeatureCollection(resumen_por_region));      
