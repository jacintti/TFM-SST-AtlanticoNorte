// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 04b: NOAA OISST v2.1 — TENDENCIAS POR SUBPERIODOS
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Complementar el análisis de tendencias del script 04 (que cubre el
//   periodo completo 1982-2024) calculando tendencias separadas para dos
//   subperiodos:
//     · 1982-2000 (~19 años): fase histórica previa al cold blob.
//     · 2001-2024 (~24 años): fase reciente donde el cold blob aparece.
//
//   Objetivo científico: detectar la firma temporal del cold blob, en
//   especial en la región CBL, comparando si su tendencia ha cambiado
//   entre ambos subperiodos.
//
// REFERENCIAS: Rahmstorf et al. (2015), Caesar et al. (2018) documentan
//   el cold blob como fenómeno especialmente acentuado desde ~2000.
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
var ANIO_CORTE = 2000;  // 1982-2000 (sub 1) y 2001-2024 (sub 2)


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


// ----------------------------------------------------------------------------
// 4. AGREGACIÓN A SERIE ANUAL
// ----------------------------------------------------------------------------
var anios_serie = ee.List.sequence(ANIO_INICIO_SERIE, ANIO_FIN_SERIE);

var oisst_anual = ee.ImageCollection.fromImages(
  anios_serie.map(function(a) {
    var meses_del_anio = oisst_mensual.filter(ee.Filter.eq('anio', a));
    return meses_del_anio.mean()
      .set('anio', a)
      .set('system:time_start', ee.Date.fromYMD(a, 7, 1).millis());
  })
);


// ----------------------------------------------------------------------------
// 5. FUNCIÓN DE CÁLCULO DE TENDENCIA POR SUBPERIODO
// ----------------------------------------------------------------------------
var calcular_tendencia_subperiodo = function(nombre, geometria,
                                              anio_min, anio_max) {
  var serie_sub = oisst_anual
    .filter(ee.Filter.gte('anio', anio_min))
    .filter(ee.Filter.lte('anio', anio_max));

  var serie_features = serie_sub.map(function(img) {
    var reduccion = img.select('sst').reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geometria,
      scale: 25000,
      maxPixels: 1e9,
      bestEffort: true
    });
    return ee.Feature(null, {
      'anio': img.get('anio'),
      'sst':  reduccion.get('sst')
    });
  }).filter(ee.Filter.notNull(['sst']))
    .sort('anio');

  var n = serie_features.size();

  var sen = serie_features.reduceColumns({
    reducer: ee.Reducer.sensSlope(),
    selectors: ['anio', 'sst']
  });

  var mk = serie_features.reduceColumns({
    reducer: ee.Reducer.kendallsCorrelation(2),
    selectors: ['anio', 'sst']
  });

  var tau = ee.Number(mk.get('tau'));

  // Cálculo manual de Z y p-valor (aprox. Abramowitz-Stegun)
  var n_num = ee.Number(n);
  var z_denominador = ee.Number(2).multiply(n_num.multiply(2).add(5));
  var z_factor = n_num.multiply(9).multiply(n_num.subtract(1)).divide(z_denominador);
  var z = tau.multiply(z_factor.sqrt());

  var z_abs = z.abs();
  var t = ee.Number(1).divide(ee.Number(1).add(z_abs.multiply(0.2316419)));
  var d = ee.Number(0.3989423).multiply(z_abs.multiply(-1).multiply(z_abs).divide(2).exp());
  var phi = ee.Number(1).subtract(d.multiply(
    t.multiply(0.3193815)
      .add(t.pow(2).multiply(-0.3565638))
      .add(t.pow(3).multiply(1.781478))
      .add(t.pow(4).multiply(-1.821256))
      .add(t.pow(5).multiply(1.330274))
  ));
  var p_valor = ee.Number(2).multiply(ee.Number(1).subtract(phi));

  return ee.Feature(null, {
    'region':                 nombre,
    'subperiodo':             ee.String(ee.Number(anio_min).int())
                                 .cat('-')
                                 .cat(ee.String(ee.Number(anio_max).int())),
    'n_anios':                n,
    'tendencia_C_por_decada': ee.Number(sen.get('slope')).multiply(10),
    'tau_kendall':            tau,
    'p_valor':                p_valor,
    'significativo_p005':     p_valor.lt(0.05)
  });
};

var regiones_lista = [
  {nombre: 'AN_limpio', geometria: geom.AN_limpio},
  {nombre: 'GSE',       geometria: geom.GSE},
  {nombre: 'LAB',       geometria: geom.LAB},
  {nombre: 'NOR',       geometria: geom.NOR},
  {nombre: 'EUR',       geometria: geom.EUR},
  {nombre: 'CBL',       geometria: geom.CBL},
  {nombre: 'T53',       geometria: geom.T53}
];


// ----------------------------------------------------------------------------
// 6. CÁLCULO Y RESULTADOS DEL SUBPERIODO 1 (1982-2000)
// ----------------------------------------------------------------------------
var resultados_sub1 = ee.FeatureCollection(regiones_lista.map(function(r) {
  return calcular_tendencia_subperiodo(r.nombre, r.geometria,
                                        ANIO_INICIO_SERIE, ANIO_CORTE);
}));

print('==================================================================');
print('SCRIPT 04b — TENDENCIAS POR SUBPERIODOS');
print('==================================================================');
print('SUBPERIODO 1: 1982-' + ANIO_CORTE);
print(resultados_sub1);


// ----------------------------------------------------------------------------
// 7. CÁLCULO Y RESULTADOS DEL SUBPERIODO 2 (2001-2024)
// ----------------------------------------------------------------------------
var resultados_sub2 = ee.FeatureCollection(regiones_lista.map(function(r) {
  return calcular_tendencia_subperiodo(r.nombre, r.geometria,
                                        ANIO_CORTE + 1, ANIO_FIN_SERIE);
}));

print('SUBPERIODO 2: ' + (ANIO_CORTE + 1) + '-' + ANIO_FIN_SERIE);
print(resultados_sub2);
