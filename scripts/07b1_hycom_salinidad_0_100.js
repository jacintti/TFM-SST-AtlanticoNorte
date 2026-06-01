// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 07b1: HYCOM — TENDENCIAS DE SALINIDAD EN SUPERFICIE Y 100 m
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Calcular tendencias Mann-Kendall + Sen's slope para la salinidad
//   en los niveles 0 m y 100 m sobre las siete regiones del trabajo,
//   durante 1993-2023 (31 años completos de HYCOM).
//
//   Objetivo científico: detectar la huella de la Gran Anomalía de
//   Salinidad reciente, especialmente en LAB y CBL, donde la entrada de
//   agua dulce ártica ha reducido la salinidad superficial en las últimas
//   décadas (Biló et al. 2022). Esta huella explica físicamente la
//   discrepancia MODIS-OISST observada en LAB en el periodo 2003-2021.
//
// SALIDA:
//   · Tabla con 14 filas (7 regiones × 2 niveles) con tendencia, tau y
//     p-valor de la salinidad.
//
// REFERENCIA:
//   Bil\u00f3, T.C., Straneo, F., Holte, J., & Le Bras, I.A.-A. (2022).
//   Arrival of new Great Salinity Anomaly weakens convection in the
//   Irminger Sea. Geophysical Research Letters, 49(11), e2022GL098857.
// ============================================================================


// ----------------------------------------------------------------------------
// 1. IMPORTAR MÓDULO DE GEOMETRÍAS
// ----------------------------------------------------------------------------
var geom = require('users/jacintod/TFM1:00_geometrias_module');


// ----------------------------------------------------------------------------
// 2. PARÁMETROS
// ----------------------------------------------------------------------------
var FECHA_INICIO = '1992-10-02';
var FECHA_FIN    = '2024-09-06';

var SST_MIN = -3;
var SST_MAX = 35;

var TEND_ANIO_INICIO = 1993;   // primer año completo HYCOM
var TEND_ANIO_FIN    = 2023;   // último año completo HYCOM

// Solo trabajamos con las bandas que necesitamos en este script
var BANDA_TEMP_REF = 'water_temp_0';  // para máscara térmica
var BANDAS_SAL = ['salinity_0', 'salinity_100'];


// ----------------------------------------------------------------------------
// 3. CARGAR Y PREPROCESAR HYCOM (réplica reducida del 07a)
// ----------------------------------------------------------------------------
// Cargamos solo las bandas que necesitamos para aliviar memoria
var hycom_raw = ee.ImageCollection('HYCOM/sea_temp_salinity')
  .filterDate(FECHA_INICIO, FECHA_FIN)
  .filterBounds(geom.AN_limpio)
  .select([BANDA_TEMP_REF].concat(BANDAS_SAL));

var hycom_diario = hycom_raw.map(function(img) {
  // Desescalado: factor 0.001, offset 20
  var temp_sup = img.select(BANDA_TEMP_REF).multiply(0.001).add(20)
                    .rename(BANDA_TEMP_REF);
  var sal = img.select(BANDAS_SAL).multiply(0.001).add(20).rename(BANDAS_SAL);

  // Máscara térmica sobre temperatura superficial
  var mascara = temp_sup.gte(SST_MIN).and(temp_sup.lte(SST_MAX));

  return sal.updateMask(mascara)
            .clip(geom.AN_limpio)
            .copyProperties(img, ['system:time_start']);
});


// ----------------------------------------------------------------------------
// 4. AGREGACIÓN MENSUAL Y ANUAL
// ----------------------------------------------------------------------------
var anio_ini = ee.Number(ee.Date(FECHA_INICIO).get('year'));
var anio_fin = ee.Number(ee.Date(FECHA_FIN).get('year'));
var meses_total = anio_fin.subtract(anio_ini).multiply(12).toInt();

// Colección mensual
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

// Agregación anual sobre el periodo de tendencias
var anios_serie = ee.List.sequence(TEND_ANIO_INICIO, TEND_ANIO_FIN);

var hycom_anual = ee.ImageCollection.fromImages(
  anios_serie.map(function(a) {
    var meses_anio = hycom_mensual.filter(ee.Filter.eq('anio', a));
    return meses_anio.mean()
      .set('anio', a)
      .set('system:time_start', ee.Date.fromYMD(a, 7, 1).millis());
  })
);


// ----------------------------------------------------------------------------
// 5. FUNCIÓN DE TENDENCIA POR REGIÓN Y BANDA
// ----------------------------------------------------------------------------
var calcular_tendencia = function(nombre_region, geometria, banda, etiqueta) {
  var serie = hycom_anual.map(function(img) {
    var v = img.select(banda).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geometria,
      scale: 25000,
      maxPixels: 1e9,
      bestEffort: true
    });
    return ee.Feature(null, {
      'anio': img.get('anio'),
      'val':  v.get(banda)
    });
  }).filter(ee.Filter.notNull(['val']))
    .sort('anio');

  var n = serie.size();

  var sen = serie.reduceColumns({
    reducer: ee.Reducer.sensSlope(),
    selectors: ['anio', 'val']
  });
  var mk = serie.reduceColumns({
    reducer: ee.Reducer.kendallsCorrelation(2),
    selectors: ['anio', 'val']
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
    'producto':                'HYCOM',
    'variable':                'salinidad',
    'profundidad':             etiqueta,
    'region':                  nombre_region,
    'n_anios':                 n,
    'tendencia_PSU_por_decada': ee.Number(sen.get('slope')).multiply(10),
    'tau_kendall':             tau,
    'p_valor':                 p_valor,
    'significativo_p005':      p_valor.lt(0.05)
  });
};


// ----------------------------------------------------------------------------
// 6. APLICAR A LAS SIETE REGIONES Y A LOS DOS NIVELES
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

// Construimos lista de features para 0 m
var feats_0m = regiones_lista.map(function(r) {
  return calcular_tendencia(r.nombre, r.geometria, 'salinity_0', '0 m');
});

// Construimos lista de features para 100 m
var feats_100m = regiones_lista.map(function(r) {
  return calcular_tendencia(r.nombre, r.geometria, 'salinity_100', '100 m');
});

// Combinamos las dos listas
var resultados = ee.FeatureCollection(feats_0m.concat(feats_100m));


// ----------------------------------------------------------------------------
// 7. INFORMACIÓN EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 07b1 — TENDENCIAS DE SALINIDAD HYCOM (0 m y 100 m)');
print('==================================================================');
print('Periodo:', TEND_ANIO_INICIO, 'a', TEND_ANIO_FIN, '(31 años)');
print('Niveles analizados: 0 m (superficie) y 100 m (subsuperficial)');
print('Unidades: PSU/década');
print('');
print('TABLA DE TENDENCIAS DE SALINIDAD HYCOM:');
print(resultados);