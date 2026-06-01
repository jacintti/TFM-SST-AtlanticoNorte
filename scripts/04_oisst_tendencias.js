// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 04: NOAA OISST v2.1 — TENDENCIAS POR MANN-KENDALL Y SEN'S SLOPE
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Cuantificar el calentamiento (o enfriamiento) del Atlántico Norte por
//   región y píxel a píxel mediante:
//     · Test de Mann-Kendall (significación estadística).
//     · Estimador de pendiente de Sen (Theil-Sen slope, en °C/década).
//
//   Trabajamos sobre series anuales, no mensuales, por consistencia con la
//   literatura climática multidecadal (Caesar et al. 2018, Biguino et al. 2023,
//   Khatib & Chaaban 2025, entre otros) y por viabilidad computacional.
//
// SALIDA:
//   · Tabla de tendencias por región (consola).
//   · Mapa píxel a píxel de tendencias (°C/década).
//   · Mapa de significación estadística (p < 0.05 destacado).
//
// REFERENCIAS METODOLÓGICAS:
//   Mann, H.B. (1945) — test no paramétrico de tendencia.
//   Kendall, M.G. (1975) — desarrollo del estadístico tau.
//   Sen, P.K. (1968) — estimador de pendiente robusto.
//   Hamed, K.H. & Rao, A.R. (1998) — corrección por autocorrelación.
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

// Para el cálculo de tendencia anual usamos años completos:
var ANIO_INICIO = 1982;   // primer año completo OISST
var ANIO_FIN    = 2024;   // último año completo

// Paleta divergente azul-blanco-rojo para tendencia (°C/década)
var PALETA_TENDENCIA = [
  '053061','2166ac','4393c3','92c5de','d1e5f0',
  'ffffff',
  'fddbc7','f4a582','d6604d','b2182b','67001f'
];


// ----------------------------------------------------------------------------
// 3. RECONSTRUIR COLECCIÓN MENSUAL OISST (réplica de scripts 01 y 02)
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
// 4. AGREGACIÓN A SERIE ANUAL (MEDIA DE LOS 12 MESES)
// ----------------------------------------------------------------------------
// Construimos una colección donde cada imagen es la SST media anual.
var anios_serie = ee.List.sequence(ANIO_INICIO, ANIO_FIN);

var oisst_anual = ee.ImageCollection.fromImages(
  anios_serie.map(function(a) {
    var meses_del_anio = oisst_mensual.filter(ee.Filter.eq('anio', a));
    var media_anual = meses_del_anio.mean()
      .set('anio', a)
      // banda 'year' como número para el reductor de tendencia
      .addBands(ee.Image.constant(a).toFloat().rename('year'));
    return media_anual.set('system:time_start',
      ee.Date.fromYMD(a, 7, 1).millis());  // 1 julio como fecha de referencia
  })
);


// ----------------------------------------------------------------------------
// 5. CÁLCULO DE TENDENCIA PÍXEL A PÍXEL (PENDIENTE DE SEN)
// ----------------------------------------------------------------------------
// Reordenamos las bandas: primero 'year' (independiente), después 'sst'
// (dependiente). sensSlope espera ese orden.
var oisst_para_sen = oisst_anual.map(function(img) {
  return img.select(['year', 'sst']);
});

// Calcular la pendiente de Sen píxel a píxel
var tendencia_decadal = oisst_para_sen
  .reduce(ee.Reducer.sensSlope())
  .select('slope')
  .multiply(10)   // de °C/año a °C/década
  .rename('tendencia_decadal');


// ----------------------------------------------------------------------------
// 6. CÁLCULO DE SIGNIFICACIÓN MANN-KENDALL PÍXEL A PÍXEL
// ----------------------------------------------------------------------------
var significancia_mk = oisst_anual.select('sst')
  .reduce(ee.Reducer.kendallsCorrelation(1));

// El reductor devuelve dos bandas: 'sst_tau' y 'sst_p-value'
// Renombramos para claridad
significancia_mk = significancia_mk.rename(['tau_kendall', 'p_valor']);

// Máscara de significación: solo píxeles con p < 0.05
var significativo = significancia_mk.select('p_valor').lt(0.05);


// ----------------------------------------------------------------------------
// 7. TENDENCIA POR REGIÓN (TABLA)
// ----------------------------------------------------------------------------
// Para cada región, calculamos:
//   · pendiente de Sen anual sobre la SST media regional anual
//   · tau de Kendall y p-valor

var calcular_tendencia_region = function(nombre, geometria) {
  // Extraer la SST media anual de la región como FeatureCollection
  var serie_anual_region = oisst_anual.map(function(img) {
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

  var n = serie_anual_region.size();

  // Sen's slope
  var sen = serie_anual_region.reduceColumns({
    reducer: ee.Reducer.sensSlope(),
    selectors: ['anio', 'sst']
  });

  // Mann-Kendall (2 inputs)
  var mk = serie_anual_region.reduceColumns({
    reducer: ee.Reducer.kendallsCorrelation(2),
    selectors: ['anio', 'sst']
  });

  var tau = ee.Number(mk.get('tau'));

  // Cálculo manual de Z y p-valor de Mann-Kendall
  // Z = tau * sqrt[9n(n-1) / (2(2n+5))]
  // (Fórmula clásica del estadístico Mann-Kendall normalizado)
  var n_num = ee.Number(n);
  var z_denominador = ee.Number(2).multiply(n_num.multiply(2).add(5));
  var z_factor = n_num.multiply(9).multiply(n_num.subtract(1)).divide(z_denominador);
  var z = tau.multiply(z_factor.sqrt());

  // p-valor a dos colas usando aproximación de la normal estándar
  // p = 2 * (1 - Φ(|Z|))
  // Aproximamos Φ(|Z|) con la fórmula de Abramowitz-Stegun (aprox. 26.2.17)
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

  // Indicador de significancia
  var significativo = p_valor.lt(0.05);

  return ee.Feature(null, {
    'region':                 nombre,
    'n_anios':                n,
    'tendencia_C_por_decada': ee.Number(sen.get('slope')).multiply(10),
    'tau_kendall':            tau,
    'z_mk':                   z,
    'p_valor':                p_valor,
    'significativo_p005':     significativo
  });
};

// Aplicar a las regiones que tenemos definidas
var resultados = ee.FeatureCollection([
  calcular_tendencia_region('AN_limpio', geom.AN_limpio),
  calcular_tendencia_region('GSE',       geom.GSE),
  calcular_tendencia_region('LAB',       geom.LAB),
  calcular_tendencia_region('NOR',       geom.NOR),
  calcular_tendencia_region('EUR',       geom.EUR),
  calcular_tendencia_region('CBL',       geom.CBL),
  calcular_tendencia_region('T53',       geom.T53)
]);

print('==================================================================');
print('SCRIPT 04 — TENDENCIAS POR REGIÓN (1982-2024, OISST anual)');
print('==================================================================');
print('Tabla de tendencias:');
print(resultados);


// ----------------------------------------------------------------------------
// 8. VISUALIZACIÓN: MAPA PÍXEL A PÍXEL
// ----------------------------------------------------------------------------
Map.centerObject(geom.AN_limpio, 3);
Map.setOptions('HYBRID');

// Mapa 1: tendencia decadal (°C/década)
Map.addLayer(
  tendencia_decadal,
  {min: -0.5, max: 0.5, palette: PALETA_TENDENCIA},
  '01. Tendencia decadal SST 1982-2024 (°C/década)',
  true
);

// Mapa 2: máscara de píxeles con tendencia significativa (p < 0.05)
Map.addLayer(
  significativo.selfMask(),
  {min: 0, max: 1, palette: ['000000']},
  '02. Píxeles con p < 0.05 (significativos)',
  false
);


// ----------------------------------------------------------------------------
// 9. INFORMACIÓN DIAGNÓSTICA
// ----------------------------------------------------------------------------
print('Banda tendencia píxel a píxel calculada.');
print('Las cifras de la tabla están en °C/década.');
print('p < 0.05 → tendencia estadísticamente significativa.');
print('Valores positivos = calentamiento; negativos = enfriamiento.');

// ----------------------------------------------------------------------------
// 10. ANÁLISIS POR SUBPERIODOS
// ----------------------------------------------------------------------------
// Calculamos las mismas tendencias pero divididas en dos subperiodos:
//   · 1982-2000 (~19 años): periodo histórico previo a la aceleración del 
//     debilitamiento del cold blob según la literatura reciente.
//   · 2001-2024 (~24 años): periodo donde el cold blob se ha documentado
//     más intensamente.
//
// El objetivo es detectar si alguna región (especialmente CBL) muestra
// una desaceleración o reversión de tendencia entre los dos subperiodos.
// ----------------------------------------------------------------------------

var ANIO_CORTE = 2000;  // último año del subperiodo 1 (inclusivo)

var calcular_tendencia_region_subperiodo = function(nombre, geometria, 
                                                     anio_min, anio_max) {
  // Filtrar la serie anual al subperiodo
  var serie_anual_sub = oisst_anual
    .filter(ee.Filter.gte('anio', anio_min))
    .filter(ee.Filter.lte('anio', anio_max));

  var serie_features = serie_anual_sub.map(function(img) {
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

  // Cálculo Z y p-valor (mismo procedimiento que la función principal)
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

  var significativo = p_valor.lt(0.05);

  return ee.Feature(null, {
    'region':                 nombre,
    'subperiodo':             ee.String(ee.Number(anio_min).int())
                                 .cat('-')
                                 .cat(ee.String(ee.Number(anio_max).int())),
    'n_anios':                n,
    'tendencia_C_por_decada': ee.Number(sen.get('slope')).multiply(10),
    'tau_kendall':            tau,
    'p_valor':                p_valor,
    'significativo_p005':     significativo
  });
};

// Aplicar a las siete regiones para ambos subperiodos
var regiones_lista = [
  {nombre: 'AN_limpio', geometria: geom.AN_limpio},
  {nombre: 'GSE',       geometria: geom.GSE},
  {nombre: 'LAB',       geometria: geom.LAB},
  {nombre: 'NOR',       geometria: geom.NOR},
  {nombre: 'EUR',       geometria: geom.EUR},
  {nombre: 'CBL',       geometria: geom.CBL},
  {nombre: 'T53',       geometria: geom.T53}
];

var resultados_sub1 = ee.FeatureCollection(regiones_lista.map(function(r) {
  return calcular_tendencia_region_subperiodo(r.nombre, r.geometria, 1982, ANIO_CORTE);
}));

var resultados_sub2 = ee.FeatureCollection(regiones_lista.map(function(r) {
  return calcular_tendencia_region_subperiodo(r.nombre, r.geometria, ANIO_CORTE + 1, 2024);
}));

print('==================================================================');
print('SUBPERIODO 1: 1982-' + ANIO_CORTE);
print('==================================================================');
print(resultados_sub1);

print('==================================================================');
print('SUBPERIODO 2: ' + (ANIO_CORTE + 1) + '-2024');
print('==================================================================');
print(resultados_sub2);
