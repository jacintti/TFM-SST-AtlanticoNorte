// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 08b: ERA5-LAND — TENDENCIAS DECADALES POR ESTACIÓN
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Complemento del script 08 ampliado. Calcula las tendencias decadales
//   por ventana costera y estación, una estación cada vez para no saturar
//   la memoria de GEE. Resultados a integrar manualmente en una tabla
//   única tras tres ejecuciones.
//
// USO:
//   1. Ejecutar con ESTACION_A_CALCULAR = 'anual'
//   2. Ejecutar con ESTACION_A_CALCULAR = 'invierno'
//   3. Ejecutar con ESTACION_A_CALCULAR = 'verano'
//
//   Cada ejecución produce una FeatureCollection de 6 filas (una por ventana).
// ============================================================================


// ----------------------------------------------------------------------------
// 0. SELECTOR DE ESTACIÓN
// ----------------------------------------------------------------------------
// Cambia este valor entre 'anual', 'invierno', 'verano' y ejecuta cada vez.
var ESTACION_A_CALCULAR = 'verano';


// ----------------------------------------------------------------------------
// 1. IMPORTAR MÓDULO DE GEOMETRÍAS
// ----------------------------------------------------------------------------
var geom = require('users/jacintod/TFM1:00_geometrias_module');


// ----------------------------------------------------------------------------
// 2. PARÁMETROS
// ----------------------------------------------------------------------------
var FECHA_INICIO = '1981-12-01';
var FECHA_FIN    = '2025-03-01';

var TEND_ANIO_INICIO = 1982;
var TEND_ANIO_FIN    = 2024;


// ----------------------------------------------------------------------------
// 3. CARGAR ERA5-LAND MENSUAL
// ----------------------------------------------------------------------------
var era5_raw = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
  .filterDate(FECHA_INICIO, FECHA_FIN)
  .select('temperature_2m');

var era5_mensual = era5_raw.map(function(img) {
  var fecha = img.date();
  var anio_natural = ee.Number(fecha.get('year'));
  var mes = ee.Number(fecha.get('month'));
  var anio_DJF = ee.Algorithms.If(
    mes.eq(12),
    anio_natural.add(1),
    anio_natural
  );
  return img.subtract(273.15)
            .rename('temp_aire_C')
            .copyProperties(img, ['system:time_start'])
            .set('anio_natural', anio_natural)
            .set('anio_DJF', ee.Number(anio_DJF))
            .set('mes', mes);
});


// ----------------------------------------------------------------------------
// 4. CONSTRUCCIÓN DE LA SERIE ANUAL CORRESPONDIENTE A LA ESTACIÓN ACTIVA
// ----------------------------------------------------------------------------
var anios = ee.List.sequence(TEND_ANIO_INICIO, TEND_ANIO_FIN);

var construir_serie = function(estacion) {
  if (estacion === 'anual') {
    return ee.ImageCollection.fromImages(
      anios.map(function(a) {
        var imgs = era5_mensual.filter(ee.Filter.eq('anio_natural', a));
        return imgs.mean()
          .set('anio', a)
          .set('system:time_start', ee.Date.fromYMD(a, 7, 1).millis());
      })
    );
  } else if (estacion === 'invierno') {
    return ee.ImageCollection.fromImages(
      anios.map(function(a) {
        var imgs = era5_mensual.filter(
          ee.Filter.and(
            ee.Filter.eq('anio_DJF', a),
            ee.Filter.inList('mes', [12, 1, 2])
          )
        );
        return imgs.mean()
          .set('anio', a)
          .set('system:time_start', ee.Date.fromYMD(a, 1, 15).millis());
      })
    );
  } else {  // 'verano'
    return ee.ImageCollection.fromImages(
      anios.map(function(a) {
        var imgs = era5_mensual.filter(
          ee.Filter.and(
            ee.Filter.eq('anio_natural', a),
            ee.Filter.inList('mes', [6, 7, 8])
          )
        );
        return imgs.mean()
          .set('anio', a)
          .set('system:time_start', ee.Date.fromYMD(a, 7, 15).millis());
      })
    );
  }
};

var serie_activa = construir_serie(ESTACION_A_CALCULAR);


// ----------------------------------------------------------------------------
// 5. FUNCIÓN DE TENDENCIA POR VENTANA
// ----------------------------------------------------------------------------
var calcular_tendencia = function(nombre_ventana, geometria, costa, latitud) {
  var serie = serie_activa.map(function(img) {
    var v = img.select('temp_aire_C').reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geometria,
      scale: 11000,
      maxPixels: 1e9,
      bestEffort: true
    });
    return ee.Feature(null, {
      'anio': img.get('anio'),
      'temp': v.get('temp_aire_C')
    });
  }).filter(ee.Filter.notNull(['temp']))
    .sort('anio');

  var n = serie.size();

  var sen = serie.reduceColumns({
    reducer: ee.Reducer.sensSlope(),
    selectors: ['anio', 'temp']
  });
  var mk = serie.reduceColumns({
    reducer: ee.Reducer.kendallsCorrelation(2),
    selectors: ['anio', 'temp']
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
    'producto':                  'ERA5-Land',
    'estacion':                  ESTACION_A_CALCULAR,
    'ventana':                   nombre_ventana,
    'costa':                     costa,
    'latitud':                   latitud,
    'n_anios':                   n,
    'tendencia_C_por_decada':    ee.Number(sen.get('slope')).multiply(10),
    'tau_kendall':               tau,
    'p_valor':                   p_valor,
    'significativo_p005':        p_valor.lt(0.05)
  });
};


// ----------------------------------------------------------------------------
// 6. APLICAR A LAS SEIS VENTANAS
// ----------------------------------------------------------------------------
var ventanas_lista = [
  {nombre: 'L40w', geometria: geom.L40w, costa: 'Oeste (América)', lat: 40},
  {nombre: 'L40e', geometria: geom.L40e, costa: 'Este (Europa)',   lat: 40},
  {nombre: 'L53w', geometria: geom.L53w, costa: 'Oeste (América)', lat: 53},
  {nombre: 'L53e', geometria: geom.L53e, costa: 'Este (Europa)',   lat: 53},
  {nombre: 'L67w', geometria: geom.L67w, costa: 'Oeste (América)', lat: 67},
  {nombre: 'L67e', geometria: geom.L67e, costa: 'Este (Europa)',   lat: 67}
];

var resultados = ee.FeatureCollection(ventanas_lista.map(function(v) {
  return calcular_tendencia(v.nombre, v.geometria, v.costa, v.lat);
}));


// ----------------------------------------------------------------------------
// 7. INFORMACIÓN EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 08b — TENDENCIAS ESTACIONALES ERA5-LAND');
print('==================================================================');
print('Estación activa en esta ejecución:', ESTACION_A_CALCULAR);
print('Periodo:', TEND_ANIO_INICIO, 'a', TEND_ANIO_FIN);
print('');
print('TABLA — TENDENCIAS (' + ESTACION_A_CALCULAR.toUpperCase() + '):');
print(resultados);