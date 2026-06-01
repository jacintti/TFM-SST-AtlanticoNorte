// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 08: ERA5-LAND — MEDIAS DECADALES DE TEMPERATURA POR VENTANA Y ESTACIÓN
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Calcular la temperatura media de cada ventana costera, por estación
//   (anual, invierno DJF, verano JJA) y por década:
//     · 1982-1989, 1990s, 2000s, 2010s, 2020-2024
//
//   El cálculo de tendencias Mann-Kendall + Sen se realiza separadamente
//   en el script 08b para evitar problemas de memoria.
//
// SALIDA:
//   · Una tabla con 90 filas (6 ventanas × 3 estaciones × 5 décadas) con
//     la temperatura media correspondiente.
// ============================================================================


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

var DECADAS = [
  {nombre: '1982-1989', ini: 1982, fin: 1989},
  {nombre: '1990s',     ini: 1990, fin: 1999},
  {nombre: '2000s',     ini: 2000, fin: 2009},
  {nombre: '2010s',     ini: 2010, fin: 2019},
  {nombre: '2020-2024', ini: 2020, fin: 2024}
];


// ----------------------------------------------------------------------------
// 3. CARGAR ERA5-LAND MENSUAL Y PREPARAR ETIQUETAS
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
// 4. CONSTRUCCIÓN DE SERIES ANUAL, INVIERNO Y VERANO
// ----------------------------------------------------------------------------
var anios = ee.List.sequence(TEND_ANIO_INICIO, TEND_ANIO_FIN);

var era5_anual = ee.ImageCollection.fromImages(
  anios.map(function(a) {
    var imgs = era5_mensual.filter(ee.Filter.eq('anio_natural', a));
    return imgs.mean()
      .set('anio', a)
      .set('system:time_start', ee.Date.fromYMD(a, 7, 1).millis());
  })
);

var era5_invierno = ee.ImageCollection.fromImages(
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

var era5_verano = ee.ImageCollection.fromImages(
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


// ----------------------------------------------------------------------------
// 5. FUNCIÓN DE CÁLCULO DE MEDIAS DECADALES POR VENTANA Y ESTACIÓN
// ----------------------------------------------------------------------------
var calcular_media_decada = function(nombre_ventana, geometria, costa, latitud,
                                      coleccion, etiqueta_estacion) {
  return DECADAS.map(function(d) {
    var imgs = coleccion
      .filter(ee.Filter.gte('anio', d.ini))
      .filter(ee.Filter.lte('anio', d.fin));
    var media_imagen = imgs.mean();
    var v = media_imagen.select('temp_aire_C').reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geometria,
      scale: 11000,
      maxPixels: 1e9,
      bestEffort: true
    });
    return ee.Feature(null, {
      'producto':       'ERA5-Land',
      'estacion':       etiqueta_estacion,
      'ventana':        nombre_ventana,
      'costa':          costa,
      'latitud':        latitud,
      'decada':         d.nombre,
      'decada_ini':     d.ini,
      'decada_fin':     d.fin,
      'temperatura_C':  v.get('temp_aire_C')
    });
  });
};


// ----------------------------------------------------------------------------
// 6. APLICAR A LAS SEIS VENTANAS Y A LAS TRES SERIES
// ----------------------------------------------------------------------------
var ventanas_lista = [
  {nombre: 'L40w', geometria: geom.L40w, costa: 'Oeste (América)', lat: 40},
  {nombre: 'L40e', geometria: geom.L40e, costa: 'Este (Europa)',   lat: 40},
  {nombre: 'L53w', geometria: geom.L53w, costa: 'Oeste (América)', lat: 53},
  {nombre: 'L53e', geometria: geom.L53e, costa: 'Este (Europa)',   lat: 53},
  {nombre: 'L67w', geometria: geom.L67w, costa: 'Oeste (América)', lat: 67},
  {nombre: 'L67e', geometria: geom.L67e, costa: 'Este (Europa)',   lat: 67}
];

var medias_anual = ventanas_lista.map(function(v) {
  return calcular_media_decada(v.nombre, v.geometria, v.costa, v.lat,
                                era5_anual, 'anual');
});
var medias_invierno = ventanas_lista.map(function(v) {
  return calcular_media_decada(v.nombre, v.geometria, v.costa, v.lat,
                                era5_invierno, 'invierno');
});
var medias_verano = ventanas_lista.map(function(v) {
  return calcular_media_decada(v.nombre, v.geometria, v.costa, v.lat,
                                era5_verano, 'verano');
});

var aplanar = function(listaDeListas) {
  return ee.List(listaDeListas).flatten();
};

var tabla_medias = ee.FeatureCollection(
  ee.List(aplanar(medias_anual))
    .cat(aplanar(medias_invierno))
    .cat(aplanar(medias_verano))
);


// ----------------------------------------------------------------------------
// 7. INFORMACIÓN EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 08 — ERA5-LAND MEDIAS DECADALES POR VENTANA Y ESTACIÓN');
print('==================================================================');
print('Periodo:', TEND_ANIO_INICIO, 'a', TEND_ANIO_FIN, '(43 años)');
print('Décadas: 1982-1989, 1990s, 2000s, 2010s, 2020-2024');
print('Estaciones: anual, invierno (DJF), verano (JJA)');
print('');
print('TABLA — MEDIAS DECADALES (90 filas: 6 ventanas × 3 estaciones × 5 décadas):');
print(tabla_medias);