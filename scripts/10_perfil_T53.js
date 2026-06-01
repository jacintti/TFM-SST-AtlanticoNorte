// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 10: PERFIL LONGITUDINAL DE SST A LO LARGO DEL TRANSECTO T53
// ----------------------------------------------------------------------------
// PROPÓSITO:
//   Generar un perfil longitudinal de la SST media a lo largo del paralelo
//   53°N (transecto Terranova-Irlanda), comparando dos décadas:
//     · Década inicial: 1982-1991
//     · Década final:   2015-2024
//
//   El gráfico resultante muestra simultáneamente:
//     (a) el cambio gradual de SST de oeste (Terranova) a este (Irlanda),
//     (b) la diferencia entre las dos décadas en cada punto del transecto.
//
//   La zona central del transecto corresponde al cold blob, donde se espera
//   un calentamiento menor o nulo entre las dos décadas.
// ============================================================================


// ----------------------------------------------------------------------------
// 1. PARÁMETROS
// ----------------------------------------------------------------------------
var LAT_TRANSECTO = 53.0;       // paralelo central del transecto
var LON_OESTE     = -55.0;      // Terranova
var LON_ESTE      = -10.0;      // Irlanda
var N_PUNTOS      = 46;         // un punto por grado de longitud aprox.

var SST_MIN = -3;
var SST_MAX = 35;


// ----------------------------------------------------------------------------
// 2. CARGAR OISST Y CONSTRUIR MEDIAS DECADALES
// ----------------------------------------------------------------------------
var oisst = ee.ImageCollection('NOAA/CDR/OISST/V2_1').select('sst');

var mediaDecada = function(anioIni, anioFin) {
  var col = oisst
    .filterDate(anioIni + '-01-01', (anioFin + 1) + '-01-01')
    .map(function(img) {
      var sst = img.select('sst').multiply(0.01);
      var mask = sst.gte(SST_MIN).and(sst.lte(SST_MAX));
      return sst.updateMask(mask);
    });
  return col.mean().rename('SST');
};

var sst_8291 = mediaDecada(1982, 1991);
var sst_1524 = mediaDecada(2015, 2024);


// ----------------------------------------------------------------------------
// 3. CONSTRUIR PUNTOS A LO LARGO DEL TRANSECTO
// ----------------------------------------------------------------------------
var puntos = ee.List.sequence(0, N_PUNTOS - 1).map(function(i) {
  i = ee.Number(i);
  var frac = i.divide(N_PUNTOS - 1);
  var lon = ee.Number(LON_OESTE).add(
    frac.multiply(ee.Number(LON_ESTE).subtract(LON_OESTE))
  );
  var punto = ee.Geometry.Point([lon, LAT_TRANSECTO]);
  return ee.Feature(punto, {'lon': lon});
});
var puntosFC = ee.FeatureCollection(puntos);


// ----------------------------------------------------------------------------
// 4. EXTRAER SST EN CADA PUNTO PARA LAS DOS DÉCADAS
// ----------------------------------------------------------------------------
var imagenDoble = sst_8291.rename('SST_1982_1991')
  .addBands(sst_1524.rename('SST_2015_2024'));

var perfil = imagenDoble.reduceRegions({
  collection: puntosFC,
  reducer: ee.Reducer.first(),
  scale: 25000
});


// ----------------------------------------------------------------------------
// 5. GRÁFICO DEL PERFIL
// ----------------------------------------------------------------------------
var grafico = ui.Chart.feature.byFeature(perfil, 'lon',
    ['SST_1982_1991', 'SST_2015_2024'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Perfil de SST a lo largo del transecto T53 (53°N): Terranova → Irlanda',
    hAxis: {
      title: 'Longitud (grados; oeste → este)',
      format: '#.#'
    },
    vAxis: {title: 'SST media (°C)'},
    series: {
      0: {color: '#1f77b4', lineWidth: 2, pointSize: 3},  // 1982-1991 azul
      1: {color: '#d62728', lineWidth: 2, pointSize: 3}    // 2015-2024 rojo
    },
    legend: {position: 'top'},
    interpolateNulls: true
  });

print('Perfil longitudinal del transecto T53:');
print(grafico);


// ----------------------------------------------------------------------------
// 6. VISUALIZACIÓN DEL TRANSECTO EN EL MAPA (verificación)
// ----------------------------------------------------------------------------
Map.centerObject(puntosFC, 4);
Map.addLayer(puntosFC, {color: 'red'}, 'Puntos del transecto T53');
Map.addLayer(
  ee.Geometry.LineString([[LON_OESTE, LAT_TRANSECTO], [LON_ESTE, LAT_TRANSECTO]]),
  {color: 'yellow'}, 'Línea del transecto'
);