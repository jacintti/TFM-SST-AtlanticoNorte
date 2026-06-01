// ============================================================================
// TFM - JACINTO DOMÍNGUEZ MORALES - 2026
// Análisis Temporal de la SST en el Atlántico Norte (1981-2024)
// ============================================================================
// SCRIPT 01: NOAA OISST v2.1 — CARGA Y PREPROCESAMIENTO 
// ----------------------------------------------------------------------------
// CORRECCIÓN v2:
//   GEE no aplica automáticamente el factor de escala 0.01 a la banda 'sst'
//   de OISST v2.1. Hay que multiplicar manualmente por 0.01 después de cargar
//   la banda. Sin esta operación, los valores estaban en centésimas de grado
//   (1750 en lugar de 17.5 °C), lo que hacía que la máscara térmica
//   descartase todos los píxeles del océano abierto.
//
// CATÁLOGO GEE: NOAA/CDR/OISST/V2_1
// REFERENCIA: Huang et al. (2021). J. Climate, 34(8), 2923-2939.
//             https://doi.org/10.1175/JCLI-D-20-0166.1
// ============================================================================


// ----------------------------------------------------------------------------
// 1. IMPORTAR MÓDULO DE GEOMETRÍAS
// ----------------------------------------------------------------------------
var geom = require('users/jacintod/TFM1:00_geometrias_module');


// ----------------------------------------------------------------------------
// 2. PARÁMETROS DEL SCRIPT
// ----------------------------------------------------------------------------
var FECHA_INICIO = '1981-09-01';
var FECHA_FIN    = '2025-01-01';

var SST_MIN = -3;
var SST_MAX = 35;

var PALETA_SST = [
  '040274','040281','0502a3','0502b8','0502ce','0502e6','0602ff',
  '235cb1','307ef3','269db1','30c8e2','32d3ef','3be285','3ff38f',
  '86e26f','3ae237','b5e22e','d6e21f','fff705','ffd611','ffb613',
  '17f80e','ff8b13','ff6e08','ff500d','ff0000','de0101','c21301','a71001','911003'
];


// ----------------------------------------------------------------------------
// 3. CARGAR PRODUCTO Y APLICAR FACTOR DE ESCALA
// ----------------------------------------------------------------------------
// La banda 'sst' viene en centésimas de grado Celsius (entero corto).
// Multiplicamos por 0.01 para obtener °C.
var oisst_raw = ee.ImageCollection('NOAA/CDR/OISST/V2_1')
  .filterDate(FECHA_INICIO, FECHA_FIN)
  .filterBounds(geom.AN_limpio)
  .select('sst');


// ----------------------------------------------------------------------------
// 4. PREPROCESAMIENTO POR IMAGEN: ESCALA + MÁSCARA TÉRMICA + RECORTE
// ----------------------------------------------------------------------------
var oisst_diario = oisst_raw.map(function(img) {
  // Aplicar factor de escala (centésimas de grado → grados centígrados)
  var sst = img.select('sst').multiply(0.01).rename('sst');
  // Máscara térmica: conservar solo píxeles con SST en rango físico válido
  var mascara = sst.gte(SST_MIN).and(sst.lte(SST_MAX));
  return sst.updateMask(mascara)
            .clip(geom.AN_limpio)
            .copyProperties(img, ['system:time_start']);
});


// ----------------------------------------------------------------------------
// 5. AGREGACIÓN A ESCALA MENSUAL (MEDIANA DE LOS DÍAS DEL MES)
// ----------------------------------------------------------------------------
var anio_inicio = ee.Number(ee.Date(FECHA_INICIO).get('year'));
var anio_fin    = ee.Number(ee.Date(FECHA_FIN).get('year'));
var meses_total = anio_fin.subtract(anio_inicio).multiply(12).toInt();

var oisst_mensual = ee.ImageCollection.fromImages(
  ee.List.sequence(0, meses_total.subtract(1)).map(function(offset) {
    var fecha_mes = ee.Date(FECHA_INICIO).advance(ee.Number(offset), 'month');
    var fecha_fin_mes = fecha_mes.advance(1, 'month');
    var imgs_del_mes = oisst_diario.filterDate(fecha_mes, fecha_fin_mes);
    var num_imagenes = imgs_del_mes.size();
    var mes_compuesto = imgs_del_mes.median()
      .set('system:time_start', fecha_mes.millis())
      .set('anio', fecha_mes.get('year'))
      .set('mes', fecha_mes.get('month'))
      .set('n_imagenes', num_imagenes);
    return mes_compuesto;
  })
)
.filter(ee.Filter.gt('n_imagenes', 0));


// ----------------------------------------------------------------------------
// 6. INFORMACIÓN DIAGNÓSTICA EN CONSOLA
// ----------------------------------------------------------------------------
print('==================================================================');
print('SCRIPT 01 — OISST v2.1 CARGA Y PREPROCESAMIENTO (v2)');
print('==================================================================');
print('Periodo solicitado:', FECHA_INICIO, '→', FECHA_FIN);
print('Imágenes diarias cargadas:', oisst_diario.size());
print('Imágenes mensuales generadas:', oisst_mensual.size());
print('Primer mes con dato:', ee.Date(
  oisst_mensual.aggregate_min('system:time_start')).format('YYYY-MM'));
print('Último mes con dato:', ee.Date(
  oisst_mensual.aggregate_max('system:time_start')).format('YYYY-MM'));

// SST media de junio 2010 sobre AN_limpio (verificación)
var junio_2010 = oisst_mensual.filterDate('2010-06-01', '2010-07-01').first();
var sst_jun2010 = junio_2010.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: geom.AN_limpio,
  scale: 25000,
  maxPixels: 1e9,
  bestEffort: true
});
print('SST media junio 2010 sobre AN_limpio (esperado ~17-18°C):', sst_jun2010);


// ----------------------------------------------------------------------------
// 7. VISUALIZACIÓN: MEDIA CLIMATOLÓGICA DEL PERIODO COMPLETO
// ----------------------------------------------------------------------------
Map.centerObject(geom.AN_limpio, 3);
Map.setOptions('HYBRID');

var sst_media_periodo = oisst_mensual.mean();

Map.addLayer(
  sst_media_periodo,
  {min: 0, max: 28, palette: PALETA_SST},
  '01. OISST — SST media 1981-2024 (°C)',
  true
);

// Contornos de las geometrías regionales (capas ocultas por defecto)
var pintar_contorno = function(g, color, nombre) {
  Map.addLayer(
    ee.Image().paint(ee.FeatureCollection([ee.Feature(g)]), 0, 2),
    {palette: [color]}, nombre, false
  );
};
pintar_contorno(geom.AN_limpio, 'white',  '02. AN_limpio (contorno)');
pintar_contorno(geom.GSE,       'red',    '03. GSE (contorno)');
pintar_contorno(geom.LAB,       'blue',   '04. LAB (contorno)');
pintar_contorno(geom.NOR,       'purple', '05. NOR (contorno)');
pintar_contorno(geom.EUR,       'orange', '06. EUR (contorno)');
