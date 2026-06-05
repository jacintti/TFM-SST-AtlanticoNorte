# TFM-SST-AtlanticoNorte
Scripts de Google Earth Engine para el análisis temporal de la SST en el Atlántico Norte (1981-2024). TFM, Máster en TIG, UCM, 2026
# Análisis temporal de la SST en el Atlántico Norte (1981-2024)

Repositorio de scripts de Google Earth Engine asociado al Trabajo Fin de Máster:

> Domínguez Morales, J. (2026). *Análisis temporal de la temperatura superficial del mar en el Atlántico Norte (1981-2024): evolución, asimetrías regionales e influencia del Atlántico sobre el clima europeo*. Máster en Tecnologías de la Información Geográfica, UCM.

## Productos analizados

- **NOAA OISST v2.1** (1982-2024): temperatura superficial del mar.
- **MODIS Aqua L3SMI** (2003-2021): validación cruzada satelital.
- **HYCOM** (1993-2023): análisis tridimensional de temperatura y salinidad.
- **ERA5-Land** (1982-2024): temperatura del aire en ventanas costeras.

## Estructura del repositorio

La carpeta `scripts/` contiene los códigos de GEE numerados según el flujo de trabajo:

- **00**: módulo de geometrías reutilizable.
- **01-04**: análisis OISST (carga, climatologías, anomalías, tendencias, subperiodos, exportación).
- **05-06**: análisis MODIS (validación cruzada).
- **07**: análisis HYCOM (temperatura y salinidad por profundidad).
- **08**: análisis ERA5-Land (ventanas costeras).
- **09**: exportación de productos cartográficos finales.
- **10**: perfil longitudinal de la SST media a lo largo del paralelo 53°N.
## Uso

Los scripts están diseñados para ejecutarse en Google Earth Engine Code Editor (https://code.earthengine.google.com). El módulo 00 se importa mediante `require()` desde los demás scripts.

## Licencia

MIT License. Si usas este código, cita el trabajo original.

## Contacto

Jacinto Domínguez Morales — jacintod@ucm.es
