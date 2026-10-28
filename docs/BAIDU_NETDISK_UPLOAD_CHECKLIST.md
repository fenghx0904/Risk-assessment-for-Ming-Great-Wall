# Raw-raster archive checklist

The raw-raster download link is pending. After the archive is uploaded, add the Baidu Netdisk URL, extraction code, archive size and SHA-256 here.

## Recommended archive contents

Create one archive named `ming_great_wall_11_factor_rasters.zip` containing exactly the model-used raster files below. The paths are relative to the local factor-data directory.

| Public archive name | Local source file | Factor |
| --- | --- | --- |
| `01_EQ_earthquake_density.tif` | `地震/研究区地震核密度.tif` | EQ |
| `02_LS_landslide_density.tif` | `滑坡/滑坡核密度1.tif` | LS |
| `03_WE_wind_erosion.tif` | `2012年/删除0重采样/cWinderosion.tif` | WE |
| `04_WAE_water_erosion.tif` | `2012年/删除0重采样/dRainerosion.tif` | WAE |
| `05_FT_freeze_thaw.tif` | `2012年/删除0重采样/eFreeze-thaw.tif` | FT |
| `06_SA_salinization.tif` | `2012年/删除0重采样/fSalinization.tif` | SA |
| `07_CL_cropland.tif` | `2012年/删除0重采样/gCultivatedland.tif` | CL |
| `08_NDVI_vegetation.tif` | `2012年/删除0重采样/hVegetationgrowth.tif` | NDVI |
| `09_RND_road_density.tif` | `2012年/删除0重采样/iRoad.tif` | RND |
| `10_GI_grazing_intensity.tif` | `2012年/删除0重采样/jGrazingintensity.tif` | GI |
| `11_HRF_heavy_rainfall.tif` | `2012年/删除0重采样/kHeavyrainfall.tif` | HRF |

Do not substitute similarly named files from the `2024年` folder. The audit of the current model input traced the nine listed factors to `2012年/删除0重采样`, while EQ and LS came from the separate density rasters listed above.

## Optional spatial-support archive

Only if sharing detailed heritage geometry is permitted, create `ming_great_wall_spatial_support.zip` with the 1 km assessment-grid shapefile components and the DEM used for elevation-sensitivity analysis. Keep this separate from the 11-factor archive. Do not publish precise heritage locations if institutional, legal or conservation restrictions apply.

## Before upload

1. Confirm redistribution permission for every third-party source product.
2. Retain source-provider citations and license files.
3. Exclude ArcGIS caches, overviews, temporary files, duplicated rasters and old candidate versions.
4. Exclude every text file containing FTP, account, token or password information.
5. Record the final archive size and SHA-256 after compression.
