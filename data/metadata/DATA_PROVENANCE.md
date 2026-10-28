# Data provenance

The public analysis table was assembled from processed factor rasters and Great Wall heritage attributes. A local audit compared a random sample of 300 assessment units with the model rasters. WE, WAE, FT, SA, CL, NDVI, RND, GI and HRF matched the audited rasters for all 300 sampled units. The earthquake and landslide fields were traced to separate density rasters rather than the older copies stored with the nine other factors.

## Source and period notes

| Factor | Period represented in the model audit | Source entry point | Note |
| --- | --- | --- | --- |
| EQ | 1990–2024; records through 2024-03-27 | https://data.earthquake.cn/index.html | Density surface derived from the study-area earthquake records |
| LS | Not confirmed from local metadata | http://www.gisrs.cn/infofordata?id=52c30b65-7214-422e-9e89-eeb101b49561 | Density surface; do not assign a year without further source evidence |
| WE | 2015 | https://www.resdc.cn/data.aspx?DATAID=363 | Wind-erosion source product |
| WAE | 2015 | https://www.resdc.cn/data.aspx?DATAID=363 | Water-erosion source product |
| FT | Baseline period not confirmed | https://data.tpdc.ac.cn/zh-hans/data/fa9263cb-70fb-4b51-9b2b-6d671b0fbeb2 | Thematic freeze–thaw susceptibility product, not treated as a single-year state |
| SA | 2005 | https://www.geodata.cn/main/face_science_detail?typeName=face_science&guid=186416065572720 | Salinization source product |
| NDVI | August 2012 | https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MOD13A2 | Monthly maximum-value composite, not an annual mean |
| HRF | 2012 | https://data.tpdc.ac.cn/zh-hans/data/4753b961-6441-4c42-b671-26a7b72e6029 | Heavy-rainfall event count |
| CL | 2012 | https://zenodo.org/records/12779975 | Cropland extracted from the land-cover product |
| RND | Strongly supported as 2014; complete generation history unavailable | https://download.geofabrik.de/asia/china.html | Derived from the road network; OpenStreetMap attribution and ODbL terms apply |
| GI | 2012 | https://www.nesdc.org.cn/sdo/detail?id=672aea607e28174998e63234 | Grazing-intensity product |

The final 1 × 1 km assessment units differ from the native raster cells. Raster processing and assignment to the assessment grid do not increase the effective resolution of the source products.

## Heritage attributes

Wall material, preservation status and protection designation were compiled from the China Great Wall Heritage website: http://www.greatwallheritage.cn/CCMCMS/html/1/index.html.

## Redistribution

The processed tabular values are released with this repository. The source and processed rasters are not part of the GitHub distribution. Before redistributing raster archives, confirm the terms of every source provider and retain all required attribution. OpenStreetMap-derived road data must comply with ODbL 1.0.
