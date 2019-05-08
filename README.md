# land-use-land-cover-classification

Preprocessing and Classification of landsat 8 satellite images using Random Forest Classifier.

Preprocessing includes:
- False color composite creation.
- DN conversion to TOA Reflectance.
- Atmospheric Correction using Dark Object Subtraction.
- Cloud Correction using time series data.
- Color Balancing using CLAHE (Contrast Local Adaptive Histogram Equalisation) and contrast stretching.
- Reprojection to WGS84 coordinates.
- Mosaicing and clipping using shapefile.
- Rasterize training shapefiles necessary for classification.
