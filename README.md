# MicaSense RedEdge-P Spectral Indices for QGIS

A QGIS plugin for calculating spectral indices from MicaSense RedEdge-P multispectral orthomosaics. It can group output rasters by index family, optionally inject spatial detail from a separate panchromatic (Panchro) raster, run standardized PCA on selected indices, and process large rasters in tiles.

**Compatibility:** This plugin works only with QGIS 4.

## Features

- Calculate selected indices from a five-band multispectral raster.
- Select bands explicitly; the default order is Blue, Green, Red, Red Edge, NIR.
- Group output rasters by index family in the QGIS Layers panel.
- Optionally create PAN-detail-enhanced index rasters from a separate PAN GeoTIFF or a raster layer already loaded in QGIS.
- Run PCA on the selected original indices and export component rasters and a loadings table.
- Use PCA-only mode to write PCA outputs without exporting every individual index raster.
- Process raster data in tiles to reduce peak memory use.
- Automatically write results into the folder containing the selected multispectral raster.

## Index families

The plugin groups indices into the categories used in its selection dialog:

- **Vegetation / stress:** NDVI, NDRE, GNDVI, CIre, CIg, MCARI, MTVI2, REIP, NPCRI, PRI, SIPI, ARI, PSRI, NPCI, GI, RVI, DVI, IPVI, WDRVI
- **Soil / bare surface:** SAVI, OSAVI, MSAVI2, BI, SI, ColourIndex, RI, SCI, GSAVI, TSAVI
- **Moisture / organic:** NDWI, MNDWI, WI, NDII
- **Red-edge specific:** NDVIRE, RENDVI, SRRE, DATT, RECI, mNDVI705, REChlIdx
- **Structural / contrast:** SR, EVI, EVI2, TVI, VARI, ExG, ExGR, GRVI, RGBVI, IKAW, VEG
- **Archaeological / cropmarks:** NPCRI, ARCVI, SNDVI, CropRE, SoilBright, NDVI_RE, SoilLine, PSRI_arch, GBR, RGR, BGI, TNDVI, NDVI_B

Some entries are proxies or adapted indices rather than direct equivalents of indices requiring bands the RedEdge-P does not record. For example, MNDWI uses Blue as a proxy because the sensor does not include a SWIR band. Check each formula and its intended interpretation before analysis.

## Requirements

- QGIS 4 with its Python plugin environment.
- GDAL and NumPy available in that QGIS environment. These are commonly available through QGIS, but availability can depend on the installation.
- A local, file-based multispectral raster that GDAL can open. The output folder is derived from this raster's file location.

The code is a single-file QGIS plugin module. Keep the plugin's `metadata.txt` in its plugin directory.

## Installation

1. Create a plugin directory named `micasense_indices` in the QGIS profile's Python plugins folder.
2. Save the plugin module as `__init__.py` inside that directory.
3. Ensure the directory also contains a valid `metadata.txt` file.
4. Restart QGIS 4, open **Plugins → Manage and Install Plugins**, and enable **MicaSense Indices**.
5. Open the plugin from the Plugins menu or its toolbar icon.

If QGIS reports a syntax error while loading the plugin, check that `__init__.py` contains Python code only. R syntax such as `<-` or `commandArgs()` does not belong in a QGIS Python plugin.

## Input data

### Multispectral raster

The default band mapping is:

| Plugin band | RedEdge-P band |
|---|---|
| B | Blue |
| G | Green |
| R | Red |
| RE | Red Edge |
| NIR | NIR |

The dialog allows band numbers to be changed. Confirm the band order in the raster's metadata or properties rather than relying only on the default.

Reflectance-calibrated input is recommended. Several formulas contain additive constants that assume reflectance values are on a 0–1 scale. If the raster contains unscaled digital numbers, index values for those formulas may not be meaningful without adapting the formulas or scaling the data. Ratios and normalized-difference indices may also be affected by data scaling when NoData, offsets, or other transformations are involved.

### Panchro raster

PAN detail enhancement is optional. Choose either a raster layer already loaded in the QGIS project or a GeoTIFF file, then select the band containing Panchro data. The PAN source and multispectral raster need defined coordinate reference systems so GDAL can align them.

## Running the plugin

1. Add the multispectral raster to the QGIS project and open the plugin.
2. Select the multispectral layer and confirm or adjust its five band numbers.
3. Select the indices to calculate. The category filter and **Select all / Select none** buttons help manage the list.
4. Optionally enable PAN-detail enhancement and select the PAN source, band, and detail strength.
5. Optionally enable PCA, choose the number of components, and select PCA-only mode if individual VI files are not needed.
6. Leave tiled processing enabled for large rasters, choose a tile size, and run the plugin.

### Output location

Outputs are written to the same folder as the selected multispectral raster. The output-folder field is read-only. The plugin requires a local file-based raster to determine this folder.

### Tiled processing

Tiling is enabled by default with a tile size of 512 × 512 pixels. Smaller tiles reduce peak memory use; larger tiles can improve throughput when sufficient memory is available. Disabling tiling processes the raster as a single window, which can require substantial memory for large rasters or many selected indices.

Tiling reduces the size of arrays held in memory, but does not guarantee a faster run. PCA requires a statistics pass followed by a scoring pass. PAN enhancement also reads the aligned index and Panchro data in tiles.

## Outputs

### Index rasters

When PCA-only mode is off, each selected index is written as a single-band Float32 GeoTIFF. Filenames include the index family to distinguish repeated short names, for example:

- `vegetation_stress_NDVI.tif`
- `archaeological_cropmarks_NPCRI.tif`

When enabled, PAN-detail-enhanced versions are written with `_PAN` added to the filename and use the PAN raster's grid. The original index raster remains at the multispectral raster's grid.

### PAN-detail enhancement

This implementation is a high-pass detail-injection procedure, not a fully radiometrically constrained pansharpening method. It bilinearly aligns the index to the PAN grid, subtracts a 3 × 3 local mean from the PAN image to form a detail layer, and injects a globally scaled version of that detail into the index. A one-pixel tile halo is used for the local filter, and global scale estimates are used to reduce visible seams between tiles.

The detail-strength control scales the injected detail. At 0%, no PAN detail is injected; the index is still aligned/resampled to the PAN grid. Increasing the setting adds stronger local PAN contrast and changes the index values. The enhanced raster should therefore be treated as a derived product, not as an unchanged spectral index. Inspect it alongside the original and test more than one strength.

### PCA rasters and loadings table

PCA uses the selected **original** indices, not PAN-enhanced copies. Inputs are standardized, so PCA is based on their correlation structure rather than their original numeric ranges. Pixels with NoData in any selected index are excluded from the PCA calculations and are NoData in component rasters.

The plugin writes:

- `PCA_PC1.tif`, `PCA_PC2.tif`, and so on, at the multispectral raster's grid. The QGIS layer names include each component's explained-variance percentage.
- `PCA_loadings.csv`, a UTF-8 CSV table that opens in spreadsheet software and can be loaded into QGIS as a non-spatial table.

The loadings table is in long format, with one row per input-index/component pair. It includes the eigenvector coefficient, correlation loading, explained variance, and cumulative explained variance. The coefficient is the component weight used to calculate scores from standardized indices. The correlation loading is the correlation between an input index and a component. Component signs are arbitrary, so interpret the pattern of loadings rather than the sign alone.

PCA is exploratory: a component or a high score is not, by itself, evidence of an archaeological feature. Highly correlated or duplicate indices can give a shared pattern additional influence.

### QGIS layer groups

When **Add raster outputs to QGIS** is enabled, index rasters are added under their index-family groups. PCA component rasters and the loadings table are added under `PCA`.

## Notes and limitations

- The plugin computes indices from the selected raster bands; it does not perform reflectance calibration.
- The PAN output is limited to the coverage and grid of the PAN raster. Areas without overlapping valid multispectral data remain NoData.
- PAN detail injection modifies index values and can emphasize PAN noise, edges, or non-vegetation features.
- PCA uses complete-case pixels: a pixel is used only if every selected index has a valid value there.
- Output files with the same names are created in the input folder. Before rerunning, close/remove existing output layers if the operating system or GDAL prevents overwriting files currently in use.
- Processing many indices or large tiles can still require significant memory and disk space.

## Troubleshooting

### Plugin fails to load with a Python syntax error

Check the exact `__init__.py` path named in the QGIS log. Replace the file with the plugin's Python module; remove any non-Python text, Markdown fences, or R code.

### A requested band is outside the raster's band count

Check the raster's actual band count and update the band selectors. A separate single-band Panchro file normally uses band 1.

### PAN alignment fails

Confirm both sources have defined CRSs and spatial overlap. For project layers, use a GDAL-readable raster source.

### Processing runs out of memory or becomes unresponsive

Keep tiling enabled and reduce the tile size. If memory is available, a larger tile size may reduce processing overhead. PCA and PAN processing require additional passes through the data.

### PAN-enhanced output is missing in PCA-only mode

PCA-only mode writes PCA products and loadings, not individual index or PAN-enhanced index rasters. Turn off PCA-only mode to export those rasters.

## Interpreting results

Index values depend on formula, reflectance scaling, and data quality. PCA component scores are not vegetation-index values and may be positive or negative. Review the original index maps, component rasters, and loadings table together, and document the selected indices and processing settings for reproducibility.

## License

This project is licensed under the Apache License 2.0. See the `LICENSE` file for the full license text.

## Contributing and issue reports

When reporting a problem, include the QGIS version, operating system, raster dimensions and band count, selected indices, tile size, whether PCA or PAN enhancement was enabled, and the full error message from the QGIS log. Avoid sharing sensitive project data.
