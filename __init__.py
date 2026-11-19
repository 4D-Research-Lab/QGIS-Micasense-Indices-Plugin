"""Single-file MicaSense RedEdge-P indices plugin for QGIS.
Replace the plugin folder's __init__.py with this file; retain metadata.txt.
"""
import os
import csv
import math
import numpy as np
from qgis.PyQt.QtCore import QUrl
from qgis.PyQt.QtWidgets import (
    QAction, QCheckBox, QComboBox, QDialog, QFileDialog, QGridLayout,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QProgressBar,
    QPushButton, QScrollArea, QSpinBox, QVBoxLayout, QWidget
)
from qgis.core import QgsProject, QgsRasterLayer, QgsVectorLayer

NODATA = -9999.0

def safe_divide(a, b, fill=np.nan):
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        return np.where(b != 0, a / b, fill).astype(np.float32)

def iter_windows(width, height, tile_size):
    for y in range(0, height, tile_size):
        h = min(tile_size, height-y)
        for x in range(0, width, tile_size):
            yield x, y, min(tile_size, width-x), h

def read_band_tile(band, x, y, width, height):
    arr = band.ReadAsArray(x, y, width, height).astype(np.float32)
    nd = band.GetNoDataValue()
    if nd is not None:
        if np.isnan(nd): arr[np.isnan(arr)] = np.nan
        else: arr[arr == nd] = np.nan
    return arr

def create_output_dataset(path, reference_ds):
    from osgeo import gdal
    ds = gdal.GetDriverByName("GTiff").Create(
        path, reference_ds.RasterXSize, reference_ds.RasterYSize, 1,
        gdal.GDT_Float32,
        options=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"]
    )
    if ds is None: raise IOError(f"Could not create output raster: {path}")
    ds.SetGeoTransform(reference_ds.GetGeoTransform())
    ds.SetProjection(reference_ds.GetProjection())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(NODATA)
    band.Fill(NODATA)
    return ds

def write_tile(ds, array, x, y):
    values = np.where(np.isfinite(array), array, NODATA).astype(np.float32)
    ds.GetRasterBand(1).WriteArray(values, xoff=x, yoff=y)

def add_grouped_raster(path, name, group_name, groups):
    layer = QgsRasterLayer(path, name)
    if not layer.isValid(): raise RuntimeError(f"Could not load output layer: {path}")
    project = QgsProject.instance()
    project.addMapLayer(layer, False)
    root = project.layerTreeRoot()
    if group_name not in groups:
        groups[group_name] = root.findGroup(group_name)
        if groups[group_name] is None: groups[group_name] = root.addGroup(group_name)
    groups[group_name].addLayer(layer)

def pan_detail_tile(pan_band, full_width, full_height, x, y, width, height):
    """PAN high-pass tile with a one-pixel halo, avoiding tile-edge seams."""
    x0, y0 = max(0,x-1), max(0,y-1)
    x1, y1 = min(full_width,x+width+1), min(full_height,y+height+1)
    expanded = read_band_tile(pan_band,x0,y0,x1-x0,y1-y0)
    valid = np.isfinite(expanded)
    h,w = expanded.shape
    vals = np.pad(np.where(valid,expanded,0).astype(np.float32),1)
    weights = np.pad(valid.astype(np.float32),1)
    total = np.zeros((h,w),np.float32); count = np.zeros((h,w),np.float32)
    for dy in range(3):
        for dx in range(3):
            total += vals[dy:dy+h,dx:dx+w]
            count += weights[dy:dy+h,dx:dx+w]
    smooth=np.full((h,w),np.nan,np.float32)
    np.divide(total,count,out=smooth,where=count>0)
    detail=expanded-smooth; detail[~valid]=np.nan
    cx,cy=x-x0,y-y0
    return detail[cy:cy+height,cx:cx+width]

def pansharpen_tile(index, detail, strength, index_sd, detail_sd):
    valid=np.isfinite(index)&np.isfinite(detail)
    out=np.full(index.shape,np.nan,np.float32); out[valid]=index[valid]
    if index_sd>0 and detail_sd>0 and np.isfinite(detail_sd):
        out[valid]+=(strength*index_sd/detail_sd*detail[valid]).astype(np.float32)
    return out

def create_pan_grid_index(index_path, pan_ds, temp_path):
    """Warp one index file to PAN grid; GDAL processes the raster in blocks."""
    from osgeo import gdal
    src=gdal.Open(index_path,gdal.GA_ReadOnly)
    if src is None: raise IOError(f"Could not reopen index raster: {index_path}")
    if not src.GetProjection() or not pan_ds.GetProjection():
        src=None
        raise ValueError("Both index and PAN rasters need a defined CRS.")
    dst=create_output_dataset(temp_path,pan_ds)
    result=gdal.ReprojectImage(src,dst,src.GetProjection(),pan_ds.GetProjection(),gdal.GRA_Bilinear)
    src=None; dst.FlushCache(); dst=None
    if result not in (None,0): raise RuntimeError("GDAL could not align index to PAN grid.")

def pan_global_scales(index_band, pan_band, tile_windows, pan_width, pan_height):
    """Global scales from the index *band* and PAN detail; avoids tile seams."""
    total_n=0
    sum_i=sum_i2=sum_d=sum_d2=0.0
    # FIX: index_band is already a GDAL RasterBand. Do not call GetRasterBand on it.
    for x,y,w,h in tile_windows:
        index=read_band_tile(index_band,x,y,w,h)
        detail=pan_detail_tile(pan_band,pan_width,pan_height,x,y,w,h)
        valid=np.isfinite(index)&np.isfinite(detail)
        if not np.any(valid): continue
        iv=index[valid].astype(np.float64); dv=detail[valid].astype(np.float64)
        total_n+=iv.size
        sum_i+=iv.sum(); sum_i2+=np.dot(iv,iv)
        sum_d+=dv.sum(); sum_d2+=np.dot(dv,dv)
    if total_n<2: return 0.0,0.0
    mean_i,mean_d=sum_i/total_n,sum_d/total_n
    sd_i=math.sqrt(max(sum_i2/total_n-mean_i*mean_i,0.0))
    sd_d=math.sqrt(max(sum_d2/total_n-mean_d*mean_d,0.0))
    return sd_i,sd_d

def compute_index(entry,bands):
    result=entry[4](**bands)
    return np.where(np.isfinite(result),result,np.nan).astype(np.float32)

class PCAAccumulator:
    """Accumulate global covariance statistics from index tiles."""
    def __init__(self,n_features):
        self.n_features=n_features; self.count=0
        self.sums=np.zeros(n_features,np.float64)
        self.cross=np.zeros((n_features,n_features),np.float64)
    def add(self,arrays):
        block=np.column_stack([a.ravel() for a in arrays]).astype(np.float64,copy=False)
        valid=np.all(np.isfinite(block),axis=1); x=block[valid]
        if x.size:
            self.count+=x.shape[0]; self.sums+=x.sum(axis=0); self.cross+=x.T@x
    def finish(self,n_components):
        if self.count<2: raise ValueError("Not enough complete pixels for PCA.")
        means=self.sums/self.count
        cov=(self.cross-np.outer(self.sums,self.sums)/self.count)/(self.count-1)
        cov=(cov+cov.T)/2.0
        stds=np.sqrt(np.maximum(np.diag(cov),0))
        if np.any(stds==0): raise ValueError("A selected index has zero variance; remove it from PCA.")
        corr=cov/np.outer(stds,stds)
        eigvals,eigvecs=np.linalg.eigh(corr)
        order=np.argsort(eigvals)[::-1]
        eigvals=np.maximum(eigvals[order],0); eigvecs=eigvecs[:,order]
        total=eigvals.sum()
        if total<=0: raise ValueError("PCA found no variance in selected indices.")
        ratios=eigvals/total
        return means,stds,ratios[:n_components],eigvecs[:,:n_components],eigvals[:n_components]

def pca_score_tiles(arrays,means,stds,coefficients):
    block=np.column_stack([a.ravel() for a in arrays]).astype(np.float64,copy=False)
    valid=np.all(np.isfinite(block),axis=1)
    scores=[np.full(block.shape[0],np.nan,np.float32) for _ in range(coefficients.shape[1])]
    if np.any(valid):
        values=((block[valid]-means)/stds)@coefficients
        for c in range(coefficients.shape[1]): scores[c][valid]=values[:,c].astype(np.float32)
    shape=arrays[0].shape
    return [s.reshape(shape) for s in scores]

def write_loadings_csv(path,keys,ratios,coefficients,eigenvalues):
    cumulative=np.cumsum(ratios)
    with open(path,"w",newline="",encoding="utf-8-sig") as f:
        writer=csv.writer(f)
        writer.writerow(["Index","Family","Component","Eigenvector coefficient","Correlation loading","Explained variance (%)","Cumulative variance (%)"])
        for row,key in enumerate(keys):
            for c in range(len(ratios)):
                writer.writerow([key[1],key[0],f"PC{c+1}",float(coefficients[row,c]),float(coefficients[row,c]*np.sqrt(eigenvalues[c])),float(ratios[c]*100),float(cumulative[c]*100)])

# Index entries: short, long, category, formula, function(B,G,R,RE,NIR)
INDICES=[
("NDVI","Normalized Difference Vegetation Index","Vegetation / stress","(NIR-R)/(NIR+R)",lambda B,G,R,RE,NIR:safe_divide(NIR-R,NIR+R)),
("NDRE","Normalized Difference Red Edge","Vegetation / stress","(NIR-RE)/(NIR+RE)",lambda B,G,R,RE,NIR:safe_divide(NIR-RE,NIR+RE)),
("GNDVI","Green NDVI","Vegetation / stress","(NIR-G)/(NIR+G)",lambda B,G,R,RE,NIR:safe_divide(NIR-G,NIR+G)),
("CIre","Red Edge Chlorophyll Index","Vegetation / stress","NIR/RE-1",lambda B,G,R,RE,NIR:safe_divide(NIR,RE)-1),
("CIg","Green Chlorophyll Index","Vegetation / stress","NIR/G-1",lambda B,G,R,RE,NIR:safe_divide(NIR,G)-1),
("MCARI","Modified Chlorophyll Absorption Ratio Index","Vegetation / stress","MCARI",lambda B,G,R,RE,NIR:((RE-R)-.2*(RE-G))*safe_divide(RE,R)),
("MTVI2","Modified Triangular Vegetation Index 2","Vegetation / stress","MTVI2",lambda B,G,R,RE,NIR:1.5*(1.2*(NIR-G)-2.5*(R-G))/np.sqrt(np.maximum((2*NIR+1)**2-(6*NIR-5*np.sqrt(np.maximum(R,0)))-.5,1e-10))),
("REIP","Red Edge Inflection Point (linear approx.)","Vegetation / stress","REIP",lambda B,G,R,RE,NIR:700+40*safe_divide((R+NIR)/2-RE,NIR-RE)),
("NPCRI","Normalized Pigment Chlorophyll Ratio Index","Vegetation / stress","(R-B)/(R+B)",lambda B,G,R,RE,NIR:safe_divide(R-B,R+B)),
("PRI","Photochemical Reflectance Index (G/R proxy)","Vegetation / stress","(G-R)/(G+R)",lambda B,G,R,RE,NIR:safe_divide(G-R,G+R)),
("SIPI","Structure Insensitive Pigment Index","Vegetation / stress","(NIR-B)/(NIR-R)",lambda B,G,R,RE,NIR:safe_divide(NIR-B,NIR-R)),
("ARI","Anthocyanin Reflectance Index (proxy)","Vegetation / stress","1/G-1/RE",lambda B,G,R,RE,NIR:safe_divide(1.,G)-safe_divide(1.,RE)),
("PSRI","Plant Senescence Reflectance Index","Vegetation / stress","(R-B)/RE",lambda B,G,R,RE,NIR:safe_divide(R-B,RE)),
("NPCI","Normalized Pigment-Chlorophyll Index","Vegetation / stress","(R-B)/(R+B)",lambda B,G,R,RE,NIR:safe_divide(R-B,R+B)),
("GI","Greenness Index","Vegetation / stress","G/R",lambda B,G,R,RE,NIR:safe_divide(G,R)),
("RVI","Ratio Vegetation Index","Vegetation / stress","NIR/R",lambda B,G,R,RE,NIR:safe_divide(NIR,R)),
("DVI","Difference Vegetation Index","Vegetation / stress","NIR-R",lambda B,G,R,RE,NIR:NIR-R),
("IPVI","Infrared Percentage Vegetation Index","Vegetation / stress","NIR/(NIR+R)",lambda B,G,R,RE,NIR:safe_divide(NIR,NIR+R)),
("WDRVI","Wide Dynamic Range Vegetation Index (a=0.1)","Vegetation / stress","WDRVI",lambda B,G,R,RE,NIR:safe_divide(.1*NIR-R,.1*NIR+R)),
("SAVI","Soil-Adjusted Vegetation Index (L=0.5)","Soil / bare surface","SAVI",lambda B,G,R,RE,NIR:1.5*safe_divide(NIR-R,NIR+R+.5)),
("OSAVI","Optimised SAVI (L=0.16)","Soil / bare surface","OSAVI",lambda B,G,R,RE,NIR:safe_divide(NIR-R,NIR+R+.16)),
("MSAVI2","Modified SAVI 2","Soil / bare surface","MSAVI2",lambda B,G,R,RE,NIR:(2*NIR+1-np.sqrt(np.maximum((2*NIR+1)**2-8*(NIR-R),0)))/2),
("BI","Brightness Index","Soil / bare surface","BI",lambda B,G,R,RE,NIR:np.sqrt((R**2+G**2+NIR**2)/3)),
("SI","Soil Index","Soil / bare surface","(R-B)/(R+B)",lambda B,G,R,RE,NIR:safe_divide(R-B,R+B)),
("ColourIndex","Colour Index of Vegetation Extraction (CIVE)","Soil / bare surface","CIVE",lambda B,G,R,RE,NIR:.441*R-.881*G+.385*B+18.787),
("RI","Redness Index","Soil / bare surface","R^2/(B*G^3)",lambda B,G,R,RE,NIR:safe_divide(R**2,B*G**3)),
("SCI","Soil Colour Index","Soil / bare surface","(R-G)/(R+G)",lambda B,G,R,RE,NIR:safe_divide(R-G,R+G)),
("GSAVI","Green Soil-Adjusted Vegetation Index (L=0.5)","Soil / bare surface","GSAVI",lambda B,G,R,RE,NIR:1.5*safe_divide(NIR-G,NIR+G+.5)),
("TSAVI","Transformed SAVI","Soil / bare surface","TSAVI",lambda B,G,R,RE,NIR:safe_divide(.5*(NIR-.5*R-.1),.5*NIR+R-.05+.08*(1+.5**2))),
("NDWI","Normalized Difference Water Index (Green/NIR)","Moisture / organic","(G-NIR)/(G+NIR)",lambda B,G,R,RE,NIR:safe_divide(G-NIR,G+NIR)),
("MNDWI","Modified NDWI (Blue/NIR proxy)","Moisture / organic","(B-NIR)/(B+NIR)",lambda B,G,R,RE,NIR:safe_divide(B-NIR,B+NIR)),
("WI","Water Index (NIR/Red ratio)","Moisture / organic","NIR/R",lambda B,G,R,RE,NIR:safe_divide(NIR,R)),
("NDII","Normalized Difference Infrared Index (NIR/RE proxy)","Moisture / organic","(NIR-RE)/(NIR+RE)",lambda B,G,R,RE,NIR:safe_divide(NIR-RE,NIR+RE)),
("NDVIRE","NDVI Red Edge variant","Red-edge specific","(RE-R)/(RE+R)",lambda B,G,R,RE,NIR:safe_divide(RE-R,RE+R)),
("RENDVI","Red Edge NDVI (NIR/RE)","Red-edge specific","(NIR-RE)/(NIR+RE)",lambda B,G,R,RE,NIR:safe_divide(NIR-RE,NIR+RE)),
("SRRE","Simple Ratio Red Edge","Red-edge specific","NIR/RE",lambda B,G,R,RE,NIR:safe_divide(NIR,RE)),
("DATT","Datt Index (canopy N)","Red-edge specific","(NIR-RE)/(NIR-R)",lambda B,G,R,RE,NIR:safe_divide(NIR-RE,NIR-R)),
("RECI","Red Edge Chlorophyll Index (Gitelson)","Red-edge specific","NIR/RE-1",lambda B,G,R,RE,NIR:safe_divide(NIR,RE)-1),
("mNDVI705","Modified NDVI at 705/750 nm","Red-edge specific","(NIR-RE)/(NIR+RE-2*B)",lambda B,G,R,RE,NIR:safe_divide(NIR-RE,NIR+RE-2*B)),
("REChlIdx","Red Edge Chlorophyll Stress Index","Red-edge specific","R/RE",lambda B,G,R,RE,NIR:safe_divide(R,RE)),
("SR","Simple Ratio NIR/Red","Structural / contrast","NIR/R",lambda B,G,R,RE,NIR:safe_divide(NIR,R)),
("EVI","Enhanced Vegetation Index","Structural / contrast","EVI",lambda B,G,R,RE,NIR:2.5*safe_divide(NIR-R,NIR+6*R-7.5*B+1)),
("EVI2","Enhanced Vegetation Index 2","Structural / contrast","EVI2",lambda B,G,R,RE,NIR:2.5*safe_divide(NIR-R,NIR+2.4*R+1)),
("TVI","Transformed Vegetation Index","Structural / contrast","sqrt(NDVI+0.5)",lambda B,G,R,RE,NIR:np.sqrt(np.maximum(safe_divide(NIR-R,NIR+R)+.5,0))),
("VARI","Visible Atmospherically Resistant Index","Structural / contrast","(G-R)/(G+R-B)",lambda B,G,R,RE,NIR:safe_divide(G-R,G+R-B)),
("ExG","Excess Green Index","Structural / contrast","2*G-R-B",lambda B,G,R,RE,NIR:2*G-R-B),
("ExGR","Excess Green minus Excess Red","Structural / contrast","ExGR",lambda B,G,R,RE,NIR:(2*G-R-B)-(1.4*R-G)),
("GRVI","Green-Red Vegetation Index","Structural / contrast","(G-R)/(G+R)",lambda B,G,R,RE,NIR:safe_divide(G-R,G+R)),
("RGBVI","RGB Vegetation Index","Structural / contrast","RGBVI",lambda B,G,R,RE,NIR:safe_divide(G**2-B*R,G**2+B*R)),
("IKAW","Kawashima Index","Structural / contrast","(R-B)/(R+B)",lambda B,G,R,RE,NIR:safe_divide(R-B,R+B)),
("VEG","Vegetative Index (a=0.667)","Structural / contrast","VEG",lambda B,G,R,RE,NIR:safe_divide(G,np.power(np.maximum(R,0),.667)*np.power(np.maximum(B,0),.333))),
("NPCRI","Normalized Pigment Chlorophyll Ratio Index","Archaeological / cropmarks","(R-B)/(R+B)",lambda B,G,R,RE,NIR:safe_divide(R-B,R+B)),
("ARCVI","Archaeological Cropmark Vegetation Index (RE-NDVI)","Archaeological / cropmarks","(RE-R)/(RE+R)",lambda B,G,R,RE,NIR:safe_divide(RE-R,RE+R)),
("SNDVI","Soil-Normalised Difference Vegetation Index","Archaeological / cropmarks","(NIR-R)/(NIR+R+0.16)",lambda B,G,R,RE,NIR:safe_divide(NIR-R,NIR+R+.16)),
("CropRE","Cropmark Red Edge Ratio","Archaeological / cropmarks","RE/R",lambda B,G,R,RE,NIR:safe_divide(RE,R)),
("SoilBright","Soil Brightness Contrast Index","Archaeological / cropmarks","SoilBright",lambda B,G,R,RE,NIR:np.sqrt((R**2+G**2+B**2)/3)),
("NDVI_RE","NDVI Red Edge Difference","Archaeological / cropmarks","NDVI-NDRE",lambda B,G,R,RE,NIR:safe_divide(NIR-R,NIR+R)-safe_divide(NIR-RE,NIR+RE)),
("SoilLine","Soil Line Deviation Index","Archaeological / cropmarks","NIR-R-0.05",lambda B,G,R,RE,NIR:NIR-R-.05),
("PSRI_arch","Plant Senescence / Crop Stress Index","Archaeological / cropmarks","(R-B)/RE",lambda B,G,R,RE,NIR:safe_divide(R-B,RE)),
("GBR","Green-Blue Ratio","Archaeological / cropmarks","G/B",lambda B,G,R,RE,NIR:safe_divide(G,B)),
("RGR","Red-Green Ratio","Archaeological / cropmarks","R/G",lambda B,G,R,RE,NIR:safe_divide(R,G)),
("BGI","Blue-Green Index","Archaeological / cropmarks","B/G",lambda B,G,R,RE,NIR:safe_divide(B,G)),
("TNDVI","Transformed NDVI","Archaeological / cropmarks","sqrt(NDVI+0.5)",lambda B,G,R,RE,NIR:np.sqrt(np.maximum(safe_divide(NIR-R,NIR+R)+.5,0))),
("NDVI_B","NDVI Blue variant","Archaeological / cropmarks","(NIR-B)/(NIR+B)",lambda B,G,R,RE,NIR:safe_divide(NIR-B,NIR+B)),
]
CATEGORIES={}
for e in INDICES: CATEGORIES.setdefault(e[2],[]).append(e)

class IndicesDialog(QDialog):
    def __init__(self,iface,parent=None):
        super().__init__(parent); self.iface=iface
        self.setWindowTitle("MicaSense RedEdge-P Spectral Indices"); self.setMinimumSize(780,900)
        self.checkboxes={}; self.cat_labels={}; self._build_ui()
    def _build_ui(self):
        main=QVBoxLayout(); grp=QGroupBox("1. Multispectral raster"); grid=QGridLayout()
        grid.addWidget(QLabel("Project raster layer:"),0,0); self.layer_combo=QComboBox(); self._populate_layers()
        self.layer_combo.currentIndexChanged.connect(self._update_output_folder); grid.addWidget(self.layer_combo,0,1)
        grid.addWidget(QLabel("Expected bands: 1=Blue, 2=Green, 3=Red, 4=Red Edge, 5=NIR"),1,0,1,2)
        grid.addWidget(QLabel("Override band numbers:"),2,0); row=QHBoxLayout()
        self.spin_blue=QSpinBox(); self.spin_blue.setRange(1,99); self.spin_blue.setValue(1)
        self.spin_green=QSpinBox(); self.spin_green.setRange(1,99); self.spin_green.setValue(2)
        self.spin_red=QSpinBox(); self.spin_red.setRange(1,99); self.spin_red.setValue(3)
        self.spin_re=QSpinBox(); self.spin_re.setRange(1,99); self.spin_re.setValue(4)
        self.spin_nir=QSpinBox(); self.spin_nir.setRange(1,99); self.spin_nir.setValue(5)
        for label,spin in (("B",self.spin_blue),("G",self.spin_green),("R",self.spin_red),("RE",self.spin_re),("NIR",self.spin_nir)):
            row.addWidget(QLabel(label)); row.addWidget(spin)
        w=QWidget(); w.setLayout(row); grid.addWidget(w,2,1)
        grid.addWidget(QLabel("Output folder (input raster folder):"),3,0); self.out_dir=QLineEdit(); self.out_dir.setReadOnly(True); grid.addWidget(self.out_dir,3,1)
        grp.setLayout(grid); main.addWidget(grp)

        pgbox=QGroupBox("2. Optional Panchro source"); pg=QGridLayout(); pg.addWidget(QLabel("Source type:"),0,0)
        self.pan_mode=QComboBox(); self.pan_mode.addItems(["Project raster layer","GeoTIFF file"]); self.pan_mode.currentIndexChanged.connect(self._update_pan_controls); pg.addWidget(self.pan_mode,0,1,1,2)
        pg.addWidget(QLabel("Project raster:"),1,0); self.pan_layer=QComboBox(); self._populate_pan_layers(); pg.addWidget(self.pan_layer,1,1,1,2)
        pg.addWidget(QLabel("GeoTIFF file:"),2,0); self.pan_path=QLineEdit(); self.pan_browse=QPushButton("Browse…"); self.pan_browse.clicked.connect(self._browse_pan)
        pg.addWidget(self.pan_path,2,1); pg.addWidget(self.pan_browse,2,2)
        pg.addWidget(QLabel("Panchro band:"),3,0); self.spin_pan=QSpinBox(); self.spin_pan.setRange(1,99); self.spin_pan.setValue(1); pg.addWidget(self.spin_pan,3,1)
        pg.addWidget(QLabel("Detail strength (%):"),4,0); self.spin_strength=QSpinBox(); self.spin_strength.setRange(0,100); self.spin_strength.setValue(25); pg.addWidget(self.spin_strength,4,1)
        self.cb_pan=QCheckBox("Write PAN-detail-enhanced copies of the indices"); pg.addWidget(self.cb_pan,5,0,1,3)
        pgbox.setLayout(pg); main.addWidget(pgbox); self._update_pan_controls()

        grpidx=QGroupBox("3. Select vegetation indices"); vb=QVBoxLayout(); fr=QHBoxLayout(); fr.addWidget(QLabel("Filter category:"))
        self.cat_filter=QComboBox(); self.cat_filter.addItem("All categories")
        for cat in CATEGORIES: self.cat_filter.addItem(cat)
        self.cat_filter.currentTextChanged.connect(self._filter_category); fr.addWidget(self.cat_filter,1)
        allb=QPushButton("Select all"); noneb=QPushButton("Select none"); allb.clicked.connect(lambda:self._set_all(True)); noneb.clicked.connect(lambda:self._set_all(False)); fr.addWidget(allb); fr.addWidget(noneb); vb.addLayout(fr)
        scroll=QScrollArea(); scroll.setWidgetResizable(True); self.scroll_widget=QWidget(); self.scroll_layout=QVBoxLayout()
        for cat,entries in CATEGORIES.items():
            label=QLabel(f"<b>{cat}</b>"); self.scroll_layout.addWidget(label); self.cat_labels[cat]=label
            for short,long_name,_,formula,_fn in entries:
                cb=QCheckBox(f"{short} — {long_name}"); cb.setToolTip(f"Formula: {formula}"); self.checkboxes[(cat,short)]=cb; self.scroll_layout.addWidget(cb)
        self.scroll_layout.addStretch(); self.scroll_widget.setLayout(self.scroll_layout); scroll.setWidget(self.scroll_widget); scroll.setMinimumHeight(300); vb.addWidget(scroll); grpidx.setLayout(vb); main.addWidget(grpidx,1)

        opts=QGroupBox("4. Output, PCA and performance"); ov=QVBoxLayout()
        self.cb_add=QCheckBox("Add raster outputs to QGIS and group by family"); self.cb_add.setChecked(True); ov.addWidget(self.cb_add)
        self.cb_pca=QCheckBox("Run PCA on selected original indices"); ov.addWidget(self.cb_pca)
        pr=QHBoxLayout(); pr.addWidget(QLabel("PCA components:")); self.spin_pc=QSpinBox(); self.spin_pc.setRange(1,99); self.spin_pc.setValue(3); pr.addWidget(self.spin_pc); pr.addStretch(); ov.addLayout(pr)
        self.cb_pca_only=QCheckBox("PCA-only: write PCA rasters and loadings, not individual VI rasters"); self.cb_pca_only.toggled.connect(self._pca_only_toggled); ov.addWidget(self.cb_pca_only)
        self.cb_tiling=QCheckBox("Process in tiles to reduce memory use"); self.cb_tiling.setChecked(True); ov.addWidget(self.cb_tiling)
        tr=QHBoxLayout(); tr.addWidget(QLabel("Tile size (pixels):")); self.spin_tile=QSpinBox(); self.spin_tile.setRange(128,8192); self.spin_tile.setSingleStep(128); self.spin_tile.setValue(512); tr.addWidget(self.spin_tile); tr.addStretch(); ov.addLayout(tr); self.cb_tiling.toggled.connect(self.spin_tile.setEnabled)
        opts.setLayout(ov); main.addWidget(opts)
        self.progress=QProgressBar(); self.progress.setVisible(False); main.addWidget(self.progress); self.status=QLabel(""); main.addWidget(self.status)
        br=QHBoxLayout(); run=QPushButton("Run"); run.setDefault(True); run.clicked.connect(self.run); close=QPushButton("Close"); close.clicked.connect(self.reject); br.addWidget(run); br.addWidget(close); main.addLayout(br)
        self.setLayout(main); self._update_output_folder()
    def _populate_layers(self):
        self.layer_combo.clear()
        for l in QgsProject.instance().mapLayers().values():
            if isinstance(l,QgsRasterLayer): self.layer_combo.addItem(l.name(),l.id())
    def _populate_pan_layers(self):
        self.pan_layer.clear(); self.pan_layer.addItem("-- Select raster layer --",None)
        for l in QgsProject.instance().mapLayers().values():
            if isinstance(l,QgsRasterLayer): self.pan_layer.addItem(l.name(),l.id())
    @staticmethod
    def _path(layer):
        uri=layer.source().split("|",1)[0]
        if uri.startswith("file:"): uri=QUrl(uri).toLocalFile()
        return os.path.abspath(os.path.expanduser(uri)) if uri else ""
    def _update_output_folder(self,*_):
        l=QgsProject.instance().mapLayer(self.layer_combo.currentData()) if self.layer_combo.count() else None
        p=self._path(l) if l else ""; d=os.path.dirname(p) if p else ""; self.out_dir.setText(d if os.path.isdir(d) else "")
    def _update_pan_controls(self,*_):
        is_layer=self.pan_mode.currentIndex()==0; self.pan_layer.setEnabled(is_layer); self.pan_path.setEnabled(not is_layer); self.pan_browse.setEnabled(not is_layer)
    def _pca_only_toggled(self,on):
        if on: self.cb_pca.setChecked(True)
    def _browse_pan(self):
        p,_=QFileDialog.getOpenFileName(self,"Select Panchro GeoTIFF","","GeoTIFF (*.tif *.tiff);;All files (*.*)")
        if p: self.pan_path.setText(p)
    def _set_all(self,state):
        for cb in self.checkboxes.values():
            if cb.isVisible(): cb.setChecked(state)
    def _filter_category(self,text):
        allc=text=="All categories"
        for cat,entries in CATEGORIES.items():
            self.cat_labels[cat].setVisible(allc or cat==text)
            for short,*_ in entries: self.checkboxes[(cat,short)].setVisible(allc or cat==text)
    def _get_pan_source(self):
        if self.pan_mode.currentIndex()==0:
            lid=self.pan_layer.currentData(); l=QgsProject.instance().mapLayer(lid) if lid else None
            if not l or not l.isValid(): raise ValueError("Select a valid Panchro raster layer.")
            return l.source()
        p=self.pan_path.text().strip()
        if not p or not os.path.isfile(p): raise ValueError("Select a valid Panchro GeoTIFF.")
        return p
    def run(self):
        if not self.layer_combo.count(): QMessageBox.warning(self,"No raster","No raster layers found."); return
        source_layer=QgsProject.instance().mapLayer(self.layer_combo.currentData())
        if not source_layer or not source_layer.isValid(): QMessageBox.warning(self,"Invalid raster","Selected input is invalid."); return
        selected=[k for k,cb in self.checkboxes.items() if cb.isChecked()]
        if not selected: QMessageBox.warning(self,"No indices","Select at least one index."); return
        out_dir=self.out_dir.text().strip()
        if not out_dir or not os.path.isdir(out_dir): QMessageBox.warning(self,"Output folder unavailable","Select a local file-based multispectral raster."); return
        pca_only=self.cb_pca_only.isChecked(); do_pca=self.cb_pca.isChecked() or pca_only; use_pan=self.cb_pan.isChecked() and not pca_only
        if do_pca and len(selected)<2: QMessageBox.warning(self,"PCA selection","PCA requires at least two indices."); return
        if do_pca and self.spin_pc.value()>len(selected): QMessageBox.warning(self,"PCA components","Component count exceeds selected index count."); return
        try: pan_source=self._get_pan_source() if use_pan else None
        except Exception as e: QMessageBox.warning(self,"Panchro source",str(e)); return
        try:
            from osgeo import gdal
            src_ds=gdal.Open(source_layer.source(),gdal.GA_ReadOnly)
            if src_ds is None: raise IOError("GDAL could not open multispectral raster.")
            spins=(self.spin_blue,self.spin_green,self.spin_red,self.spin_re,self.spin_nir)
            if any(s.value()>src_ds.RasterCount for s in spins): raise ValueError("A requested spectral band exceeds the input band's count.")
            ms_bands=[src_ds.GetRasterBand(s.value()) for s in spins]
            pan_ds=None; pan_band=None
            if use_pan:
                pan_ds=gdal.Open(pan_source,gdal.GA_ReadOnly)
                if pan_ds is None: raise IOError("GDAL could not open selected Panchro source.")
                if self.spin_pan.value()>pan_ds.RasterCount: raise ValueError("Panchro band exceeds source band count.")
                pan_band=pan_ds.GetRasterBand(self.spin_pan.value())
        except Exception as e: QMessageBox.critical(self,"Input error",str(e)); return

        width,height=src_ds.RasterXSize,src_ds.RasterYSize
        tile_size=self.spin_tile.value() if self.cb_tiling.isChecked() else max(width,height)
        windows=list(iter_windows(width,height,tile_size))
        if not windows: QMessageBox.critical(self,"Empty raster","Input raster has no pixels."); return
        pca_acc=PCAAccumulator(len(selected)) if do_pca else None
        lookup={(e[2],e[0]):e for e in INDICES}; failed=set(); errors=[]; raw_paths={}; output_datasets={}; groups={}
        progress_value=0; expected=len(windows)*(1+(1 if do_pca else 0)+(2*len(selected) if use_pan else 0))
        self.progress.setVisible(True); self.progress.setMaximum(max(expected,1)); self.progress.setValue(0)
        def advance(message):
            nonlocal progress_value
            progress_value+=1; self.progress.setValue(min(progress_value,self.progress.maximum())); self.status.setText(message)
            from qgis.PyQt.QtWidgets import QApplication
            QApplication.processEvents()
        def token(text):
            s="".join(c.lower() if c.isalnum() else "_" for c in text)
            while "__" in s: s=s.replace("__","_")
            return s.strip("_")
        try:
            if not pca_only:
                for cat,short in selected:
                    path=os.path.join(out_dir,f"{token(cat)}_{short}.tif"); raw_paths[(cat,short)]=path; output_datasets[(cat,short)]=create_output_dataset(path,src_ds)
            # First tiled pass: compute indices, write original VI files, accumulate PCA moments.
            for ti,(x,y,w,h) in enumerate(windows,1):
                bands={name:read_band_tile(band,x,y,w,h) for name,band in zip(("B","G","R","RE","NIR"),ms_bands)}
                values=[]
                for key in selected:
                    try:
                        arr=compute_index(lookup[key],bands); values.append(arr)
                        if not pca_only: write_tile(output_datasets[key],arr,x,y)
                    except Exception as e:
                        failed.add(key); errors.append(f"{key[0]} — {key[1]}: {e}"); values.append(np.full((h,w),np.nan,np.float32))
                if do_pca:
                    if failed: raise RuntimeError("An index failed; PCA stopped to avoid incomplete inputs.")
                    pca_acc.add(values)
                advance(f"Calculating input tile {ti}/{len(windows)}")
            for ds in output_datasets.values(): ds.FlushCache()
            output_datasets.clear()
            if not pca_only:
                for key,path in raw_paths.items():
                    if key not in failed: add_grouped_raster(path,key[1],key[0],groups)

            # PAN stage: GDAL aligns each completed VI to PAN grid; processing remains tiled.
            if use_pan:
                pan_windows=list(iter_windows(pan_ds.RasterXSize,pan_ds.RasterYSize,tile_size)); strength=self.spin_strength.value()/100.0
                for key in selected:
                    if key in failed: continue
                    temp_path=os.path.join(out_dir,f"._{token(key[0])}_{key[1]}_pan_aligned_tmp.tif")
                    out_path=os.path.join(out_dir,f"{token(key[0])}_{key[1]}_PAN.tif")
                    temp_ds=out_ds=None
                    try:
                        create_pan_grid_index(raw_paths[key],pan_ds,temp_path)
                        temp_ds=gdal.Open(temp_path,gdal.GA_ReadOnly); out_ds=create_output_dataset(out_path,pan_ds)
                        ib=temp_ds.GetRasterBand(1)
                        sd_i,sd_d=pan_global_scales(ib,pan_band,pan_windows,pan_ds.RasterXSize,pan_ds.RasterYSize)
                        for ti,(x,y,w,h) in enumerate(pan_windows,1):
                            it=read_band_tile(ib,x,y,w,h)
                            detail=pan_detail_tile(pan_band,pan_ds.RasterXSize,pan_ds.RasterYSize,x,y,w,h)
                            write_tile(out_ds,pansharpen_tile(it,detail,strength,sd_i,sd_d),x,y)
                            advance(f"PAN sharpening {key[1]}: tile {ti}/{len(pan_windows)}")
                        out_ds.FlushCache(); out_ds=None; temp_ds=None
                        add_grouped_raster(out_path,key[1]+" [PAN]",key[0],groups)
                    finally:
                        out_ds=None; temp_ds=None
                        try:
                            if os.path.exists(temp_path): os.remove(temp_path)
                        except OSError: pass

            # PCA second pass: score rasters use global statistics from pass one.
            if do_pca:
                means,stds,ratios,coefficients,eigvals=pca_acc.finish(self.spin_pc.value())
                pc_paths=[os.path.join(out_dir,f"PCA_PC{i+1}.tif") for i in range(self.spin_pc.value())]
                pc_datasets=[create_output_dataset(p,src_ds) for p in pc_paths]
                try:
                    for ti,(x,y,w,h) in enumerate(windows,1):
                        bands={name:read_band_tile(band,x,y,w,h) for name,band in zip(("B","G","R","RE","NIR"),ms_bands)}
                        arrays=[compute_index(lookup[k],bands) for k in selected]
                        scores=pca_score_tiles(arrays,means,stds,coefficients)
                        for ds,score in zip(pc_datasets,scores): write_tile(ds,score,x,y)
                        advance(f"PCA scoring tile {ti}/{len(windows)}")
                    for ds in pc_datasets: ds.FlushCache()
                finally: pc_datasets.clear()
                for i,path in enumerate(pc_paths): add_grouped_raster(path,f"PC{i+1} ({ratios[i]*100:.1f}%)","PCA",groups)
                load_path=os.path.join(out_dir,"PCA_loadings.csv")
                write_loadings_csv(load_path,selected,ratios,coefficients,eigvals)
                if self.cb_add.isChecked():
                    table=QgsVectorLayer(load_path,"PCA loadings (CSV table)","ogr")
                    if table.isValid():
                        QgsProject.instance().addMapLayer(table,False)
                        if "PCA" not in groups:
                            root=QgsProject.instance().layerTreeRoot(); groups["PCA"]=root.findGroup("PCA") or root.addGroup("PCA")
                        groups["PCA"].addLayer(table)
                    else: errors.append("PCA_loadings.csv was saved but could not be loaded as a QGIS table.")
        except Exception as e:
            errors.append(str(e))
        finally:
            output_datasets.clear(); src_ds=None
            if pan_ds is not None: pan_ds=None
        self.progress.setValue(self.progress.maximum()); self.status.setText("Done.")
        msg=f"Output folder:\n{out_dir}"
        if pca_only: msg+="\nPCA-only mode: individual VI rasters were not written."
        if use_pan: msg+="\nPAN-enhanced rasters use the selected PAN grid."
        if do_pca: msg+="\nPCA score rasters and PCA_loadings.csv were written."
        if errors:
            msg+="\n\nErrors:\n"+"\n".join(errors); QMessageBox.warning(self,"Completed with errors",msg)
        else: QMessageBox.information(self,"Done",msg)
        self.progress.setVisible(False); self.status.setText("")

class MicaSenseIndicesPlugin:
    def __init__(self,iface): self.iface=iface; self.action=None
    def initGui(self):
        self.action=QAction("MicaSense Spectral Indices…",self.iface.mainWindow()); self.action.triggered.connect(self.run)
        self.iface.addPluginToMenu("MicaSense Indices",self.action); self.iface.addToolBarIcon(self.action)
    def unload(self):
        if self.action:
            self.iface.removePluginMenu("MicaSense Indices",self.action); self.iface.removeToolBarIcon(self.action)
    def run(self): IndicesDialog(self.iface,self.iface.mainWindow()).exec()

def classFactory(iface): return MicaSenseIndicesPlugin(iface)
