from osgeo import gdal
import pandas as pd, xarray as xr, os
import scipy.misc as sm, subprocess
import geopandas as gpd
from facets import facets
from grab_meta import grab_meta
import rasterio, numpy as np, glob
from rasterio.merge import merge 
from rasterio import plot 
from l8qa.qa import *
from PIL import Image
import matplotlib.pyplot as plt
import fiona, rasterio.mask
from skimage import exposure 
from rasterio.plot import reshape_as_raster, reshape_as_image


def norm(band):
    band = band.astype('float')
    band = (band-band.min())/(band.max()-band.min())
    return band

def reproject_dataset(filename):
    """Project a GeoTIFF to the WGS84 coordinate reference system.
    See https://mapbox.github.io/rasterio/topics/reproject.html"""

    # We want to project the GeoTIFF coordinate reference system (crs)
    # to WGS84 (e.g. into the familiar Lat/Lon pairs). WGS84 is analogous
    # to EPSG:4326
    dst_crs = 'EPSG:4326'

    with rasterio.open(filename) as src:
        transform, width, height = rasterio.warp.calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height,
            'compress': 'lzw'
        })

        out_path = filename.split('.tif')[0] + '.wgs84.tif'
        with rasterio.open(out_path, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                rasterio.warp.reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=rasterio.warp.Resampling.nearest)

        return out_path

class Landsat():
    def __init__(self, files):
        self.infiles = files
    
    def sort_files(self):
        files = self.infiles
        #x = int(files[0].split('.')[0].split('/')[-1].split('_')[-1][1:])
        files = sorted(files, key = lambda x: x.split('.')[1].\
                       split('/')[-1].split('_')[-1])
        lst = []
        
        for filename in files:
            try:
                x = int(filename.split('.')[0].\
                        split('/')[-1].split('_')[-1][1:])
                lst.append(filename)
            except:
                #print filename
                pass
            
        self.filenames = sorted(lst, key = lambda x: int(x.split('.')[0].\
                                split('/')[-1].split('_')[-1][1:]))
        return self
    
    def layer_stacking(self, ofile='stack.tif'):
        file_list = self.filenames
        self.stackfile = ofile
        # Read metadata of first file
        with rasterio.open(file_list[0]) as src0:
            meta = src0.meta

        # Update meta to reflect the number of layers
        meta.update(count = len(file_list), compress='lzw')

        # Read each layer and write it to stack
        with rasterio.open(ofile, 'w', **meta) as dst:
            for id, layer in enumerate(file_list, start=1):
                with rasterio.open(layer) as src1:
                    dst.write_band(id, src1.read(1))
        return self

    def composite(self, band = [5, 4, 3], ofile='fcc.tif', cmin=2, cmax=98):
        self.composite_file = ofile
        dataset = rasterio.open(self.stackfile)
        meta = dataset.meta
        meta.update(count = len(band), compress='lzw')
        #print meta

        red = dataset.read([band[0]])
        red = (255*self.norm(red)).astype('int')
        green = dataset.read([band[1]])
        green = (255*self.norm(green)).astype('int')
        blue = dataset.read([band[2]])
        blue = (255*self.norm(blue)).astype('int')

        rgb = np.vstack((red, green, blue)).astype('uint16')
        
        with rasterio.open(self.composite_file, "w", **meta) as dest:
                dest.write(rgb)
        return self
    
    @staticmethod                
    def norm(band):
        band = band.astype('float')
        band = (band-band.min())/(band.max()-band.min())
        return band
    
    @staticmethod
    def imshow(filename):
        temp = 'temp.tif'
        src = rasterio.open(filename)
        data = np.squeeze(src.read())
        sm.toimage(data).save(temp)
        
        Image.MAX_IMAGE_PIXELS = data.size + 100
        im = Image.open(temp)
        im.show()
               
class AtmosCorrection():
    """ Correct for atmospheric scatter (haze) using Dark Object Subtraction

    Given a red, green and blue band and model of atmospheric
    scattering, correct the bands to remove this scatter component.

    The Rayleigh model states that relative scattering is
    inversely proportional to the fourth power of the wavelength (-4)

    The Mie model uses the first power for moderate atmospheric
    conditions (-1)

    Assumes input is TOA reflectance, float32, 0..1
    Outputs a uint16 RGB, 0..55000
    """
    def __init__(self, infile):
        self.infile = infile        
        
    def correct(self, ofile, model = -2, \
                scale_factor = 55000, dtype = 'float32'):
        """ model ranges from -4 (Very Clear) to -0.5 (Very Hazy)
        """
        with rasterio.open(self.infile) as src:
            bands = src.read()
            profile = src.profile
            
        profile.update(dtype=dtype, count=3, compress='lzw')            
        self.dark = self.find_dark_object_value(bands[0])
        hazes = self.calculate_haze_reduction(self.dark, 'red', model) 
        
        corrected = []        
        for band, haze in zip(bands, hazes):
            corr = band - haze
            corr[corr < 0] = 0
            corr = norm(corr * scale_factor)
            corrected.append(corr)
            
        with rasterio.open(ofile, "w", **profile) as dst:
            for id, layer in enumerate(corrected, start=1):
                dst.write_band(id, layer.astype(profile['dtype']))
        return self        	
            
    @staticmethod        
    def find_dark_object_value(arr):
        """Find the value of the dark object
        ie the first dark value that is not an outlier
        """
        preval = None
        step = arr.max() #/ 255.0
        for val in np.unique(arr)[:100]:
            if val == 0:
                continue
            if preval is not None and (val - preval) < step:
                break
            else:
                preval = val
        return preval
    
    @staticmethod     
    def calculate_haze_reduction(dark, band, model, mult=1.0):
        """Calculate the haze reduction for red, blue and green bands
        given the a visible band's TOA reflectance of a dark object
        and a scattering model (inverse proportion)

        Returns a three element tuple containing th TOA values to
        substract from the R, G and B bands to account for scatter
        """
        # wavelengths, micrometers
        colors = ['red', 'green', 'blue']
        wavelengths = {'blue': 0.48, 'green': 0.56, 'red': 0.655}
        haze = [dark * ((wavelengths[color] ** model) / \
                        (wavelengths[band]  ** model)) * mult for color in colors]
        return haze
        
class CloudCorrection():
    def __init__(self, file1, mask, file2):
        self.file1 = file1
        self.mask1 = mask
        self.file2 = file2

    def cloud_correction(self, band=1):
        orig = rasterio.open(self.file1).read(band)
        mask = rasterio.open(self.mask1).read(1)
        othr = rasterio.open(self.file2).read(band)
        corr = orig*(~mask) + othr*mask
        return corr/255.0  

    def correct(self, out = 'ex.tif', dtype = 'float32'):   
        cloud_r = self.cloud_correction(band=1)
        cloud_g = self.cloud_correction(band=2)
        cloud_b = self.cloud_correction(band=3)
        rgb = [cloud_r, cloud_g, cloud_b]

        meta = rasterio.open(self.file1).meta
        meta.update(dtype=dtype, count=3, compress='lzw')
        with rasterio.open(out, "w", **meta) as dst:
            for id, layer in enumerate(rgb, start=1):
                dst.write_band(id, layer.astype(meta['dtype']))
        return self

class ColorCorrection():
    def __init__(self):
        pass
        
    def enhance(self, im):
    	im = exposure.equalize_adapthist(im, clip_limit=0.03)
    	im = exposure.rescale_intensity(im)
    	return im
    	 
    def color_balancing(self, infile, outfile, dtype = 'float32', ndv=0):
    	with rasterio.open(infile) as src:
	     im = src.read()
	     profile = src.profile
    	im = reshape_as_raster(self.enhance(reshape_as_image(im)))
    	profile.update(dtype=dtype, compress='lzw', nodata=ndv)
    	with rasterio.open(outfile, "w", **profile) as dst:
    		for i, layer in enumerate(im, start=1):  	
  	     		dst.write_band(i, layer.astype(dtype))
        
        
class Analysis():
    def __init__(self):
        pass
     
    @staticmethod    
    def get_mosaic(file_list, mosaic_file, srcnodata=0, dstnodata=0, method='maxband'):
    	fils = '-i '+ " -i ".join(file_list)
    	cmd = 'pkcomposite ' + fils + ' -o ' + mosaic_file + \
    		                  ' -srcnodata ' + str(srcnodata) + \
		                  ' -dstnodata ' + str(dstnodata) + \
		                  ' -cr ' + method
    	if os.path.exists(mosaic_file):
    		os.remove(mosaic_file)    	
    	subprocess.check_call(cmd, shell=True)
    
    @staticmethod
    def get_clipped(mosaic_file, crop_file, shapefile):
    	cmd = 'gdalwarp -cutline ' + shapefile + ' -crop_to_cutline -dstalpha ' + \
    	                             mosaic_file + ' ' + crop_file
    	if os.path.exists(crop_file):
    		os.remove(crop_file)    	
    	os.system(cmd)     
    	
class DNConversion():
    def __init__(self, base):
        self.files = glob.glob(base + '/*.TIF')
        self.basename = self.files[0].split('.TIF')[0].split('_T1_')[0]
        self.meta_file = glob.glob(base + '/*MTL*')[0]
        self.meta = grab_meta(self.meta_file)       
        
    def to_radiance(self, ofile, dtype = 'float32', bands = np.arange(1, 12)):
        corrected = []
        for id in bands:
            filename = self.basename + '_T1_' + 'B' + str(id) + '.TIF'
            with rasterio.open(filename) as src:
                band = src.read()
                profile = src.profile 
            M = getattr(self.meta, "RADIANCE_MULT_BAND_{0}".format(id))
            A = getattr(self.meta, "RADIANCE_ADD_BAND_{0}".format(id))
            rad = band * M + A
            corrected.append(np.squeeze(rad))
            
        profile.update(dtype=dtype, count=len(bands), compress='lzw')    
        with rasterio.open(ofile, "w", **profile) as dst:
            for i, layer in enumerate(corrected, start=1):
                dst.write_band(i, layer.astype(profile['dtype']))
        return self
    
    def to_reflectance(self, ofile, dtype = 'float32', bands = np.arange(1, 10)):
        corrected = []  
        for id in bands:
            filename = self.basename + '_T1_' + 'B' + str(id) + '.TIF'
            with rasterio.open(filename) as src:
                band = src.read()
                profile = src.profile
            M = getattr(self.meta, "REFLECTANCE_MULT_BAND_{0}".format(id))
            A = getattr(self.meta, "REFLECTANCE_ADD_BAND_{0}".format(id))
            S = getattr(self.meta, "SUN_ELEVATION")*(np.pi/180)
            rad = band * M + A
            rad = rad/np.sin(S)
            corrected.append(np.squeeze(rad))
            
        profile.update(dtype=dtype, count=len(bands), compress='lzw')    
        with rasterio.open(ofile, "w", **profile) as dst:
            for i, layer in enumerate(corrected, start=1):
                dst.write_band(i, layer.astype(profile['dtype']))
        return self
    
    @staticmethod
    def get_composite(ofile='fcc.tif', stackfile = 'ex.tif', \
                      band = [5,4,3]):
        dtype='float32'
        dataset = rasterio.open(stackfile)
        meta = dataset.meta
        meta.update(count = len(band), compress='lzw', dtype = dtype)
        corr = dataset.read(band)    
        if corr.max()<0:   
            corr *= -1
        corr = np.array([norm(layer) for layer in corr])
        
        # contrast stretching between 2 and 98 percentiles
        low, high = np.percentile(corr, (2, 98))
        corr = exposure.rescale_intensity(corr, in_range=(low, high))
        with rasterio.open(ofile, "w", **meta) as dest:
                dest.write(corr.astype(dtype)) 
