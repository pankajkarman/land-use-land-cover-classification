
# coding: utf-8

# In[1]:

import warnings; warnings.filterwarnings('ignore')


# In[2]:


from lulc_utils import *
from skimage import exposure
import rasterio, numpy as np, glob, matplotlib.pyplot as plt
import pandas as pd, xarray as xr, scipy.misc as sm
import salem, geopandas as gpd
from facets import facets
from rasterio import plot 


# In[3]:


# cloud correction
band = [5, 4, 3]

# for fold in sorted(glob.glob('/mnt/ext1/data/sur/lulc/cloud/*T1*')):
#     tiles = sorted(glob.glob(fold + '/*tile*'))
#     for tile in tiles:
#         stack_file = tile + '/stack.tif'
#         fcc_file = tile + '/fcc.tif'
#         ofile = fcc_file.split('.tif')[0] + '.acorr.tif'
#         print ofile
#         conv = DNConversion(tile).to_reflectance(stack_file)
#         conv.get_composite(ofile = fcc_file, stackfile = stack_file, band = band)
#         atm = AtmosCorrection(fcc_file).correct(ofile, model = -1)


# In[4]:


# cloud correction
files_to_merge = []
fcc_list = []
base = '/mnt/ext1/data/sur/lulc/'
folds = glob.glob('/mnt/ext1/data/sur/lulc/cloud/*T1*')
for fold in folds:
    fcc_file = fold + '/tile1/fcc.tif'
    file1 = fold + '/tile1/fcc.acorr.tif'
    mask1 = fold + '/tile1/mask.tif'
    file2 = fold + '/tile2/fcc.acorr.tif'
    ofile = fold + '/fcc.both.tif'    
    files_to_merge.append(ofile)
    fcc_list.append(fcc_file)
    print ofile
    
#     try:
#         corr = CloudCorrection(file1, mask1, file2).correct(ofile)
#     except:
#         file3 = fold + '/tile2/fcc1.tif'
#         !rio warp --overwrite $file2 $file3 --like $file1
#         corr = CloudCorrection(file1, mask1, file3).correct(ofile)


# In[5]:


# No cloud correction
for fold in sorted(glob.glob('/mnt/ext1/data/sur/lulc/cloudless/*')):
    stack_file = fold + '/stack.tif'
    fcc_file = fold + '/fcc.tif'
    ofile = fcc_file.split('.tif')[0] + '.acorr.tif'
    print ofile
#     conv = DNConversion(fold).to_reflectance(stack_file)
#     conv.get_composite(ofile = fcc_file, stackfile = stack_file, band = band)
#     atm = AtmosCorrection(fcc_file).correct(ofile, model = -1)
    files_to_merge.append(ofile)
    fcc_list.append(fcc_file)


# In[6]:


final = []
for afile in files_to_merge:
    #print afile
    cfile = afile.split('.tif')[0] + '.color.tif'
    wgs = afile.split('.tif')[0] + '.color.wgs84.tif'
    print wgs
#     clr = ColorCorrection().color_balancing(afile, cfile)
#     wgs1 = reproject_dataset(cfile)
    final.append(wgs)


# In[7]:


# mosaicing of files
shapefile = '/home/pankaj/phd/practice/sur/data/neshpfile/NE_utm.shp'
base = '/mnt/ext1/data/sur/lulc/'
mosfile = base + 'mosaic.tif'
clip_file = base + 'clipped.tif'
print clip_file
print mosfile
anl = Analysis()
# anl.get_mosaic(final, mosfile)
anl.get_clipped(mosfile, clip_file, shapefile)

