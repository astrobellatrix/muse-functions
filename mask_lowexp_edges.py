#!/usr/bin/env python
#
# FILE:    mask_lowexp_edges.py
# AUTHORS: Tanya Urrutia, Lutz Wisotzki
# DESCR.:  Compute and apply a low exposure edge mask to a cube
#

import argparse
import sys,time
import astropy.io.fits as fits
import numpy as np
from scipy.ndimage import filters

import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)

parser = argparse.ArgumentParser(description="""
Compute and apply a low exposure edge mask to a cube, the edge mask is dilated twice to thicken the edges that are masked.""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i","--input",
                    required=True,
                    type=str,
                    help="Base name of cube that is modified. Example: 'DATACUBE_SF_M0416-SF2_03.fits', then the base name is 'SF_M0416-SF_03'. It is expected that the datacube is zero-indexed and has flux in the first extension, the variance in the second extension and the whitelight image in the third extension (the default of the pipeline).")
parser.add_argument("-o","--output",
                    type=str,
                    default=None,
                    help="""
                    Name of the output masked FITS file. The output FITS file will contain the exact same HDUs with the corresponding voxels masked. [Default: `DATACUBE_[base_name].mask.fits`, where base name is the input base name used.]
                    """)
parser.add_argument("-m","--outmask",
                    type=str,
                    default=None,
                    help="""
                    Name of the final edgemask FITS file, saved as 8-bit integer file, where 1 is the pixel masked and 0 remains untouched. [Default: `edgemask_[base_name].fits`, where base name is the input base name used.]
                    """)
parser.add_argument("--nnan_wave_thresh",
                    type=int,
                    default=3000,
                    help="""
                    Threshold wavelength layers below which the pixels will get masked, this corresponds to about 83 percent for AO and 81 percent for non-AO data.
                    """)
parser.add_argument("--tmpmask",
                    action='store_true',
                    help="""Save the various intermediate masks with the default as mask[n]_[base_name].fits.
                    """)

args = parser.parse_args()

name = args.input
nnan_wave_thresh = args.nnan_wave_thresh

if args.output == None:
    out_fits = 'DATACUBE_' + name + '.mask.fits'
else:
    out_fits = args.output

if args.outmask == None:
    out_mask = 'edgemask_' + name + '.fits'
else:
    out_mask = args.outmask

save_tmp = args.tmpmask

starttime = time.time()

################

print('Generating and applying a low exposure mask for %s' % name)
print('...')

muse_hdu = fits.open('DATACUBE_%s.fits' % name)
mdata = muse_hdu[1].data
sdata = muse_hdu[2].data
wdata = muse_hdu[3].data
mask = np.ones(wdata.shape,dtype=bool)
newdata = np.zeros(mdata.shape,dtype=np.float32)
newstat = np.zeros(sdata.shape,dtype=np.float32)
newwhite = np.zeros(wdata.shape,dtype=np.float32)

print('Create a masking threshold based on a collapsed weight cube,')
print('in which NaN values receive a 0 and everything else a 1.')   
passtime='('+str(round(time.time()-starttime,3))+'s)' 
print('... '+passtime)

# create a masking threshold based on the weight cube
weight_cube = np.ones(mdata.shape,dtype=np.int64)
weight_cube[np.isnan(mdata)] = 0
white_weight = np.sum(weight_cube,axis=0)

print('Creating the mask.')
print('One erosion filtering to get rid of the inner low exposure pixels.')
print('Two dilation filterings to expand / thicken the edge mask.')
passtime='('+str(round(time.time()-starttime,3))+'s)'
print('... '+passtime)
   
mask[white_weight > nnan_wave_thresh] = False
mask = mask.astype(float)
if save_tmp:
    mosko1 = mask.astype(bool)
    mosko1[white_weight==0] = False
    fits.writeto('mask1_%s.fits' % name,data=mosko1.astype(np.uint8),header=muse_hdu[3].header,overwrite=True)

# use erosion to get rid of individual pixels inside
np.place(mask, filters.uniform_filter(mask, size=4)*16<4.5, 0.)
if save_tmp:
    mosko2 = mask.astype(bool)
    mosko2[white_weight==0] = False
    fits.writeto('mask2_%s.fits' % name,data=mosko2.astype(np.uint8),header=muse_hdu[3].header,overwrite=True)

# two dilations to expand the edge mask
np.place(mask, filters.uniform_filter(mask, size=3)*9>1.5, 1.)
if save_tmp:
    mosko3 = mask.astype(bool)
    mosko3[white_weight==0] = False
    fits.writeto('mask3_%s.fits' % name,data=mosko3.astype(np.uint8),header=muse_hdu[3].header,overwrite=True)
np.place(mask, filters.uniform_filter(mask, size=3)*9>1.5, 1.)
if save_tmp:
    mosko4 = mask.astype(bool)
    mosko4[white_weight==0] = False
    fits.writeto('mask4_%s.fits' % name,data=mosko4.astype(np.uint8),header=muse_hdu[3].header,overwrite=True)
    
# save the mask
mask_hdu = fits.PrimaryHDU(data=mask.astype(np.uint8),header=muse_hdu[3].header)
mask_hdu.writeto(out_mask,overwrite=True)
# convert to a boolean mask
mask = mask.astype(bool)

print('Now mask the corresponding pixels in the three extensions')
print('and write it out.')
passtime='('+str(round(time.time()-starttime,3))+'s)'
print('... '+passtime)

# do the actual edge masking
for wave in range(mdata.shape[0]):
    imdata = mdata[wave]
    imstat = sdata[wave]                 
    imdata[mask] = np.nan
    imstat[mask] = np.nan
    newdata[wave] = imdata
    newstat[wave] = imstat
imwhite = wdata
imwhite[mask] = np.nan
newwhite =  imwhite

# write out the masked data
prihdu = fits.PrimaryHDU(header=muse_hdu[0].header)
datan = fits.ImageHDU(data=newdata,header=muse_hdu[1].header)
statn = fits.ImageHDU(data=newstat,header=muse_hdu[2].header)
whiten = fits.ImageHDU(data=newwhite,header=muse_hdu[3].header)
hdulist = fits.HDUList([prihdu,datan,statn,whiten])
hdulist.writeto(out_fits,overwrite=True)
hdulist.close()

passtime='('+str(round(time.time()-starttime,3))+'s)'
print('Done. '+passtime)
