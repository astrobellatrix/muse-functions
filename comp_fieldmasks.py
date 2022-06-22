#! /usr/bin/env python
#
# FILE:    comp_fieldmasks.py
# AUTHORS: Lutz Wisotzki, Tanya Urrutia
# DESCR.:  Compute a sky mask image for the datacube
#

import argparse
import math as m
import numpy as np
from astropy.io import fits
import numpy.ma as ma
from scipy.ndimage import filters
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)

parser = argparse.ArgumentParser(description="""
Compute a sky mask image for the datacube using erosion and dilation.""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i","--input",
                    required=True,
                    type=str,
                    help="Name of the input FITS datacube for which the various calculations will be done. The datacube must have a flux and variance extension and optionally also an exposure cube extension.")
parser.add_argument("-o","--output",
                    type=str,
                    default="blankskymask.fits",
                    help="Name of the output blank sky mask FITS file.")
parser.add_argument("-S","--SHDU",
                    type=int,
                    default=1,
                    help="HDU number (0-indexed) in the input FITS file containing the flux data.")
parser.add_argument("-N","--NHDU",
                    type=int,
                    default=2,
                    help="HDU number (0-indexed) in the input FITS file containing the variance data.")
parser.add_argument("-E","--EHDU",
                    type=int,
                    default=3,
                    help="HDU number (0-indexed) in the exposure or input FITS file containing the exposure cube.")
parser.add_argument("-W","--WHDU",
                    type=int,
                    default=4,
                    help="HDU number (0-indexed) in the input FITS file containing the whitelight image.")
parser.add_argument("--expcube",
                    type=str,
                    default='None',
                    help="Name of the exposure FITS datacube. The default is that it is contained in the input datacube as extension 3.")
parser.add_argument("--whiteima",
                    type=str,
                    default='None',
                    help="Name of the whitelight FITS image. The default is that it is contained in the input datacube as extension 4.")
parser.add_argument("--medexpfrac",
                    type=float,
                    default=0.4,
                    help="Fraction of median exposure for field-of-view threshold.")
parser.add_argument("--bbl1",
                    type=int,
                    default=200,
                    help="Lower limit layer index of broad band continuum window.")
parser.add_argument("--bbl2",
                    type=int,
                    default=1800,
                    help="Upper limit layer index of broad band continuum window.")
parser.add_argument("--kappa_blanksky",
                    type=float,
                    default=3.0,
                    help="Kappa-sigma clipping threshold for blank sky mask.")
parser.add_argument("--tmpmasks",
                    action='store_true',
                    help="""Save three intermediate masks, (i) a FoV mask (fov_[base_name].fits), (ii) the broadband image (bb_ima_[base_name].fits) and (iii) the initial field mask before erosion and dilation (fmask_[base_name].fits).
                    """)
parser.add_argument("--base_name",
                    type=str,
                    default=None,
                    help="""Base name for the various temporary masks. The default expects the input cube to have the name 'DATACUBE_[base_name]_v2.0.fits'
                    """)


args = parser.parse_args()
in_datacube = args.input

if (args.expcube == 'None'):
    in_expcube = in_datacube
else:
    in_expcube = args.expcube

if (args.whiteima == 'None'):
    in_white = in_datacube
else:
    in_white = args.whiteima

if args.base_name == None:
    base_name = in_datacube[9:-10]
else:
    base_name = args.base_name

ihdu_d = args.SHDU
ihdu_v = args.NHDU
ihdu_e = args.EHDU
ihdu_w = args.WHDU
outputname = args.output

medexpfrac = args.medexpfrac
bbl1 = args.bbl1
bbl2 = args.bbl2
kappa_blanksky = args.kappa_blanksky

save_tmp = args.tmpmasks

# Open cube HDUs
cubeHDU = fits.open(in_datacube)
dhead = cubeHDU[ihdu_d].header
npix = ( dhead['naxis1'], dhead['naxis2'], dhead['naxis3'] )

ecubeHDU = fits.open(in_expcube)
e_cube = ecubeHDU[ihdu_e].data

wimaHDU = fits.open(in_white)
whead = wimaHDU[ihdu_w].header

# Create mean exposure image as average through all layers
exp_ima = np.nanmean(e_cube, axis=0)
del e_cube
ecubeHDU.close()

print("Opening cubes. Creating Field-of-View mask image")
# Create Field-of-View mask image, defined as all pixels in mean exposure 
# image above some minimum number of contributing exposures
# (1 = inside FoV, 0 = outside FoV)
medianexp = np.nanmedian(exp_ima[exp_ima > 0.])
threshold = medianexp * medexpfrac
fmask_fov = (exp_ima > threshold).astype(np.int16)
if save_tmp:
    fits.writeto('fov_%s.fits' % base_name,data=fmask_fov,header=whead,overwrite=True)
  
print("Creating a broad-band image of the cube (not full white-light)")
# Compute a broad-band image of the cube (rather than full white-light), 
# with zero-order background correction
d_cube = cubeHDU[ihdu_d].data
newdata = np.zeros(d_cube.shape,dtype=np.float32)
bb_ima = np.nanmean(d_cube[bbl1:bbl2+1,:,:], axis=0)
del d_cube
bb_ima[fmask_fov < 1] = np.NaN
bgrcor = np.nanmedian(bb_ima)
bb_ima = bb_ima - bgrcor
if save_tmp:
    fits.writeto('bb_ima_%s.fits' % base_name,data=bb_ima,header=whead,overwrite=True)

print("Masking all pixels above (kappa * bgrnoise)")
# Mask all pixels in broadband image above threshold = kappa * bgrnoise, 
# where bgrnoise is estimated from variance cube
v_cube = cubeHDU[ihdu_v].data
# 1D array with median variances per layer
vmed_vec = np.nanmedian(v_cube[bbl1:bbl2+1,:,:], axis=(1,2))  
del v_cube
# approximate factor correcting for covariance losses due to resampling 
corrfac = 1.7  
# bgrnoise in bb image = sqrt of mean variance 
sigbb = m.sqrt(np.nanmean(vmed_vec)/npix[2]) * corrfac      
threshold = kappa_blanksky * sigbb
tmpmask = (np.logical_and(bb_ima<threshold,fmask_fov>0)).astype(float)
if save_tmp:
    fits.writeto('fmask_%s.fits' % base_name,data=tmpmask,header=whead,overwrite=True)

print("Creating the final blank-sky mask with erosion and dilation operations.")
# Create mask of "blank sky" pixels (1 = blank sky within FoV, 
# 0 = everything else)
# by modifying pixels of initial mask in sequence of mask shrinking and 
# growing operations
# shrink mask to eliminate solitary pixels
np.place(tmpmask, filters.uniform_filter(tmpmask, size=3)*9>7.5, 1.)     
# grow mask unconditionally
np.place(tmpmask, filters.uniform_filter(tmpmask, size=3)*9<8.5, 0.)     
# grow mask multiple times with round edges
for i in range(0,5):
    np.place(tmpmask, filters.uniform_filter(tmpmask, size=3)*9<6.5, 0.)    
fmask_blanksky = tmpmask.astype(np.int16) 
# ZAP wants exactly the inverse kind of mask, so I could program it here, but decided against it
#fmask_blanksky_zap = np.where((fmask_blanksky==0)|(fmask_blanksky==1), fmask_blanksky^1, fmask_blanksky)

print("Writing...")
fits.writeto(outputname, data=fmask_blanksky, header=cubeHDU[ihdu_d].header, overwrite=True)
#fits.writeto(outputname[:-5] + '.zap.fits', data=fmask_blanksky_zap, header=cubeHDU[ihdu_d].header, overwrite=True)
del fmask_blanksky
    
print("Done.")



