#! /usr/bin/env python
#
# FILE:    comp_fieldmasks.py
# AUTHORS: Lutz Wisotzki, Tanya Urrutia
# DESCR.:  Compute various mask images for the datacube
#

import argparse
import math as m
import numpy as np
from astropy.io import fits
import numpy.ma as ma
from scipy.ndimage import filters

import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)

parser = argparse.ArgumentParser(description="""
Compute mask images for the datacube. 
Blurb about erosion and dilation.""",
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
                    help="HDU number (0-indexed) or name in the input FITS file containing the flux data.")
parser.add_argument("-N","--NHDU",
                    type=int,
                    default=2,
                    help="HDU number (0-indexed) or name in the input FITS file containing the variance data.")
parser.add_argument("-E","--EHDU",
                    type=int,
                    default=3,
                    help="HDU number (0-indexed) or name in the exposuer or input FITS file containing the exposure cube.")
parser.add_argument("--expcube",
                    type=str,
                    default='default',
                    help="Name of the exposure FITS datacube. The default is that it is contained in the input datacube as extension 3.")
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

args = parser.parse_args()
in_datacube = args.input
exp_in_input = args.expcube

if (exp_in_input == 'default'):
  in_expcube = in_datacube
else:
  in_expcube = exp_in_input

ihdu_d = args.SHDU
ihdu_v = args.NHDU
ihdu_e = args.EHDU
outputname = args.output

medexpfrac = args.medexpfrac
bbl1 = args.bbl1
bbl2 = args.bbl2
kappa_blanksky = args.kappa_blanksky


#def comp_fieldmasks( path, in_datacube, ihdu_d, ihdu_v, in_expcube, ihdu_e, parameters, out_rootname ):
#
#    print("Compute various mask images for cube %s" % in_datacube )
#    print("(please ignore 'empty slice' and 'invalid value' runtime warnings)." )

# Open cube HDUs
cubeHDU = fits.open(in_datacube)
npix = ( cubeHDU[ihdu_d].header['naxis1'],  cubeHDU[ihdu_d].header['naxis2'],  cubeHDU[ihdu_d].header['naxis3'] )
ecubeHDU = fits.open(in_expcube)

# Create mean exposure image as average through all layers
e_cube = ecubeHDU[ihdu_e].data
exp_ima = np.nanmean(e_cube, axis=0)
del e_cube
ecubeHDU.close()
#fits.writeto('expima.fits', exp_ima, header=ecubeHDU.header, overwrite=True)
#fits.writeto(path + '/' + out_rootname + '_expima.fits', exp_ima, header=ecubeHDU.header, overwrite=True)
#    print(" ... writing output mean exposure image %s" % out_rootname + '_expima.fits' )

# Create Field-of-View mask image, defined as all pixels in mean exposure 
# image above some minimum number of contributing exposures
# (1 = inside FoV, 0 = outside FoV)
medianexp = np.nanmedian(exp_ima[exp_ima > 0.])
threshold = medianexp * medexpfrac
fmask_fov = (exp_ima > threshold).astype(np.int16)

#    fits.writeto('fovmask.fits', fmask_fov, header=dcubeHDU.header, overwrite=True)
#    fits.writeto(path + '/' + out_rootname + '_fovmask.fits', fmask_fov, header=dcubeHDU.header, overwrite=True)
#    print(" ... writing output FoV mask %s for expima threshold = %6.3f" % (out_rootname + '_fovmask.fits', threshold) )

    # Compute a broad-band image of the cube (rather than full white-light), with zero-order background correction
d_cube = cubeHDU[ihdu_d].data
newdata = np.zeros(d_cube.shape,dtype=np.float32)
bb_ima = np.nanmean(d_cube[bbl1:bbl2+1,:,:], axis=0)
del d_cube
bb_ima[fmask_fov < 1] = np.NaN
bgrcor = np.nanmedian(bb_ima)
bb_ima = bb_ima - bgrcor
    #fits.writeto('bbima.fits', bb_ima, header=dcubeHDU.header, overwrite=True)
    #fits.writeto(path + '/' + out_rootname + '_bbima.fits', bb_ima, header=dcubeHDU.header, overwrite=True)
    #print(" ... writing broadband image %s for layer range %d - %d" % (out_rootname + '_bbima.fits', bbl1, bbl2) )

    # Mask all pixels in broadband image above threshold = kappa * bgrnoise, where bgrnoise is estimated from variance cube
v_cube = cubeHDU[ihdu_v].data
vmed_vec = np.nanmedian(v_cube[bbl1:bbl2+1,:,:], axis=(1,2))    
# 1D array with median variances per layer
del v_cube
corrfac = 1.7  
# approximate factor correcting for covariance losses due to resampling 
sigbb = m.sqrt(np.mean(vmed_vec)/npix[2]) * corrfac      
# bgrnoise in bb image = sqrt of mean variance 
threshold = kappa_blanksky * sigbb
tmpmask = (np.logical_and(bb_ima<threshold,fmask_fov>0)).astype(np.float)

#    fits.writeto('bbmask1.fits', tmpmask, header=dcubeHDU.header, overwrite=True)

    # Create mask of "blank sky" pixels (1 = blank sky within FoV, 0 = everything else)
    # by modifying pixels of initial mask in sequence of mask shrinking and growing operations
np.place(tmpmask, filters.uniform_filter(tmpmask, size=3)*9>7.5, 1.)     
# shrink mask to eliminate solitary pixels
#    fits.writeto('bbmask2.fits', tmpmask, header=dcubeHDU.header, overwrite=True)
np.place(tmpmask, filters.uniform_filter(tmpmask, size=3)*9<8.5, 0.)     
# grow mask unconditionally
#    fits.writeto('bbmask3.fits', tmpmask, header=dcubeHDU.header, overwrite=True)
for i in range(0,5):
    np.place(tmpmask, filters.uniform_filter(tmpmask, size=3)*9<6.5, 0.)    
# grow mask multiple times with round edges
fmask_blanksky = tmpmask.astype(np.int16) 
fits.writeto(outputname, data=fmask_blanksky, header=cubeHDU[ihdu_d].header, overwrite=True)
#    fits.writeto(path + '/' + out_rootname + '_blankskymask.fits', fmask_blanksky, header=dcubeHDU.header, overwrite=True)
#    print(" ... writing 'blank sky' mask %s " % (out_rootname + '_blankskymask.fits') )
del fmask_blanksky
    
print("Done.")
print()
#return



