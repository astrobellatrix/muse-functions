#! /usr/bin/env python
#
# FILE:    create_cubes.py
# AUTHORS: Lutz Wisotzki, Tanya Urrutia
# DESCR.:  Create DC-subtracted, effnoised DATA- and MFS-cubes
#

import argparse
import math as m
import numpy as np
from astropy.io import fits
import numpy.ma as ma
from scipy.ndimage import filters
from scipy.ndimage import morphology

parser = argparse.ArgumentParser(description="""
Create DC-subtracted, effnoised DATA- and MFS-cubes.
-- subtract spectrally median filtered cube from original data,
-- subtract DC background correction in each layer,
-- compute effective variance from effective noise and exposure cube.
""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i","--input",
                    required=True,
                    type=str,
                    help="Name of the input FITS datacube for which the various calculations will be done. The datacube must have a flux and variance extension and optionally also an exposure cube extension.")
parser.add_argument("-f","--medfilt",
                    required=True,
                    type=str,
                    help="Name of the median filtered datacube. This is most likely computed with LSDCat's median-filter-cube.py routine.")
parser.add_argument("-bg","--bgrstat",
                    required=True,
                    type=str,
                    help="Name of the background statistic FITS file created with comp_brgstat.py")
parser.add_argument("-od","--outdatacube",
                    required=True,
                    type=str,
                    default='DATACUBE_dcsub_effnoised.fits',
                    help="Name of the effnoised DC-subtracted FITS datacube")
parser.add_argument("-om","--outmfscube",
                    required=True,
                    type=str,
                    default='DATACUBE_MFS_dcsub_effnoised.fits',
                    help="Name of the effnoised DC-subtracted median-filterd subtracted FITS datacube")
parser.add_argument("--expcube",
                    type=str,
                    default='default',
                    help="Name of the exposure FITS datacube. The default is that it is contained in the input datacube as extension 3.")
parser.add_argument("-S","--SHDU",
                    type=int,
                    default=1,
                    help="HDU number (0-indexed) or name in the input FITS file containing the flux data.")
parser.add_argument("-N","--NHDU",
                    type=int,
                    default=2,
                    help="HDU number (0-indexed) or name in the input FITS file containing the variance data.")
parser.add_argument("-MF","--MFHDU",
                    type=int,
                    default=2,
                    help="HDU number (0-indexed) or name in the median-filtered FITS file containing the median-filtered data.")
parser.add_argument("-E","--EHDU",
                    type=int,
                    default=3,
                    help="HDU number (0-indexed) or name in the exposuer or input FITS file containing the exposure cube.")


args = parser.parse_args()
in_datacube = args.input
mf_datacube = args.medfilt
bgrstat = args.bgrstat
outdatacube = args.outdatacube
outmfscube = args.outmfscube
ihdu_d = args.SHDU
ihdu_v = args.NHDU
mfhdu_mf = args.MFHDU
ihdu_e = args.EHDU
exp_in_input = args.expcube

if (exp_in_input == 'default'):
  in_expcube = in_datacube
else:
  in_expcube = exp_in_input

print('Reading in data...')
# Read original datacube
cubeHDU = fits.open(in_datacube)
dhead = cubeHDU[ihdu_d].header
d_cube = cubeHDU[ihdu_d].data
npix = ( dhead['naxis1'], dhead['naxis2'], dhead['naxis3'] )

# Read median filtered cube
mfcubeHDU = fits.open(mf_datacube)
mf_cube = mfcubeHDU[mfhdu_mf].data
mfs_cube = d_cube - mf_cube
del mf_cube

# Read exposure cube
ecubeHDU = fits.open(in_expcube)
exp_cube = ecubeHDU[ihdu_e].data

# Read table data with DC correction and effective noise per layer
intab = fits.getdata(bgrstat)
mdccor = intab['mdc_cor']
dccor = intab['dc_cor']
effsig = intab['effsig']

print('Creating a new variance cube and subtracting DC correction')
# Create new variance cube, equal to square of effective noise weighted by 
# the relative exposure.
medexp = np.nanmedian(exp_cube[exp_cube>0.])
effvar_cube = np.zeros(exp_cube.shape, np.float32)
varlayer = np.empty((npix[1], npix[0]), np.float32)

errorhandle = np.seterr(all='ignore')
for l in range(npix[2]):
    explayer = exp_cube[l,:,:]
    varlayer.fill(effsig[l]**2)
    effvar_layer = np.divide(varlayer*medexp,explayer)
    effvar_cube[l,:,:] = effvar_layer
    mfs_cube[l,:,:] -= mdccor[l]
    d_cube[l,:,:] -= dccor[l]
effvar_cube[np.isinf(effvar_cube)] = np.nan

# Create a new whitelight image from DC-subtracted data
ma_d_cube = np.ma.masked_array(d_cube,np.isnan(d_cube))
mean = np.mean(ma_d_cube,axis=0)
dc_white = mean.filled(np.nan)

# Convert to 32-bit to save space
effvar_cube = np.float32(effvar_cube)
mfs_cube = np.float32(mfs_cube)
d_cube = np.float32(d_cube)
dc_white = np.float32(dc_white)

print(" ... writing output cubes %s and %s" % (outdatacube, outmfscube))
# Write to output file
# MFS-cube
out_MFSHDU_0 = fits.PrimaryHDU()
out_MFSHDU_1 = fits.ImageHDU(mfs_cube, header=mfcubeHDU[mfhdu_mf].header, name='DATA(MFS)')
out_MFSHDU_2 = fits.ImageHDU(effvar_cube, header=cubeHDU[ihdu_v].header, name='EFFVAR')
OUTMFS = fits.HDUList([out_MFSHDU_0, out_MFSHDU_1, out_MFSHDU_2])
OUTMFS.writeto(outmfscube, output_verify='ignore', overwrite=True)
# DATA-cube
out_dHDU_0 = fits.PrimaryHDU(header=cubeHDU[0].header)
out_dHDU_1 = fits.ImageHDU(d_cube, header=cubeHDU[ihdu_d].header, name='DATA (DCcor)')
out_dHDU_2 = fits.ImageHDU(effvar_cube, header=cubeHDU[ihdu_d].header, name='EFFVAR')
out_dHDU_3 = fits.ImageHDU(exp_cube, header=ecubeHDU[ihdu_e].header)
out_dHDU_4 = fits.ImageHDU(dc_white, header=cubeHDU[4].header)
OUTD = fits.HDUList([out_dHDU_0,out_dHDU_1,out_dHDU_2,out_dHDU_3,out_dHDU_4])
OUTD.writeto(outdatacube, output_verify='ignore', overwrite=True)

print("Done.")
