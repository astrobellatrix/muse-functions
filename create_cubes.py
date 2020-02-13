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

def create_mfs_and_effvar_cube(path, in_datacube, ihdu_d, in_mfcube, ihdu_m, in_bgrstattab, in_expcube, ihdu_e, outcube ):

    print("Create a combined MFS and EFFVAR datacube: " )
    print(" -- subtract spectrally median filtered cube from original data," )
    print(" -- subtract DC background correction in each layer," )
    print(" -- compute effective variance from effective noise and exposure cube," )
    print(" -- set data to zero outside of field of view." )

in_datacube = args.input
mf_datacube = args.medfilt
ihdu_d = args.SHDU
ihdu_v = args.NHDU
mfhdu_mf = args.MFHDU

# Read original datacube
cubeHDU = fits.open(in_datacube)
dhead = cubeHDU[ihdu_d].header
d_cube = cubeHDU[ihdu_d].data
npix = ( dhead['naxis1'], dhead['naxis2'], dhead['naxis3'] )
del d_cube

# Read median filtered cube
mfcubeHDU = fits.open(mf_datacube)
mf_cube = mfcubeHDU[mfhdu_mf].data
mfs_cube = d_cube - mcubeHDU.data
del mcubeHDU.data

    # Read exposure cube
    ecubeHDU = fits.open(path + '/' + in_expcube)[ihdu_e]
    exp_cube = ecubeHDU.data

    # Read table data with DC correction and effective noise per layer
    intab = fits.open(path + '/' + in_bgrstattab)[1].data
    mdccor = intab['mdc_cor']
    effsig = intab['effsig']

    # Create new variance cube, equal to square of effective noise weighted by the relative exposure.
    # Everything outside the field-of-view is set to zero.
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
    effvar_cube[np.isinf(effvar_cube)] = 0
    effvar_cube[np.isnan(effvar_cube)] = 0
    effvar_cube[np.isnan(mfs_cube)] = 0
    mfs_cube[np.isnan(mfs_cube)] = 0
    mfs_cube[effvar_cube==0] = 0

    # Write to output file
    out_HDU_0 = fits.PrimaryHDU()
    out_HDU_1 = fits.ImageHDU(mfs_cube, header=dcubeHDU.header, name='DATA(MFS)')
    out_HDU_2 = fits.ImageHDU(effvar_cube, header=dcubeHDU.header, name='EFFVAR')
    OUT = fits.HDUList([out_HDU_0, out_HDU_1, out_HDU_2])
    OUT.writeto(path + '/' + outcube, output_verify='ignore', overwrite=True)
    print(" ... writing output cube %s" % outcube )

    print("Done.")
    





