#! /usr/bin/env python
#
# FILE:    comp_bgrstat.py
# AUTHORS: Lutz Wisotzki, Tanya Urrutia
# DESCR.:  Compute background statistics for the cube
#

import argparse
import math as m
import numpy as np
from astropy.io import fits
from astropy.table import Table
import numpy.ma as ma
from scipy.ndimage import filters
from scipy import interpolate

parser = argparse.ArgumentParser(description="""
Compute background statistics for the datacube. 
""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i","--input",
                    required=True,
                    type=str,
                    help="Name of the input FITS datacube for which the various calculations will be done. The datacube must have a flux and variance extension and optionally also an exposure cube extension.")
parser.add_argument("-m","--skymask",
                    required=True,
                    type=str,
                    help="Name of the input sky mask computed from a previous step, most likely with comp_fieldmask.py. 0 = masked, 1 = unmasked.")
parser.add_argument("-f","--medfilt",
                    required=True,
                    type=str,
                    help="Name of the median filtered datacube. This is most likely computed with LSDCat's median-filter-cube.py routine.")
parser.add_argument("-o","--output",
                    type=str,
                    default="bgrstat.fits",
                    help="Name of the output background stat FITS table.")
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
parser.add_argument("--statmethod",
                    type=str,
                    default="fits",
                    help="Way of calculating the wavelength dependent resampling corrections deduced from a random cube. The options are: 'fits' in which the rancubestat.fits file available on the repository must be present in the directory and 'poly' in which a 2nd order polynomial was fitted through the rancubestat and can extend to lower wavelengths, e.g. for AO data.")
parser.add_argument("--poly0",
                    type=float,
                    default=5.57934039e-01,
                    help="0th order of polynomial if statmethod='poly'.")
parser.add_argument("--poly1",
                    type=float,
                    default=-2.11008327e-06,
                    help="1st order of polynomial if statmethod='poly'.")
parser.add_argument("--poly2",
                    type=float,
                    default=2.24603677e-10,
                    help="2nd order of polynomial if statmethod='poly'.")
parser.add_argument("--dc_window_size",
                    type=int,
                    default=301,
                    help="Size of the median filter window for the DC-subtraction.")

args = parser.parse_args()
in_datacube = args.input
ihdu_d = args.SHDU
ihdu_v = args.NHDU
mfhdu_mf = args.MFHDU
outputname = args.output
mfcube = args.medfilt
inmask = args.skymask
statmethod = args.statmethod
poly0 = args.poly0
poly1 = args.poly1
poly2 = args.poly2
dc_window_size = args.dc_window_size

if (statmethod == 'fits'):
    in_rancubestattab = 'rancubestat.fits'
if (statmethod == 'poly'):
    # gathered from fitting rancubestat.fits above
    fitpoly = np.poly1d([ poly2, poly1, poly0 ])

print("Compute background statistics for cube %s" % in_datacube )

print("Opening datacubes...")
cubeHDU = fits.open(in_datacube)
dhead = cubeHDU[ihdu_d].header
d_cube = cubeHDU[ihdu_d].data
npix = ( dhead['naxis1'], dhead['naxis2'], dhead['naxis3'] )
start = ( dhead['crval1'], dhead['crval2'], dhead['crval3'] )
step = ( dhead['cd1_1'], dhead['cd2_2'], dhead['cd3_3'] )
v_cube = cubeHDU[ihdu_v].data

mfcubeHDU = fits.open(mfcube)
mf_cube = mfcubeHDU[mfhdu_mf].data
mfs_cube = d_cube - mf_cube
del mf_cube

# Input mask image
bs_mask = fits.getdata(inmask)
np_unmasked = bs_mask[bs_mask > 0].size

print("Initializing arrays...")
# Initialise 1D arrays, 1 element per cube layer
wlen = npix[2]
# wavelength
wave = np.zeros(wlen,dtype=np.float32)
# mean of unmasked data, no continuum subtraction
meandata = np.zeros(wlen,dtype=np.float32) 
# median of unmasked data, no continuum subtraction
mediandata = np.zeros(wlen,dtype=np.float32)
# standard deviation of unmasked data 
sdevdata = np.zeros(wlen,dtype=np.float32) 
# robust estimate of standard deviation from quartile distance in unmasked data
qdsigdata = np.zeros(wlen,dtype=np.float32)
# mean of unmasked data after subtraction of spectrally median filtered 
# cube (MFS)
meanmfs = np.zeros(wlen,dtype=np.float32)  
# median of unmasked data after subtraction of spectrally median filtered 
# cube (MFS)
medianmfs = np.zeros(wlen,dtype=np.float32) 
# standard deviation of unmasked MFS data 
sdevmfs = np.zeros(wlen,dtype=np.float32) 
# robust estimate of standard deviation from quartile distance in 
# unmasked MFS data    
qdsigmfs = np.zeros(wlen,dtype=np.float32)   
# median propagated variance of unmasked pixels
medpvar = np.zeros(wlen,dtype=np.float32)      

print("Calculating statistics per layer images...")
# Loop over all layers in cubes, calculate statistics per layer images
for l in range(0,npix[2]):
    wave[l] = start[2] + step[2] * l
    datalayer = d_cube[l,:,:]    # reduced cube data
    pvarlayer = v_cube[l,:,:]    # formally propagated variances
    mfslayer = mfs_cube[l,:,:]   # cube data after subtraction of spectrally median filtered cube
    data_masked = ma.masked_where(1-bs_mask, datalayer)
    pvar_masked = ma.masked_where(1-bs_mask, pvarlayer)
    mfs_masked = ma.masked_where(1-bs_mask, mfslayer)
    mediandata[l] = np.nanmedian(ma.compressed(data_masked))   
    meandata[l] = np.nanmean(ma.compressed(data_masked))
    sdevdata[l] = np.nanstd(ma.compressed(data_masked))
    q25 = np.nanpercentile(ma.compressed(data_masked),25.)
    q75 = np.nanpercentile(ma.compressed(data_masked),75.)
    # 0.7413 = ratio between standard deviation and quartile distance 
    # for a normal distribution 
    qdsigdata[l] = 0.7413 * (q75 - q25)     
    medianmfs[l] = np.nanmedian(ma.compressed(mfs_masked))
    meanmfs[l] = np.nanmean(ma.compressed(mfs_masked))
    sdevmfs[l] = np.nanstd(ma.compressed(mfs_masked))
    q25 = np.nanpercentile(ma.compressed(mfs_masked),25.)
    q75 = np.nanpercentile(ma.compressed(mfs_masked),75.)
    qdsigmfs[l] = 0.7413 * (q75 - q25) 
    medpvar[l] = np.nanmedian(ma.compressed(pvar_masked))
meandata[np.isnan(meandata)] = 0
meanmfs[np.isnan(meanmfs)] = 0
mediandata[np.isnan(mediandata)] = 0
medianmfs[np.isnan(medianmfs)] = 0

#  Below are ad-hoc modifications for AO, which we need to properly implement for MUSCATEL.
#  Let's postpone this to later. (but we should not forget them!) 
#    in1, in2 = 850, 1050
#    md_in1 = np.median(meandata[in1-25:in1])
#    md_in2 = np.median(meandata[in2:in2+25])
#    inc = (md_in2-md_in1)/(in2 - in1)
#    meandata[in1:in2] = md_in1 + inc * np.arange(in2-in1)
#    mm_in1 = np.median(meanmfs[in1-25:in1])
#    mm_in2 = np.median(meanmfs[in2:in2+25])
#    inc = (mm_in2-mm_in1)/(in2 - in1)
#    meanmfs[in1:in2] = mm_in1 + inc * np.arange(in2-in1)

print("Calculating DC background correction")
# Calculate DC background correction for original datacube 
# step 1: Produce a spectrally smoothed version of the mean values per layer 
medfil_fulwid = dc_window_size
gaufil_sigma = 20.
mediandata_mf = filters.median_filter(mediandata, size=medfil_fulwid, mode='nearest')
mediandata_gmf = filters.gaussian_filter1d(mediandata_mf, sigma=gaufil_sigma, mode='nearest')
# step 2: Per default we use this smoothed version for the DC correction
dc_cor = mediandata_gmf
# step 3: If the difference between layer mean and the smoothed version 
# deviates significantly from random then use the mean of that layer for 
# the DC correction
# NOTE: Jan 12 2022 - This step has been discontinued temporarily, since the 
# newest ZAP version doesn't overcorrect the sky lines as much anymore, 
# but I need to test this still (Tanya Urrutia)
# expected standard error of the mean = sdev within each 
# layer / sqrt(number of pixels)
#sigmean = qdsigdata/m.sqrt(np_unmasked)
# ad-hoc value for kappa-sigma clipping
#kappa = 4. 
# replace default DC correction values when indicated 
#np.copyto(dc_cor, meandata, where=np.fabs(meandata-meandata_gmf) > kappa*sigmean ) 

# No DC background correction anymore!!!
# Calculate DC background correction for MFS cube (same process as above)
#meanmfs_mf = filters.median_filter(meanmfs, size=medfil_fulwid, mode='nearest')
#meannmfs_gmf = filters.gaussian_filter1d(meanmfs_mf, sigma=gaufil_sigma, mode='nearest')
mdc_cor = np.zeros(wlen,dtype=np.float32)
# standard deviations of the layer means in MFS cube               
#sigmeanmfs = qdsigmfs/m.sqrt(np_unmasked)
#np.copyto(mdc_cor, meanmfs, where=np.fabs(meanmfs-meanmfs_gmf) > kappa*sigmeanmfs )

print("Calculating effective noise from MFS cube statistics...")
# Calculate effective noise from statistics of MFS cube
if (statmethod == 'fits'):
    # step 1: Open reference table with layer-by-layer correction factors for 
    # loss of covariances due to resampling
    ranstat_tab = Table.read(in_rancubestattab)
    ranstatwave = ranstat_tab['wave']
    ranstatcorrfac = ranstat_tab['srsig']
    # step 2: Define an interpolator for correction factors in wavelength
    # this provides a function
    corrfac_interpolate = interpolate.interp1d(ranstatwave, ranstatcorrfac, fill_value='extrapolate')  
    # step 3: Calculate default "effective noise" as layer-by-layer standard 
    # deviation, corrected for resampling
    effsigmfs = qdsigmfs/corrfac_interpolate(wave)
    # step 4: Calculate ratios of layer-by-layer standard deviations and 
    # formally propagated errors
    ratio_qdsig_pvar = qdsigmfs/np.sqrt(medpvar)
    # step 5: Rescale formally propagated errors by the mean of this ratio
    medpvar_rescaled = np.nanmean(ratio_qdsig_pvar) * np.sqrt(medpvar)/corrfac_interpolate(wave)
    # step 6: If the layer-by-layer standard deviation is below the 
    # formally propagated error then the effective noise gets replaced by 
    # the rescaled propagated error value
    np.copyto(effsigmfs, medpvar_rescaled, where=(ratio_qdsig_pvar < 1.))
if (statmethod == 'poly'):
    # Step 1 - 3 can quickly get the correction factor from the polynomial.
    # No need for interpolation
    # Steps 4 - 6 are are the same as above
    corrfac_poly = fitpoly(wave-7000)
    effsigmfs = qdsigmfs/corrfac_poly
    ratio_qdsig_pvar = qdsigmfs/np.sqrt(medpvar)
    medpvar_rescaled = np.nanmean(ratio_qdsig_pvar) * np.sqrt(medpvar)/corrfac_poly
    np.copyto(effsigmfs, medpvar_rescaled, where=(ratio_qdsig_pvar < 1.))


# Write output table
outtab = Table( [ wave, meandata, mediandata, sdevdata, qdsigdata, meanmfs, medianmfs, sdevmfs, qdsigmfs, medpvar, dc_cor, mdc_cor, ratio_qdsig_pvar, effsigmfs ],
    names=('wave', 'meandata', 'mediandata', 'sdevdata', 'qdsigdata', 'meanmfs', 'medianmfs', 'sdevmfs', 'qdsigmfs', 'medpvar', 'dc_cor', 'mdc_cor', 'r', 'effsig') )
outtab.write(outputname, overwrite=True)

print("Done.")

