### MUSE Analysis Functions developed at AIP ###

This is a handy list of functions used at the AIP for analysis of deep field data, such as MUSE-Wide and MUSCATEL in addition to the well known LSDCat, QtClassify and TDOSE.

## Requirements ##

Python 3.5+
* astropy (>=2.0) - http://www.astropy.org/
* NumPy - http://www.numpy.org/
* SciPy - http://www.scipy.org/

## Install ##

Simply clone and run the python functions. You might want to add the muse-functions folder to your $PATH and $PYTHONPATH

## Commands specific to running effective noise calculation ##

As detailed in 3.2.4 of the MUSE-Wide data release paper (Urrutia et al. 2019), the pipeline estimated variances are underestimated due to some of the power shifting to the covariances in the resampling process. We have developed a method to calculate the variance using empirically calibrated functions. The procedure involves 3 steps

* comp_fieldmasks.py - Creates a skymask

Creates a sky mask image for the datacube. It filters bright sources (with kappa-sigma thresholding) using erosion and dilation methods. It thereby ignores single bright pixels and ensures a proper object mask. This function is also used in other cases, e.g. to create a sky mask for superflat cubes. The sky mask created by this function is further used by comp_bgrstat.py. Note that this function required an exposure map, in MUSE-Wide data this is usually found in extension 3, which is the default.

* comp_bgrstat.py

Computes various background statistics for the datacube, in particular DC offset and effective noise. Requires a median filtered datacube, which can be obtained using LSDCat's median-filter-cube.py routine. The default also requires the supplied rancubestat.fits file. For detailed description on how the data are generated, please refer to Section 3.2.4 of the MUSE-Wide Data Release paper.

* create_cubes.py

Create DC-subtracted, effnoised DATA- and MFS-cubes.
-- subtract DC background correction in each layer,
-- compute effective variance from effective noise and exposure cube,
-- apply these to the datacubes and median-filtered datacubes

* Example run using MUSE-Wide DR2 data:

comp_fieldmasks.py -i DATACUBE_candels-cdfs-11_v2.0.fits -o candels-cdfs-11_v2.0_blankskymask.fits

comp_bgrstat.py -i DATACUBE_candels-cdfs-11_v2.0.fits -m candels-cdfs-11_v2.0_blankskymask.fits -f median_filtered_DATACUBE_candels-cdfs-11_v2.0.fits -o candels-cdfs-11_v2.0_bgrstat.fits

create_cubes.py -i DATACUBE_candels-cdfs-11_v2.0.fits -f median_filtered_DATACUBE_candels-cdfs-11_v2.0.fits -bg candels-cdfs-11_v2.0_bgrstat.fits -od DATACUBE_candels-cdfs-11_v2.0_dcsub_effnoised.fits -om DATACUBE_MFS_candels-cdfs-11_v2.0_dcsub_effnoised.fits

## Commands specific to running slice sky subtraction calculation ##

As described in the MUSE-Wide data release paper, for the first version of data reduction a previous step of sky subtraction was done before using ZAP. This worked on the basis that the LSF in a slice was self similar and also adjusted IFU sky background based on the strength of these sky lines. The routine is obsolete as its main goal of creating a zero-level sky background are now implemented in the autocal option in the pipeline. Hence, now the pipeline sky subtraction is the preferred option. The routines work on the individual IFU PIXTABLE basis. A FOV image in PIXTABLE (not sky) coordinates is required.

Note that these routines were written in 2014/5 and used Python 2.7 and a pyfits instead of astropy.

* fit_levels.py

Masks out the brightest pixels and creates 48 sky spectra with a 0.2A sampling (one per slice) using mean with a 5 sigma clip. Calculates the relative sky levels based on a gauss fit to three bright, isolated sky lines.

* subtract_av_slice.py

Subtracts the sky on the PIXTABLES on a slice per slice basis scaling by the relative flux deduced in the previous step by the fit to the sky lines.