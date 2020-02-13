import math as m
import numpy as np
from astropy.io import fits
import numpy.ma as ma
from scipy.ndimage import filters
from scipy.ndimage import morphology


def create_mfs_and_effvar_cube(path, in_datacube, ihdu_d, in_mfcube, ihdu_m, in_bgrstattab, in_expcube, ihdu_e, outcube ):

    print("Create a combined MFS and EFFVAR datacube: " )
    print(" -- subtract spectrally median filtered cube from original data," )
    print(" -- subtract DC background correction in each layer," )
    print(" -- compute effective variance from effective noise and exposure cube," )
    print(" -- set data to zero outside of field of view." )

    # Read original datacube
    dcubeHDU = fits.open(path + '/' + in_datacube)[ihdu_d]
    npix = ( dcubeHDU.header['naxis1'],  dcubeHDU.header['naxis2'],  dcubeHDU.header['naxis3'] )
    d_cube = dcubeHDU.data

    # Read spectrally median filtered cube
    mcubeHDU = fits.open(path + '/' + in_mfcube)[ihdu_m]
    mfs_cube = d_cube - mcubeHDU.data
    del d_cube
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
    print()
    return


##########################
if __name__ == "__main__":
    import sys

    # create_mfs_and_effvar_cube( path, in_datacube, ihdu_d, in_mfcube, ihdu_m, in_bgrstattab, in_expcube, ihdu_e outcube )
#    create_mfs_and_effvar_cube( 'candels-cdfs-47', 'DATACUBE_candels-cdfs-47_v1.0.fits', 1, 'median_filtered_DATACUBE_candels-cdfs-47_v1.0.fits', 2, 'cdfs-47_bgrstat.fits', 'DATACUBE_candels-cdfs-47_v1.0.fits', 3, 'cdfs-47_mfs-and-effvar-cube.fits' )
#    create_mfs_and_effvar_cube( '../candels-cosmos-17', 'DATACUBE_candels-cosmos-17_v1.0.fits', 1, 'median_filtered_DATACUBE_candels-cosmos-17_v1.0.fits', 2, 'cosmos-17_bgrstat.fits', 'DATACUBE_candels-cosmos-17_v1.0.fits', 3, 'cosmos-17_mfs-and-effvar-cube.fits' )
#    create_mfs_and_effvar_cube( '../udf-01', 'DATACUBE_udf-01_v1.0.fits', 1, 'median_filtered_DATACUBE_udf-01_v1.0.fits', 2, 'udf-01_bgrstat.fits', 'DATACUBE_udf-01_v1.0.fits', 3, 'udf-01_mfs-and-effvar-cube.fits' )

create_mfs_and_effvar_cube( '.', 'DATACUBE_QUBES_1422.fits', 1, 'median_filtered_DATACUBE_QUBES_1422.fits', 2, 'QUBES_1422_bgrstat.fits', 'DATACUBE_QUBES_1422_exp.fits', 0, 'QUBES_1422_mfs-and-effvar-cube.fits' )

#    create_mfs_and_effvar_cube ( sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), sys.argv[6], sys.argv[7], int(sys.argv[8]), sys.argv[9] )




