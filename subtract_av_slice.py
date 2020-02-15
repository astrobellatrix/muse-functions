#!/usr/bin/python

import sys, shutil
import astropy.io.fits as pyfits
import numpy as np
import scipy.interpolate as interp
import datetime

def main():
  # define base names to read the pixel tables, set up...
  print 'Setting up, copying...'
  base_name = sys.argv[1:]
  if not base_name:
    print 'Please provide a base name, e.g. PIXTABLE_TNJ1338_12'
    sys.exit(1)
  base_name = str(base_name[0])
  relflux_filename = base_name + '_relflux.txt'
  rel_fluxes = open(relflux_filename,'r')
  relflux_table = []
  for line in rel_fluxes:
    line = line.replace('[','')
    line = line.replace(']','')
    line = line.rstrip('\n')
    lspl = line.split()
    for each in lspl:
      relflux_table.append(each)

  for nifu in range(1,25):
    name = base_name + '_%02d.fits' % nifu
    nname = base_name + '_%02d_slice.fits' % nifu
    shutil.copy(name,nname)

  for slice_number in range(1,49):
    # read in the previously created sky spectrum
    skyspec_filename = base_name + '_sky_slice%02d.fits' % slice_number
    skyspec_hdu = pyfits.open(skyspec_filename)
    wav_start = skyspec_hdu[0].header['CRVAL1']
    wav_delta = skyspec_hdu[0].header['CDELT1']
    wav_length = skyspec_hdu[0].header['NAXIS1']
    sky_wave = np.arange(wav_start,wav_start+(wav_length*wav_delta),wav_delta)
    sky_flux = skyspec_hdu[0].data
    func = interp.interp1d(sky_wave,sky_flux)

    # do the sky subtraction
    for ifu in range(1,25):
      print 'Subtracting slice %02d in IFU %02d' % (slice_number, ifu)
      pixname = base_name + '_%02d_slice.fits' % ifu
      hdu = pyfits.open(pixname)
      origin = hdu[7].data
      flux_big = hdu[4].data
      slice = origin & 0x3f
      wave = hdu[3].data[slice==slice_number]
      flux = hdu[4].data[slice==slice_number]
      flux_new = np.zeros(flux.shape,dtype=np.float32)
      reli = ifu-1
      relative = float(relflux_table[reli])

      for i in range(len(wave)):        
        flux_new[i] = flux[i] - (func(wave[i])*relative)
      flux_big[slice==slice_number] = flux_new    

      # modify and write out the Pixtable accordingly
      hdu0 = pyfits.PrimaryHDU(header=hdu[0].header)
      hdu1 = pyfits.ImageHDU(data=hdu[1].data,header=hdu[1].header)
      hdu2 = pyfits.ImageHDU(data=hdu[2].data,header=hdu[2].header)
      hdu3 = pyfits.ImageHDU(data=hdu[3].data,header=hdu[3].header)
      hdu4 = pyfits.ImageHDU(data=flux_big,header=hdu[4].header)
      hdu5 = pyfits.ImageHDU(data=hdu[5].data,header=hdu[5].header)
      hdu6 = pyfits.ImageHDU(data=hdu[6].data,header=hdu[6].header)
      hdu7 = pyfits.ImageHDU(data=hdu[7].data,header=hdu[7].header)
      hdulist = pyfits.HDUList([hdu0,hdu1,hdu2,hdu3,hdu4,hdu5,hdu6,hdu7])
      newpixname = base_name + '_%02d_slice.fits' % ifu
      hdulist.writeto(newpixname,clobber=True)
      hdu.close()
  
    print 'Finished with slice %02d' % slice_number
  
  hduup = pyfits.open(newpixname)
  header0=hduup[0].header
  header0['history'] = 'This file was updated using the relative slice'
  header0['history'] = 'subtraction by T. Urrutia'
  header0['history'] = 'datetime.datetime.now()'
  rel_fluxes.close()
  print 'Totally finished!'

if __name__ == '__main__':
    main()
