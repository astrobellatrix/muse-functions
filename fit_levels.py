#!/usr/bin/python

import sys, datetime
import astropy.io.fits as pyfits
import astropy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gauss(x, *p):
    A, mu, sigma, b = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))+b

def avsigc(indata,sigma):
   a = stats.sigma_clip(indata,sigma=sigma,iters=3)
   mn = np.median(a[-a.mask])
   lena = np.ma.MaskedArray.count(a)
   return mn, lena

def main():
  # define base names to read the pixel tables, set up...
  base_name = sys.argv[1:]
  if not base_name:
    print 'Please provide a base name, e.g. PIXTABLE_TNJ1338_12'
    sys.exit(1)
  base_name = str(base_name[0])
  relflux_filename = base_name + '_relflux.txt'
  rel_fluxes = open(relflux_filename,'w')
  big_rel_fluxes = np.zeros((48,24),dtype=np.float32)

  # mask brightest pixels for sky creation
  fov_name = '../mid/IMAGE_FOV_' + base_name[9:] + '.fits'
  fovhdu = pyfits.open(fov_name)
  fov = fovhdu[1].data
  fovhead = fovhdu[1].header
  xstart = fovhead['CRVAL1']
  ystart = fovhead['CRVAL2']
  lowfov = np.percentile(fov,5)
  hifov = np.percentile(fov,85)
  fovmask = np.ones(fov.shape,dtype=np.int32)
  fovmask[fov < lowfov] = 0
  fovmask[fov > hifov] = 0
  fovfile = pyfits.PrimaryHDU(data=fovmask,header=fovhead)
  fovname = base_name + '_fovmask.fits'
  fovfile.writeto(fovname,clobber=True)

  for slice_number in range(1,49):
    # set up variables, this is not necessary in python, but nice...
    s_wave = np.array([],dtype=np.float32)
    s_flux = np.array([],dtype=np.float32)
    s_xp = np.array([],dtype=np.float32)
    s_yp = np.array([],dtype=np.float32)
    rel_flux = np.array([],dtype=np.float32)
    ifu_int5577 = np.array([],dtype=np.float32)
    ifu_int6300 = np.array([],dtype=np.float32)
    ifu_int8943 = np.array([],dtype=np.float32)

    print 'Reading in IFU info for slice %02d' % slice_number
    for ifu in range(1,25):
      # set up the figure for 3 emission line (5577,6300,8942) plots.
      # In the end there will be 48x24 figures
      fignumber = int(str(slice_number)+str(ifu))
      plt.figure(fignumber,dpi=300)
      plt.axis([-4,4,0,1300])
      
      # read in the ifu info
      pixname = base_name + '_%02d.fits' % ifu
      hdu = pyfits.open(pixname)
      xpos = hdu[1].data-xstart
      ypos = hdu[2].data-ystart
      xp = xpos.astype(int)
      yp = ypos.astype(int)
      wave = hdu[3].data
      flux = hdu[4].data
      origin = hdu[7].data
      slice = origin & 0x3f
      
      # append all the flux for the relevant slice. 
      # This is for the average slice spectrum later.
      s_wave = np.append(s_wave,wave[slice==slice_number])
      s_flux = np.append(s_flux,flux[slice==slice_number])
      s_xp = np.append(s_xp,xp[slice==slice_number])
      s_yp = np.append(s_yp,yp[slice==slice_number])

      # let's just deal with the current ifu and the current slice
      si_wave = wave[slice==slice_number]
      si_flux = flux[slice==slice_number]
      si_xp = xp[slice==slice_number]
      si_yp = yp[slice==slice_number]

      # Make the bright object mask
      #print 'Making out brightest objects...'
      imask = np.zeros(len(si_wave),dtype=bool)

      # Masking out the brightest and dimmest regions of X-Y space
      for iele in range(len(si_wave)):
        if (fovmask[si_yp[iele]-1][si_xp[iele]-1] == 1):
          imask[iele] = True
   
      # Now create all the masks for the 3 different averages we have
      # of the emission line for integration.
      # I tried to use only one strong line, no blend.
      mask_for_5577_o = (si_wave > 5572) & (si_wave < 5582)
      mask_for_6300_o = (si_wave > 6295) & (si_wave < 6305)
      mask_for_8943_o = (si_wave > 8938) & (si_wave < 8948)
      mask_for_5577 = (si_wave > 5572) & (si_wave < 5582) & imask
      mask_for_6300 = (si_wave > 6295) & (si_wave < 6305) & imask
      mask_for_8943 = (si_wave > 8938) & (si_wave < 8948) & imask

      # The sky lines will now be fitted with a Gaussian. 
      # Initial conditions are approx. relative flux, theoretical center 
      # of skyline and 1.2 Angstrom sigma (should be fine)
      # I will then use the integrated flux over that Gaussian
      # Best would be to model the LSF, of course, but I fear it is too 
      # sparse points on only one line/ifu/slice. The integrated flux should 
      # account for the LSF effects to zeroth order as flux is conserved. 
      # Anything that broadens the wings over a Gaussian will then flatten at 
      # the maximum, so using the integrated flux is good.
      pini5577 = [1250.,5577.338,1.2,30.]
      pini6300 = [850.,6300.304,1.2,30.]
      pini8943 = [400.,8943.395,1.2,30.]
      plt.plot(si_wave[mask_for_5577]-5577.338,si_flux[mask_for_5577],'b.')
      plt.plot(si_wave[mask_for_6300]-6300.304,si_flux[mask_for_6300],'g.')
      plt.plot(si_wave[mask_for_8943]-8943.395,si_flux[mask_for_8943],'r.')

      if (len(si_wave[mask_for_5577]) < 50):
        coeff5577, cov5577 = curve_fit(gauss, si_wave[mask_for_5577_o], si_flux[mask_for_5577_o],p0=pini5577)
        print 'warning: slice %02d, IFU %02d had too many masked out values (5577) (%02d)' % (slice_number,ifu,len(si_wave[mask_for_5577]))
      else:
        coeff5577, cov5577 = curve_fit(gauss, si_wave[mask_for_5577], si_flux[mask_for_5577],p0=pini5577)
      if (len(si_wave[mask_for_6300]) < 50):
        coeff6300, cov6300 = curve_fit(gauss, si_wave[mask_for_6300_o], si_flux[mask_for_6300_o],p0=pini6300)
        print 'warning: slice %02d, IFU %02d had too many masked out values (6300) (%02d)' % (slice_number,ifu,len(si_wave[mask_for_6300]))
      else:
        coeff6300, cov6300 = curve_fit(gauss, si_wave[mask_for_6300], si_flux[mask_for_6300],p0=pini6300)
      if (len(si_wave[mask_for_8943]) < 50):
        coeff8943, cov8943 = curve_fit(gauss, si_wave[mask_for_8943_o], si_flux[mask_for_8943_o],p0=pini8943) 
        print 'warning: slice %02d, IFU %02d had too many masked out values (8943) (%02d)' % (slice_number,ifu,len(si_wave[mask_for_8943]))
      else:
        coeff8943, cov8943 = curve_fit(gauss, si_wave[mask_for_8943], si_flux[mask_for_8943],p0=pini8943)

      wavplot = np.arange(-4.5,4.5, 0.01)
      plt.plot(wavplot,gauss(wavplot+5577.338,*coeff5577),color='b',linestyle='-',label='5577 line')
      plt.plot(wavplot,gauss(wavplot+6300.304,*coeff6300),color='g',linestyle='-',label='6300 line')
      plt.plot(wavplot,gauss(wavplot+8943.395,*coeff8943),color='r',linestyle='-',label='8943 line')

      int_flux5577 = sum([gauss((5573+(0.1*w)),*coeff5577)/10 for w in range(80)])
      int_flux6300 = sum([gauss((6296+(0.1*w)),*coeff6300)/10 for w in range(80)])
      int_flux8943 = sum([gauss((8939+(0.1*w)),*coeff8943)/10 for w in range(80)])
      ifu_int5577 = np.append(ifu_int5577,int_flux5577)
      ifu_int6300 = np.append(ifu_int6300,int_flux6300)
      ifu_int8943 = np.append(ifu_int8943,int_flux8943)
      plt.ylabel('Counts')
      plt.xlabel('Difference in Angstrom rel. to theoretical line center')
      plt.title('Line fit for Slice %02d IFU %02d' % (slice_number,ifu))
      plt.legend(bbox_to_anchor=(0.9,0.95),borderaxespad=0.,prop={'size':10})
      figname1 = 'line_slice%02d_ifu%02d.png' % (slice_number,ifu)
      figname2 = 'png/' + base_name + "_" + figname1
      plt.savefig(figname2,format='png')
      plt.close(fignumber)
      hdu.close()
    
    # now we can take the relative values in comparison to the other ifus
    rel_int5577 = ifu_int5577/np.median(ifu_int5577)
    rel_int6300 = ifu_int6300/np.median(ifu_int6300)
    rel_int8943 = ifu_int8943/np.median(ifu_int8943)
    print rel_int5577
    print rel_int6300
    print rel_int8943
    
    # take the median of these relative values for each ifu 
    # note that median in 4 elements is the average of the middle 2
    for nifu in range(24):
      ifu_rel = np.median(np.array([rel_int5577[nifu],rel_int6300[nifu],rel_int8943[nifu]]))
      rel_flux = np.append(rel_flux,ifu_rel)
    print rel_flux
    big_rel_fluxes[slice_number-1] = rel_flux

    # Make the bright object mask
    print 'Making the bright object mask for sky subtraction...'
    mask = np.zeros(len(s_flux),dtype=bool)

    # Masking out the brightest and dimmest regions of X-Y space
    for ele in range(len(s_flux)):
      if (fovmask[int(s_yp[ele])-1][int(s_xp[ele])-1] > 0):
        mask[ele] = True

    # define bins
    print 'The relevant pixels are set up! Length %d' % len(s_flux)
    print 'Now on to create the spectrum for slice %02d' % slice_number
    binsize = 0.2
    global_wave = np.arange(4750.+(binsize/2),9350.+(binsize/2),binsize)
    global_flux = np.array([],dtype=np.float32)
    global_len = np.array([],dtype=np.float32)
    s2_flux = s_flux[mask]
    s2_wave = s_wave[mask]
    print 'Masked down to %d pixels' % len(s2_flux)

    # and now create the average flux in each bin
    minute = 4/binsize
    print 'This takes a while...%2d minutes' % minute
    for each in global_wave:
      mask2 = (s2_wave > each-(binsize/2)) & (s2_wave < each+(binsize/2))
      values = s2_flux[mask2]
      av, lena = avsigc(values,2.0)
      global_len = np.append(global_len, lena)
      if (lena != 0):
        global_flux = np.append(global_flux,av)
      else:
        global_flux = np.append(global_flux,global_flux[-1])
    
    # plot the average function (this will be our sky)
    plt.figure(1000+slice_number,dpi=300)
    plt.plot(global_wave,global_flux,'k-')
    figname = 'png/' + base_name + '_sky_slice%02d.png' % slice_number
    plt.savefig(figname,format='png')
    plt.close(1000+slice_number)

    # add points above and below to sample spectrum later fully from 4750-9300
    global_flux = np.insert(global_flux,0,global_flux[0])
    global_flux = np.append(global_flux,global_flux[-1])
    global_len = np.insert(global_len,0,global_len[0])
    global_len = np.append(global_len,global_flux[-1])
    global_wave = np.insert(global_wave,0,4750-(binsize/2))
    global_wave = np.append(global_wave,9450+(binsize/2))
    
    # write sky to a fits file
    header = pyfits.Header()
    header['SIMPLE'] = (True,'Fits standard')
    header['BITPIX'] = (-32,'Bits per pixel')
    header['NAXIS'] = (1,'Number of axes')
    header['NAXIS1'] = (len(flux),'Axis length')
    header['BUNIT'] ='10**(-20)*erg/s/cm**2/Angstrom'
    spec_title = 'Sky Spectrum Slice ' + str(slice_number) + ' ' + base_name
    header['OBJECT'] = (spec_title,'Name of the spectrum')
    header['DATE'] = (str(datetime.datetime.now()),'Date FITS file was generated')
    header['CRPIX1'] = 1.
    header['CRVAL1'] = 4750.-(binsize/2)
    header['CDELT1'] = binsize
    header['CUNIT1'] = 'Angstrom'
    header['DC-FLAG'] = 0
    hduout = pyfits.PrimaryHDU(data=global_flux,header=header)
    outspec = base_name + '_sky_slice%02d.fits' % slice_number
    hduout.writeto(outspec,clobber=True)

    # write numspec used to a fits file
    header = pyfits.Header()
    header['SIMPLE'] = (True,'Fits standard')
    header['BITPIX'] = (-32,'Bits per pixel')
    header['NAXIS'] = (1,'Number of axes')
    header['NAXIS1'] = (len(flux),'Axis length')
    header['BUNIT'] ='10**(-20)*erg/s/cm**2/Angstrom'
    spec_title = 'Len Spectrum Slice ' + str(slice_number) + ' ' + base_name
    header['OBJECT'] = (spec_title,'Name of the spectrum')
    header['DATE'] = (str(datetime.datetime.now()),'Date FITS file was generated')
    header['CRPIX1'] = 1.
    header['CRVAL1'] = 4750.-(binsize/2)
    header['CDELT1'] = binsize
    header['CUNIT1'] = 'Angstrom'
    header['DC-FLAG'] = 0
    hduout = pyfits.PrimaryHDU(data=global_len,header=header)
    outspec = base_name + '_sky_slen%02d.fits' % slice_number
    hduout.writeto(outspec,clobber=True)

    print 'Finished with slice %02d' % slice_number

  # now that we have a big 48x24 table of relative fluxes, we median  
  # the fluxes per ifu and write to a file
  per_fluxes = np.percentile(big_rel_fluxes,50,axis=0)
  rel_fluxes.write(str(per_fluxes))
  rel_fluxes.close()
  print 'Totally finished!'

if __name__ == '__main__':
    main()
