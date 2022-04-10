#!/usr/bin/python

import astropy.io.fits as fits
import numpy as np
from astropy.table import Table
import glob as g

def main():
  print 'reading'
  data01 = fits.getdata('data_sfc_M0416-SF2_comb_a.fits')
  data02 = fits.getdata('data_sfc_M0416-SF2_comb_b.fits')
  data03 = fits.getdata('data_sfc_M0416-SF2_comb_c.fits')
  data04 = fits.getdata('data_sfc_M0416-SF2_comb_d.fits')
  data05 = fits.getdata('data_sfc_M0416-SF2_comb_e.fits')
  data06 = fits.getdata('data_sfc_M0416-SF2_comb_f.fits')

  stat01 = fits.getdata('stat_sfc_M0416-SF2_comb_a.fits')
  stat02 = fits.getdata('stat_sfc_M0416-SF2_comb_b.fits')
  stat03 = fits.getdata('stat_sfc_M0416-SF2_comb_c.fits')
  stat04 = fits.getdata('stat_sfc_M0416-SF2_comb_d.fits')
  stat05 = fits.getdata('stat_sfc_M0416-SF2_comb_e.fits')
  stat06 = fits.getdata('stat_sfc_M0416-SF2_comb_f.fits')

  wht01 = fits.getdata('wht_sfc_M0416-SF2_comb_a.fits')
  wht02 = fits.getdata('wht_sfc_M0416-SF2_comb_b.fits')
  wht03 = fits.getdata('wht_sfc_M0416-SF2_comb_c.fits')
  wht04 = fits.getdata('wht_sfc_M0416-SF2_comb_d.fits')
  wht05 = fits.getdata('wht_sfc_M0416-SF2_comb_e.fits')
  wht06 = fits.getdata('wht_sfc_M0416-SF2_comb_f.fits')

  clipfiles = sorted(g.glob('clip_sfc_M0416-SF2_comb_?.fits'))
  cliplist = []
  for c in clipfiles:
    clip = fits.getdata(c)
    cliplist.append(clip)
  f_name = []
  tot_vox = np.zeros(cliplist[0].shape[0])
  clip_vox = np.zeros(cliplist[0].shape[0])
  per_vox = np.zeros(cliplist[0].shape[0])

  hdu1 = fits.open('DATACUBE_SFC_M0416-SF2_01.fits')
  
  print 'combining'
  data_comb = np.concatenate((data01,data02,data03,data04,data05,data06),axis=0)
  stat_comb = np.concatenate((stat01,stat02,stat03,stat04,stat05,stat06),axis=0)
  wht_comb = np.concatenate((wht01,wht02,wht03,wht04,wht05,wht06),axis=0)
   
  for i in range(cliplist[0].shape[0]):
    f_name.append(cliplist[0][i][0][2:-7])
    for k in range(len(clipfiles)):
      tot_vox[i] += cliplist[k][i][1]
      clip_vox[i] += cliplist[k][i][2]
    per_vox[i] = 100 * clip_vox[i] / tot_vox[i]
  f_name = f_name[:-1] + ['Total averages']
  tot_vox[-1] = np.mean(tot_vox[:-1])
  clip_vox[-1] = np.mean(clip_vox[:-1])
  per_vox[-1] = np.mean(per_vox[:-1])
  t = Table([f_name,tot_vox,clip_vox,per_vox], names=('FIELD_NAME','TOTAL_VOXELS','CLIPPED_VOXELS','PERCENTAGE_CLIPPED'))

  print 'whitelight'
  ma_data_comb = np.ma.masked_array(data_comb,np.isnan(data_comb))
  mean = np.mean(ma_data_comb,axis=0)
  white_comb = mean.filled(np.nan)

  print 'writing'
  prihdu = fits.PrimaryHDU(header=hdu1[0].header)
  datac = fits.ImageHDU(data=data_comb,header=hdu1[1].header)
  statc = fits.ImageHDU(data=stat_comb,header=hdu1[2].header)
  head3 = hdu1[2].header
  head3.set('BITPIX',8)
  head3.set('EXTNAME','EXP    ','This extension contains the exposure number')
  del head3['HDUCLAS2']
  del head3['HDUCLAS3']
  del head3['BUNIT']
  head3.set('OBJECT','M0416-SF2 (EXP)')
  weighc = fits.ImageHDU(data=wht_comb,header=head3)
  whitec = fits.ImageHDU(data=white_comb,header=hdu1[3].header)
  clipc = fits.BinTableHDU(data=t)
  hdulist = fits.HDUList([prihdu,datac,statc,weighc,whitec,clipc])
  hdulist.writeto('DATACUBE_SFC_M0416-SF2_imcomb.new.fits')
  whit = fits.PrimaryHDU(data=white_comb,header=hdu1[3].header)
  whit.writeto('white_sfc_M0416-SF2_imcomb.new.fits')
  

if __name__ == '__main__':
  main()
