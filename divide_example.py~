#!/usr/bin/python

import astropy.io.fits as fits
import glob as g

def main():
  filelist = g.glob('./DATACUBE_SFC_M0416-SF2_??.fits')
  for f in filelist:
    # read in stuff
    hdu = fits.open(f)
    data = hdu[1].data
    stat = hdu[2].data
    whdu = fits.open(f[:-5] + '.wht.fits')
    wht = whdu[0].data
    wav = data.shape[0]
    divi = wav/6 # divide by 6

    # first one
    data_a = data[:divi,:,:]
    stat_a = stat[:divi,:,:]
    wht_a = wht[:divi,:,:]
    dnew_a = fits.PrimaryHDU(data=data_a)
    snew_a = fits.PrimaryHDU(data=stat_a)
    wnew_a = fits.PrimaryHDU(data=wht_a)
    dname_a = 'data_sfc_M0416-SF2_' + f[25:27] + '_a.fits'
    sname_a = 'stat_sfc_M0416-SF2_' + f[25:27] + '_a.fits'
    wname_a = 'wht_sfc_M0416-SF2_' + f[25:27] + '_a.fits'
    dnew_a.writeto(dname_a)
    snew_a.writeto(sname_a)
    wnew_a.writeto(wname_a)

    # second one
    data_b = data[divi:divi*2,:,:]
    stat_b = stat[divi:divi*2,:,:]
    wht_b = wht[divi:divi*2,:,:]
    dnew_b = fits.PrimaryHDU(data=data_b)
    snew_b = fits.PrimaryHDU(data=stat_b)
    wnew_b = fits.PrimaryHDU(data=wht_b)
    dname_b = 'data_sfc_M0416-SF2_' + f[25:27] + '_b.fits'
    sname_b = 'stat_sfc_M0416-SF2_' + f[25:27] + '_b.fits'
    wname_b = 'wht_sfc_M0416-SF2_' + f[25:27] + '_b.fits'
    dnew_b.writeto(dname_b)
    snew_b.writeto(sname_b)
    wnew_b.writeto(wname_b)

    # third one
    data_c = data[divi*2:divi*3,:,:]
    stat_c = stat[divi*2:divi*3,:,:]
    wht_c = wht[divi*2:divi*3,:,:]
    dnew_c = fits.PrimaryHDU(data=data_c)
    snew_c = fits.PrimaryHDU(data=stat_c)
    wnew_c = fits.PrimaryHDU(data=wht_c)
    dname_c = 'data_sfc_M0416-SF2_' + f[25:27] + '_c.fits'
    sname_c = 'stat_sfc_M0416-SF2_' + f[25:27] + '_c.fits'
    wname_c = 'wht_sfc_M0416-SF2_' + f[25:27] + '_c.fits'
    dnew_c.writeto(dname_c)
    snew_c.writeto(sname_c)
    wnew_c.writeto(wname_c)

    # fourth one
    data_d = data[divi*3:divi*4,:,:]
    stat_d = stat[divi*3:divi*4,:,:]
    wht_d = wht[divi*3:divi*4,:,:]
    dnew_d = fits.PrimaryHDU(data=data_d)
    snew_d = fits.PrimaryHDU(data=stat_d)
    wnew_d = fits.PrimaryHDU(data=wht_d)
    dname_d = 'data_sfc_M0416-SF2_' + f[25:27] + '_d.fits'
    sname_d = 'stat_sfc_M0416-SF2_' + f[25:27] + '_d.fits'
    wname_d = 'wht_sfc_M0416-SF2_' + f[25:27] + '_d.fits'
    dnew_d.writeto(dname_d)
    snew_d.writeto(sname_d)
    wnew_d.writeto(wname_d)

    # fifth one
    data_e = data[divi*4:divi*5,:,:]
    stat_e = stat[divi*4:divi*5,:,:]
    wht_e = wht[divi*4:divi*5,:,:]
    dnew_e = fits.PrimaryHDU(data=data_e)
    snew_e = fits.PrimaryHDU(data=stat_e)
    wnew_e = fits.PrimaryHDU(data=wht_e)
    dname_e = 'data_sfc_M0416-SF2_' + f[25:27] + '_e.fits'
    sname_e = 'stat_sfc_M0416-SF2_' + f[25:27] + '_e.fits'
    wname_e = 'wht_sfc_M0416-SF2_' + f[25:27] + '_e.fits'
    dnew_e.writeto(dname_e)
    snew_e.writeto(sname_e)
    wnew_e.writeto(wname_e)    

    # sixth one
    data_f = data[divi*5:,:,:]
    stat_f = stat[divi*5:,:,:]
    wht_f = wht[divi*5:,:,:]
    dnew_f = fits.PrimaryHDU(data=data_f)
    snew_f = fits.PrimaryHDU(data=stat_f)
    wnew_f = fits.PrimaryHDU(data=wht_f)
    dname_f = 'data_sfc_M0416-SF2_' + f[25:27] + '_f.fits'
    sname_f = 'stat_sfc_M0416-SF2_' + f[25:27] + '_f.fits'
    wname_f = 'wht_sfc_M0416-SF2_' + f[25:27] + '_f.fits'
    dnew_f.writeto(dname_f)
    snew_f.writeto(sname_f)
    wnew_f.writeto(wname_f)

    hdu.close()

if __name__ == '__main__':
  main()
