#! /usr/bin/env python

import sys
import numpy as np
import astropy.io.fits as fits
import astropy.stats as stats
import glob as g
from astropy.table import Table

def main():
  #segments = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
  segments = ['a','b','c','d','e','f']
  for se in segments:
    print('Going with %s' % se)
    datafiles = sorted(g.glob('./data_sfc_M0416-SF2_??_%s.fits' % se))
    statfiles = sorted(g.glob('./stat_sfc_M0416-SF2_??_%s.fits' % se))
    whtfiles = sorted(g.glob('./wht_sfc_M0416-SF2_??_%s.fits' % se))
    datalist = []
    statlist = []
    whtlist = []
  
    for d in datafiles:
      data = fits.getdata(d)
      datalist.append(data)
  
    for s in statfiles:
      stat = fits.getdata(s)
      statlist.append(stat)

    for w in whtfiles:
      wht = fits.getdata(w)
      whtlist.append(wht)

    # check if I can do the combination
    print('combining %02d files' % len(datalist))
    if (len(datalist) != len(statlist)) or (len(datalist) != len(whtlist)):
      print('Not the same amount of files')
      sys.exit(1)
    for i in range(len(datalist)):
      if (datalist[0].shape != datalist[i].shape):
        print('Data arrays do not match in size')
        sys.exit(1)
      if (datalist[0].shape != statlist[i].shape):
        print('Data arrays do not match in size with stat arrays')
        sys.exit(1)
      if (datalist[0].shape != whtlist[i].shape):
        print('Data arrays do not match in size with weight arrays')
        sys.exit(1)

    # initial weight cube sum for number of voxels
    weigh_comb = np.sum(whtlist,axis=0)
    vox_num = np.zeros(len(whtlist),dtype=np.int64)
    for v in range(len(vox_num)):
      vox_num[v] = len(whtlist[v][whtlist[v]==1])

    # setting big outliers in data to np.nan
    print('Defining threshold cube for clipping')
    kappa = 5.0
    kappa2 = 20.0
    # median of data to find out how far outliers are
    datmedian = np.nanmedian(datalist,axis=0)
    # get corrected sigma at each voxel
    # approx. correction factor due to resampling variance cube
    #corr_factor = 1.7
    corr_factor = 1.0  # individual voxels don't need correction
    stdlist = np.sqrt(statlist) * corr_factor
    stdmedian = np.nanmedian(stdlist,axis=0)
    # from this kappa2*sigma get the secondary threshold cube
    thresh = kappa2 * stdmedian

    print('kappa-sigma clipping')
    cliplist = []
    statcliplist = []
    whtcliplist = []
    numclip = np.zeros(len(datalist),np.int64)
    for i in range(len(datalist)):
        clip = np.zeros(datalist[i].shape,dtype=np.float32)
        # mask criterion
        mask1 = (datalist[i] < datmedian + (kappa*stdlist[i])) & (datalist[i] > datmedian - (kappa*stdlist[i]))
        mask2 = (datalist[i] < datmedian + thresh) & (datalist[i] > datmedian - thresh)
        critmask = mask1 & mask2
        # mask data
        clip[critmask] = datalist[i][critmask]
        clip[np.logical_not(critmask)] = np.nan
        d = len(clip[np.isnan(clip)]) - len(datalist[i][np.isnan(datalist[i])])
        numclip[i] = d
        # mask stat
        statclip = np.zeros(statlist[i].shape,dtype=np.float32)
        statclip[critmask] = statlist[i][critmask]
        statclip[np.logical_not(critmask)] = np.nan
        # new weight cubes
        wht = np.zeros(whtlist[i].shape,dtype=np.uint8)
        wht[critmask] = 1
        whtcliplist.append(wht)
        print('cube %d had %d pixels masked' % ((i+1),d))
        cliplist.append(clip)
        statcliplist.append(statclip)

    # make weighted average of the data and of stat
    # division by zero gets assigned a nan
    print('combining cubes')
    # new weight cubes
    weigh_comb = np.sum(whtcliplist,axis=0)
    stat_weigh = np.square(weigh_comb)
    # data combine
    data_comb = np.nanmean(cliplist,axis=0)
    # stat combine
    sumstat = np.sum(statcliplist,axis=0)
    stat_comb = np.divide(sumstat,stat_weigh)
    stat_comb[stat_comb==np.inf] = np.nan
    
    # for no sigma clipping uncomment below and comment above paragraph
    #print('combining cubes')
    #weigh_comb = np.sum(whtlist,axis=0)
    #stat_weigh = np.square(weigh_comb)
    #data_comb = np.nanmean(datalist,axis=0)
    #sumstat = np.sum(statlist,axis=0)
    #stat_comb = np.divide(sumstat,stat_weigh)
    #stat_comb[stat_comb==np.inf] = np.nan

    # Table of clipped voxel statistics for QC
    per = 100 * numclip / vox_num
    datafiles = np.append(datafiles,'Total averages')
    vox_num = np.append(vox_num,np.mean(vox_num))
    numclip = np.append(numclip,np.mean(numclip))
    per = np.append(per,np.mean(per))
    t = Table([datafiles,vox_num,numclip,per], names=('FIELD_NAME','TOTAL_VOXELS','CLIPPED_VOXELS','PERCENTAGE_CLIPPED'))

    print('writing')
    data_comb = np.float32(data_comb)
    stat_comb = np.float32(stat_comb)
    weigh_comb = np.uint8(weigh_comb)
    dcombhdu = fits.PrimaryHDU(data=data_comb)
    scombhdu = fits.PrimaryHDU(data=stat_comb)
    wcombhdu = fits.PrimaryHDU(data=weigh_comb)
    
    dcombhdu.writeto('data_sfc_M0416-SF2_comb_%s.fits' % se,overwrite=True)
    scombhdu.writeto('stat_sfc_M0416-SF2_comb_%s.fits' % se,overwrite=True)
    wcombhdu.writeto('wht_sfc_M0416-SF2_comb_%s.fits' % se,overwrite=True)
    t.write('clip_sfc_M0416-SF2_comb_%s.fits' % se,overwrite=True)

if __name__ == '__main__':
  main()
