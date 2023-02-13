#! /usr/bin/env python
#
# FILE: cubecomb_divided.py
# AUTHORS: Tanya Urrutia, Lutz Wisotzki
# DESCR: Combine cubes previously drizzled onto the same grid, before average
#        ombining there are two kappa-sigma clipping cycles
#

import argparse
import sys,time
import numpy as np
import astropy.io.fits as fits
import astropy.stats as stats
import glob as g
from astropy.table import Table

import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)

parser = argparse.ArgumentParser(description="""
Combine cubes previously drizzled onto the same grid. The cubes are divided into 6, 26 or 60 segments along the wavelength axis. Before average combining each voxel there is are two kappa-sigma clipping routines, the kappa values being based once on the individual variances and once on the median variances - the number of clipped voxel information is written out to a table.""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i","--input",
                    required=True,
                    type=str,
                    help="""
                    Base name of segmented cubes that are to be combined. Example: 'data_sfc_M0416-SF2_01_a.fits', then the base name is 'sfc_M0416-SF2'. Corresponding stat and wht cubes are expected. The combined segmented data cube would then be 'data_sfc_M0416-SF2_comb_a.fits'.
                    """)
parser.add_argument("--kappa1",
                    type=float,
                    default=6.0,
                    help="""
                    First kappa threshold for kappa-sigma clipping based on individual variances.
                    """)
parser.add_argument("--kappa2",
                    type=float,
                    default=20.0,
                    help="""
                    Second kappa threshold for kappa-sigma clipping based on the median of all variances.
                    """)
parser.add_argument("--seg_no",
                    type=int,
                    default=6,
                    help="""
                    Number of segments to be combined, same number that were divided. Possible values are 6, 26 and 60.
                    """)

args = parser.parse_args()

name = args.input
kappa1 = args.kappa1
kappa2 = args.kappa2
seg_no = args.seg_no

if (seg_no == 6):
    segments = ['a','b','c','d','e','f']
elif (seg_no == 26):
    segments = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
elif (seg_no == 60):
    segments = ['aa','ab','ac','ad','ae','af','ag','ah','ai','aj','ak','al','am','an','ao','ap','aq','ar','as','at','au','av','aw','ax','ay','az','ba','bb','bc','bd','be','bf','bg','bh','bi','bj','bk','bl','bm','bn','bo','bp','bq','br','bs','bt','bu','bv','bw','bx','by','bz','ca','cb','cc','cd','ce','cf','cg','ch']
else:
   print('Have not programmed for this number of segments,')
   print('possible values are 6, 26 or 60.')
   sys.exit(1)

# I still leave this in case the program crashes and I want to resume, this must be done manually...
#segments = ['a','b','c','d','e','f']
#segments = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']  
#segments = ['aa','ab','ac','ad','ae','af','ag','ah','ai','aj','ak','al','am','an','ao','ap','aq','ar','as','at','au','av','aw','ax','ay','az','ba','bb','bc','bd','be','bf','bg','bh','bi','bj','bk','bl','bm','bn','bo','bp','bq','br','bs','bt','bu','bv','bw','bx','by','bz','ca','cb','cc','cd','ce','cf','cg','ch']

starttime = time.time()

################

print('Combining %d cubes for %s' % (seg_no,name))
print('...')

for se in segments:
    print('Going with %s' % se)
    passtime='('+str(round(time.time()-starttime,3))+'s)' 
    print('... '+passtime)
    datafiles = sorted(g.glob('./data_%s_??_%s.fits' % (name,se)))
    statfiles = sorted(g.glob('./stat_%s_??_%s.fits' % (name,se)))
    whtfiles = sorted(g.glob('./wht_%s_??_%s.fits' % (name,se)))
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
    passtime='('+str(round(time.time()-starttime,3))+'s)' 
    print('... '+passtime)
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
    passtime='('+str(round(time.time()-starttime,3))+'s)' 
    print('... '+passtime)
    cliplist = []
    statcliplist = []
    whtcliplist = []
    numclip = np.zeros(len(datalist),np.int64)
    for i in range(len(datalist)):
        clip = np.zeros(datalist[i].shape,dtype=np.float32)
        # mask criterion
        mask1 = (datalist[i] < datmedian + (kappa1*stdlist[i])) & (datalist[i] > datmedian - (kappa1*stdlist[i]))
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
    passtime='('+str(round(time.time()-starttime,3))+'s)' 
    print('... '+passtime)
    # new weight cubes
    weigh_comb = np.sum(whtcliplist,axis=0)
    stat_weigh = np.square(weigh_comb)
    # data combine
    data_comb = np.nanmean(cliplist,axis=0)
    # stat combine
    sumstat = np.nansum(statcliplist,axis=0)
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
    passtime='('+str(round(time.time()-starttime,3))+'s)' 
    print('... '+passtime)
    data_comb = np.float32(data_comb)
    stat_comb = np.float32(stat_comb)
    weigh_comb = np.uint8(weigh_comb)
    dcombhdu = fits.PrimaryHDU(data=data_comb)
    scombhdu = fits.PrimaryHDU(data=stat_comb)
    wcombhdu = fits.PrimaryHDU(data=weigh_comb)
    
    dcombhdu.writeto('data_%s_comb_%s.fits' % (name,se),overwrite=True)
    scombhdu.writeto('stat_%s_comb_%s.fits' % (name,se),overwrite=True)
    wcombhdu.writeto('wht_%s_comb_%s.fits' % (name,se),overwrite=True)
    t.write('clip_%s_comb_%s.fits' % (name,se),overwrite=True)

passtime='('+str(round(time.time()-starttime,3))+'s)'
print('Done. '+passtime)

