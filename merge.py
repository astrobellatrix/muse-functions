#!/usr/bin/env python
# 
# FILE:   merge.py
# AUTHOR: Tanya Urrutia
# DESCR.: Merge all combined segmented datacubes into one.
#

import argparse
import time
import astropy.io.fits as fits
import numpy as np
from astropy.table import Table
import glob as g

parser = argparse.ArgumentParser(description="""
Merge all combined segmented datacubes into one.""",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i","--input",
                    required=True,
                    type=str,
                    help="""
                    Base name of cube that is combined. Example: 'data_sf_M0416-SF_comb_[x].fits', then the base name is 'SF_M0416-SF'.
                    """)
parser.add_argument("-z","--zero",
                    required=True,
                    type=str,
                    help="""
                    Name of the DATACUBE for which we will take the headers.
                    """)
parser.add_argument("-o","--output",
                    type=str,
                    default=None,
                    help="""
                    Name of the output combined FITS file. [Default: `DATACUBE_[base_name]_imcomb.fits`, where base name is the input base name used.]
                    """)
parser.add_argument("-c","--clipout",
                    type=str,
                    default=None,
                    help="""
                    Name of the output clipstats table FITS file. [Default: `clipstats_[base_name].fits`, where base name is the input base name used.]
                    """)

args = parser.parse_args()

name = args.input
zero = args.zero

if args.output == None:
    out_fits = 'DATACUBE_%s_imcomb.fits' % name
else:
    out_fits = args.output

if args.clipout == None:
    clipout = 'clipstats_%s.fits' % name
else:
    clipout = args.clipout

starttime = time.time()

################

print('reading')
datafiles = sorted(g.glob('./data_%s_comb_*.fits' % name))
statfiles = sorted(g.glob('./stat_%s_comb_*.fits' % name))
whtfiles = sorted(g.glob('./wht_%s_comb_*.fits' % name))
clipfiles = sorted(g.glob('clip_%s_comb_*.fits' % name))
datalist = []
statlist = []
whtlist = []
cliplist = []

for d in datafiles:
    data = fits.getdata(d)
    datalist.append(data)
  
for s in statfiles:
    stat = fits.getdata(s)
    statlist.append(stat)

for w in whtfiles:
    wht = fits.getdata(w)
    whtlist.append(wht)

for c in clipfiles:
    clip = fits.getdata(c)
    cliplist.append(clip)
  
f_name = []
tot_vox = np.zeros(cliplist[0].shape[0])
clip_vox = np.zeros(cliplist[0].shape[0])
per_vox = np.zeros(cliplist[0].shape[0])

hdu1 = fits.open(zero)
data_comb = datalist[0]
stat_comb = statlist[0]
wht_comb = whtlist[0]
  
print('combining')
for i in range(1,len(datafiles)):
    data_comb = np.concatenate((data_comb,datalist[i]),axis=0)
    stat_comb = np.concatenate((stat_comb,statlist[i]),axis=0)
    wht_comb = np.concatenate((wht_comb,whtlist[i]),axis=0)
   
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

print('whitelight')
ma_data_comb = np.ma.masked_array(data_comb,np.isnan(data_comb))
mean = np.mean(ma_data_comb,axis=0)
white_comb = mean.filled(np.nan)

print('writing')
prihdu = fits.PrimaryHDU(header=hdu1[0].header)
datac = fits.ImageHDU(data=data_comb,header=hdu1[1].header)
statc = fits.ImageHDU(data=stat_comb,header=hdu1[2].header)
head3 = hdu1[2].header
head3.set('BITPIX',8)
head3.set('EXTNAME','EXP    ','This extension contains the exposure number')
del head3['HDUCLAS2']
del head3['HDUCLAS3']
del head3['BUNIT']
head3.set('OBJECT','%s (EXP)' % name)
weighc = fits.ImageHDU(data=wht_comb,header=head3)
whitec = fits.ImageHDU(data=white_comb,header=hdu1[3].header)
hdulist = fits.HDUList([prihdu,datac,statc,weighc,whitec])
hdulist.writeto(out_fits)

whit = fits.PrimaryHDU(data=white_comb,header=hdu1[3].header)
whit.writeto('white_%s_imcomb.fits' % name)
  
prihdr = fits.Header()
prihdr['COMMENT'] = 'Clip statistics for %s' % name
clipc = fits.BinTableHDU(data=t)
thdulist = fits.HDUList([prihdu,clipc])
thdulist.writeto(clipout,overwrite=True)

passtime='('+str(round(time.time()-starttime,3))+'s)'
print('Done. '+passtime)
