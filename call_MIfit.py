# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:40:53 2018

@author: kalender
"""
from __future__ import absolute_import
import numpy as np
from astropy.table import vstack
from astropy.io import ascii
from MI_nestlefit import MI_model, nest_lc
import sncosmo
import os
import copy
import filters
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.integrate import quad
phot_d = None
###########################################################################
# DATA
sn = '16geu'    # ground data
fname = '%s_clean.dat'%(sn)   # Cleaned light curve

##############################################################################
phot_d = sncosmo.read_lc(fname) # read the ground data
phot_d['imageid'] = 0   # ground data have imageid = 0!

##############################################################################
# Here the error bars of the grond data can be adjusted. 
lenserr = 0.00  # set to 0.05 when using gravlens model
#phot_d['fluxerr'] = np.sqrt(phot_d['fluxerr']**2+(phot_d['flux']*lenserr)**2)
"""
for i in range(0,len(phot_d['fluxerr'])):   # scale error bars 
    if phot_d['fluxerr'][i] > phot_d['flux'][i]*0.2:
        phot_d['fluxerr'][i] = phot_d['flux'][i]*0.2
"""


##############################################################################
# Set up the data of the HST to be read in
sumdata = False
lcfname = { # the names of the data files
    'f625w' : 'lightcurve_f625w_b00_fitalljoint.csv',
    'f814w' : 'lightcurve_f814w_b00_fitalljoint.csv',
    'f110w' : 'lightcurve_f110w_b00_fitwfc3ir.csv',
    'f160w' : 'lightcurve_f160w_b00_fitwfc3ir.csv',
    'f390w' : 'lightcurve_f390w_b00_fitalljoint.csv',
    'f475w' : 'lightcurve_f475w_b00_fitalljoint.csv'
}
# These are the bands that are cosidered!
for n,band in enumerate(['f625w','f814w','f390w','f475w','f160w']) :
#for n,band in enumerate(['f625w']) :  
    fname = lcfname[band]                                    
    c = ascii.read(os.path.join('photometry/lightcurve',fname))
    flux,ferr2 = np.zeros(len(c)),np.zeros(len(c))
    c['zpsys']   = 'ab'
    if band in ['f625w','f814w']:
        c['band'] = 'uv%s'%(band)
    else :
        c['band'] = band
    c.rename_column('date','time')
    
    errfloor = 0.1 # Set the error floor
    
    for snid in range(1,5):
         
        if not sumdata:
            c['flux']    = c['f%d'%(snid)]
            c['fluxerr'] = np.sqrt((c['f%derr'%(snid)])**2 + (c['f%d'%(snid)]*errfloor)**2+(c['f%d'%(snid)]*lenserr)**2)
            if snid == 4:
                c['fluxerr'] = 2.
            #if (band == 'f814w' and snid == 3):
            #    c['fluxerr'] = c['f%d'%(snid)]*10*errfloor
            c['imageid'] = snid            
            tc = c[('time','band','flux','fluxerr','zp','zpsys','imageid')]
            if phot_d is None:
                phot_d = tc
            else:
                phot_d = vstack([phot_d,tc])
        else:
            flux  += c['f%d'%(snid)]
            ferr2 += c['f%derr'%(snid)]**2 + (c['f%d'%(snid)]*errfloor)**2

    if sumdata:
        c['flux'] = flux
        c['fluxerr'] = np.sqrt(ferr2)
        c['imageid'] = 0
        c = c[('time','band','flux','fluxerr','zp','zpsys','imageid')]
        #sncosmo.write_lc(c,'tables/'+band+'_lc.dat')

        
        if phot_d is None:
            phot_d = c
        else:
            phot_d = vstack([phot_d,c])
# Set redshift, initial t0 and MW extinction
z,lensz,it0 = 0.409,0.216,57653.2

# Galactic dust model, assuming that the photometry has not been
# corrected for MW extinction.
dust = sncosmo.CCM89Dust()

###########################################################################
# LIGHTCURVE MODEL

# setup Hsiao SN Ia stretch-model

p,w,f = sncosmo.read_griddata_fits('/usr/local/lib/python2.7/dist-packages/snpy/typeIa/Hsiao_SED_V3.fits')
mHs = sncosmo.StretchSource(p,w,f,name='hsiao-stretch')
mH = sncosmo.Model(source=mHs,
                  effects=[dust,dust,dust],
                  effect_names=['mw','host','lens'],
                  effect_frames=['obs','rest','free'])
                  
ref = copy.copy(mH)
mH.set(mwr_v=3.1)
mH.set(z=z,lensz=lensz)
mH.set(t0=it0)

##############################################################################
if __name__ == '__main__':
    fit_param = ['s','hostebv','hostr_v',
    't01','t02','t03','t04',
    'lensr_v1','lensr_v2','lensr_v3','lensr_v4',
    'lensebv1','lensebv2','lensebv3','lensebv4',
    'amplitude1','amplitude2','amplitude3','amplitude4']
    
    # the bounds for the fitting parameters
    fit_bounds = {'t01':(it0-5.,it0+5),'t02':(it0-5.,it0+5.),
                  't03':(it0-5.,it0+5),'t04':(it0-5.,it0+5),
                  's':(0.9,1.1),'hostebv':(-0.2,1.0),'hostr_v':(1.0,3.0),
                  'lensr_v1':(1.0,3.0),'lensr_v2':(1.0,3.0),
                  'lensr_v3':(1.0,3.0),'lensr_v4':(1.0,3.0),
                  'lensebv1': (-0.2,1.5), 'lensebv2': (-0.2,1.5), 
                  'lensebv3': (-0.2,1.5), 'lensebv4': (-0.2,1.5),
                  'amplitude1': (1.e-9,1.e-7), 'amplitude2': (1.e-9,1.e-7),
                  'amplitude3': (1.e-9,1.e-7), 'amplitude4': (1.e-10,1.e-8),
                  'lensr_v':(1.0,4.0)} 
    mH.set(amplitude=3.e-8,hostebv=0.1,hostr_v=1.8,lensr_v=2.0)
    MI_mH = MI_model(mH,4)
    model = nest_lc(phot_d,MI_mH,fit_param,fit_bounds)
    myfile = open('nestfitparam.dat','w')
    for name in fit_param:
        myfile.write(name+' = '+str(model.get(name)))
    myfile.close()
    model.plot(phot_d)