# -*- coding: utf-8 -*-
"""
This Program reads the fitted data from files, masks it and starts the fitting
in MI_nestlefit.py


Created on Tue Jan 23 10:40:53 2018

@author: kalender
"""
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from custom_source import sn16geu
from astropy.table import vstack
from astropy.io import ascii
from MI_nestlefit import MI_model, nest_lc
from scipy.interpolate import interp1d
from stats import histogram,samples
import sncosmo
import os
import corner
import copy
import filters

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
# These are the bands that are cosidered, use these for hsiao-stretch!
for n,band in enumerate(['f625w','f814w','f390w','f475w','f160w']) :
# And these for the custom sn16geu source
#for n,band in enumerate(['f625w','f814w','f160w']) :
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
z,lensz,it0 = 0.409,0.216,57651.2

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
# Set the initial model, construct a MI_model with it
mH.set(amplitude=3.e-8,hostebv=0.1,hostr_v=1.4,lensr_v=3.1)


# Setup model based on SNoopy fit for 16geu 
#geu_source = sn16geu(it0,1.e+9)

##############################################################################
if __name__ == '__main__':
    # The parameters that are varied in the fit
    fit_param = ['s','hostebv','hostr_v',
    't01','t02','t03','t04',
    'f',
    'lensebv2','lensebv3','lensebv4',
    'amplitude1','amplitude2','amplitude3','amplitude4']
    
    # With the custom sn16geu-model, we are unable to fit for the extinction    
    fit_param_16geu = ['s',
    't01','t02','t03','t04',
    'amplitude1','amplitude2','amplitude3','amplitude4']
    
    # the bounds for the fitting parameters
    fit_bounds = {'t01':(it0-2.,it0+2),'t02':(-3,+3),
                  't03':(-3.,+3),'t04':(-3,+3),
                  's':(0.95,1.05),'hostebv':(0.0,1.0),'hostr_v':(1.0,3.0),
                  'lensebv1': (0.0,1.0), 'lensebv2': (-0.2,1.0), 
                  'lensebv3': (-0.2,1.0), 'lensebv4': (-0.2,1.0),
                  #'amplitude1': (1.e+7,1.e+10), 'amplitude2': (1.e+7,1.e+10),
                  #'amplitude3': (1.e+7,1.e+10), 'amplitude4': (1.e+6,1.e+9)}
                  'amplitude1': (1.e-9,1.e-7), 'amplitude2': (1.e-9,1.e-7),
                  'amplitude3': (1.e-9,1.e-7), 'amplitude4': (1.e-10,1.e-8),
                  'lensr_v':(1.0,4.0), 'f':(1.4,2.0), 'lensr_v1':(1.0,3.0)} 
 
    # Start the nestlefit with the hsiao-stretch or the custom_model
    # Save the parameters in 'nestfitparam.dat' and plot the model

    MI_mH = MI_model(mH,4,f=1.7,useprior=True, samerv=True)
    model, res = nest_lc(phot_d,MI_mH,fit_param,fit_bounds)
    #MI_geu_source = MI_model(geu_source,4)
    #model,res = nest_lc(phot_d,MI_geu_source,fit_param_16geu,fit_bounds)
    
    # Draw the correlations between the parameters that are degenerated
    samp, params  =  res.samples, res.vparam_names
    samples(samp,params,['hostebv','hostr_v','lensebv2','f','lensebv3','lensebv4'])

    # Draw the probability density function
    df,data = res.ndof, model.tot_amp
    fmin,minamp = min(model.chi), data[np.argmin(model.chi)]
    histogram(data)
    
    # Plot the lightcurves for the model and write the parameters into a file
    model.plot(phot_d)
    myfile = open('nestfitparam.dat','w')
    for name in fit_param:
        myfile.write(name+' = '+str(model.get(name))+' ('+str(res.errors[name])+')\n')
    myfile.write("chiqsquare = "+str(fmin)+'\n')
    myfile.write('total amplification = '+str(minamp))
    myfile.close()
    
