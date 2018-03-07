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
errfloor = 0.0
##############################################################################
# Here the error bars of the grond data can be adjusted. 
lenserr = 0.00  # set to 0.05 when using gravlens model
#phot_d['fluxerr'] = np.sqrt(phot_d['fluxerr']**2+(phot_d['flux']*lenserr)**2)

phot_d['fluxerr'] = np.sqrt((phot_d['flux']*errfloor)**2+(phot_d['fluxerr'])**2)   # scale error bars
errfloor = 0.08
##############################################################################
# Set up the data of the HST to be read in
sumdata = False
lcfname = { # the names of the data files
    'uvf625w' : 'lc_uvf625w_resolved.csv',
    'uvf814w' : 'lc_uvf814w_resolved.csv',
    'f110w' : 'lc_f110w_resolved.csv',
    'f160w' : 'lc_f160w_resolved.csv',
    'f390w' : 'lc_f390w_resolved.csv',
    'f475w' : 'lc_f475w_resolved.csv'
}
# These are the bands that are cosidered, use these for hsiao-stretch!
for band in ['uvf625w','uvf814w','f390w','f475w','f160w','f110w'] :
# And these for the custom sn16geu source
#for n,band in enumerate(['f625w','f814w','f160w']) :  
    fname = lcfname[band]                                    
    c = ascii.read(os.path.join('photometry/lightcurve',fname))
    flux,ferr2 = np.zeros(len(c)),np.zeros(len(c)) 
    c['fluxerr'] = np.sqrt((c['fluxerr'])**2 + (c['flux']*errfloor)**2+(5)**2)
    tc = c[('time','band','flux','fluxerr','zp','zpsys','imageid')]
    if phot_d is None:
        phot_d = tc
    else:
        phot_d = vstack([phot_d,tc])
        
phot_d.sort('time')
            
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
mH.set(amplitude=3.e-8,hostebv=0.1,hostr_v=3.1,lensr_v=2.0)


# Setup model based on SNoopy fit for 16geu 
#geu_source = sn16geu(it0,1.e+9)

##############################################################################
if __name__ == '__main__':
    # The parameters that are varied in the fit
    fit_param = ['s','t01','t02','t03','t04',
    'amplitude1','amplitude2','amplitude3','amplitude4',
    'hostebv','lensr_v',
    'f',
    'lensebv2','lensebv3','lensebv4'
    ]
    
    # With the custom sn16geu-model, we are unable to fit for the extinction    
    fit_param_16geu = ['s',
    't01','t02','t03','t04',
    'amplitude1','amplitude2','amplitude3','amplitude4']
    
    # the bounds for the fitting parameters
    fit_bounds = {'t01':(it0-5.,it0+5),'t02':(-3,+3),
                  't03':(-3.,+3),'t04':(-3,+3),
                  's':(0.8,1.2),'hostebv':(0.0,1.0),'hostr_v':(1.0,4.0),
                  'lensebv1': (0.0,1.0), 'lensebv2': (-0.2,1.0), 
                  'lensebv3': (0.5,1.5), 'lensebv4': (-0.2,1.0),
                  #'amplitude1': (1.e+7,1.e+10), 'amplitude2': (1.e+7,1.e+10),
                  #'amplitude3': (1.e+7,1.e+10), 'amplitude4': (1.e+6,1.e+9)}
                  'amplitude1': (1.e-9,1.e-7), 'amplitude2': (1.e-9,1.e-7),
                  'amplitude3': (1.e-9,1.e-7), 'amplitude4': (1.e-10,1.e-8),
                  'lensr_v':(2.0,4.0), 'f':(1.4,2.0), 'lensr_v1':(1.0,3.0)} 
 
    # Start the nestlefit with the hsiao-stretch or the custom_model
    # Save the parameters in 'nestfitparam.dat' and plot the model

    MI_mH = MI_model(mH,4,f=1.7,useprior=True, samerv=True)
    model, res = nest_lc(phot_d,MI_mH,fit_param,fit_bounds)
    #MI_geu_source = MI_model(geu_source,4)
    #model,res = nest_lc(phot_d,MI_geu_source,fit_param_16geu,fit_bounds)
    
    # Draw the correlations between the parameters that are degenerated
    samp, params  =  res.samples, res.vparam_names
    samples(samp,params,fit_param[9:])

    # Draw the probability density function
    df,data = res.ndof, model.tot_amp
    fmin,minamp = min(model.chi), data[np.argmin(model.chi)]
    histogram(data)
    
    # Plot the lightcurves for the model and write the parameters into a file
    model.plot(phot_d)
    myfile = open('nestfitparam.dat','w')
    for name in fit_param:
        myfile.write(name+' = '+str(model.get(name))+' ('+str(res.errors[name])+')\n')
    myfile.write("chiqsquare = "+str(fmin)+' Dof = ' +str(df)+'\n')
    myfile.write('total amplification = '+str(minamp)+'\n')
    myfile.write('Bayesian evidence z = '+str(np.exp(res.logz)))
    myfile.close()
    
