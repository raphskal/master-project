# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 09:44:46 2017

@author: kalender
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
import sncosmo
import os
from astropy.io import ascii


class sn16geu(sncosmo.Source):

    _param_names = ['amplitude','s','t0']
    param_names_latex = ['A','stretch','t_0']   # used in plotting display

    def __init__(self,t0, name='sn16geu', version=1.0):
        self.name = name
        self.version = version
        self.fixt0 = t0
        self._parameters = np.array([1.,1.,t0])
        
        self.bands = ['f110w','f160w','f625w','f814w','P60g','P60i','P60r','P48R']
        
        self._model_flux = {}
        
        for name in self.bands:
            c = ascii.read('custom/iPTF16geu_lc_'+name+'_model.dat')
            t = c['time']
            flux = 10**(-0.4*c['magnitude'])
            interp = interp1d(t,flux)
            tgrid = np.linspace(t[0],t[-1],num=100)
            max_time = tgrid[np.argmax(interp(tgrid))]
            tgrid = tgrid - max_time
            self._model_flux[name] = interp1d(tgrid,interp(tgrid+max_time))
            
        print 'Custom source for 16geu was initiated'
            
    def _flux(self,wave):
        # This is a dummy
        return

    
    def bandflux(self, wave, phase,zp,zpsys):
        A,s,t0 = self._parameters
        fluxes = np.empty((len(wave)))
        for i in range(len(fluxes)):
            try:
                fluxes[i] = A*self._model_flux[wave[i].name]((phase[i]-t0)/s)
            except ValueError:
                raise ValueError('Phase is not in sn16geu-model')
            
        return fluxes
                


