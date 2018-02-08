# -*- coding: utf-8 -*-
"""
This prorgram creates a hard-coded template based on fits to the 16geu data
with SNooPy. It can be used in call_MIfit.py for the fitting process.


Created on Mon Nov  6 09:44:46 2017

@author: kalender
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
import sncosmo
import os
import matplotlib.pyplot as plt
from astropy.io import ascii


class sn16geu(sncosmo.Source):

    _param_names = ['amplitude','s','t0']
    param_names_latex = ['A','stretch','t_0']   # used in plotting display

    def __init__(self,t0,amplitude, name='sn16geu', version=1.0):
        self.name = name
        self.version = version
        self.fixt0 = t0
        self._parameters = np.array([amplitude,1.,t0])
        
        self.bands = ['f110w','f160w','uvf625w','uvf814w','P60g','P60i','P60r','P48R']
        
        self._model_flux = {}
        
        for name in self.bands:
            c = ascii.read('custom/model/iPTF16geu_lc_'+name+'_model.dat')
            t = c['time']
            flux = 10**(-0.4*c['magnitude'])
            interp = interp1d(t,flux)
            tgrid = np.linspace(t[0],t[-1],num=100)
            flux = interp(tgrid)
            #max_time = tgrid[np.argmax(interp(tgrid))]
            plt.plot(tgrid,flux)
            plt.savefig('custom/plots/'+name+'_lc.png')
            plt.close()
            self._model_flux[name] = [tgrid,flux]
            
        print 'Custom source for 16geu was initiated'
            
    def _flux(self,wave):
        # This is a dummy
        return

    
    def bandflux(self, wave, phase,zp,zpsys):
        A,s,t0 = self._parameters
        fluxes = np.empty((len(wave))) 
        if type(wave) is np.string_:
                interp = interp1d(self._model_flux[wave][0]-t0,self._model_flux[wave][1])
                fluxes = A*interp((phase-t0)/s)
        else:
            for i in range(len(fluxes)):
                if type(wave[i]) is str:
                    interp = interp1d(self._model_flux[wave[i]][0]-t0,self._model_flux[wave[i]][1])
                    fluxes = A*interp((phase[i]-t0)/s)
                else:
                    try:
                        interp = interp1d(self._model_flux[wave[i].name][0]-t0,self._model_flux[wave[i].name][1])
                        fluxes[i] = A*interp((phase[i]-t0)/s)
                    except ValueError:
                        raise ValueError('Phase = '+str(phase[i]-t0)+' for wave = '+str(wave[i].name)+' is not in sn16geu-model')
                
        return fluxes
                


