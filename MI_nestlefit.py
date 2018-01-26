# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 09:58:58 2018

@author: kalender
"""
from __future__ import absolute_import
import copy
import time
import matplotlib
matplotlib.use('Agg')
from collections import OrderedDict
import warnings
from sncosmo.photdata import photometric_data
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
from astropy.extern import six
from mpl_toolkits.axes_grid1 import make_axes_locatable


from sncosmo.utils import Result, Interp1D, ppf

class MI_model(object):
    """
    Constructs a model that consists of n-sncosmo models in order to represent 
    a n-imaged supernova. This is used for the nestlefit.
    """
    def __init__(self, model, n):
        """
        Initializes the MI_model
    
        Parameters
        ----------
        model: the sncosmo-model that is used for the multiple images
        n: number of images
    
        returns: MI_model containing n sncosmo-models
        """
        self.nimg = n
        model = copy.copy(model)
        self.models = []
        
        for i in range(n):
            self.models.append(copy.copy(model))
        
        self.param_names = self._get_parameters()
            
    def _get_parameters(self):
        """
        Gets the parameter names of the MI_model
        
        Parameters
        ----------
        self : this model
        
        returns : all parameters used in this model
        """
        params = ['hostr_v','hostebv','s']
        for i in range(self.nimg):
            params.append('amplitude'+str(i+1))
            params.append('t0'+str(i+1))
            params.append('lensr_v'+str(i+1))
            params.append('lensebv'+str(i+1))
        return params
        
    def _get_values(self):
        """
        Gets the parameter values of the MI_model
        
        Parameters
        ----------
        self : this model
        
        returns : all parameter values used in this model
        """
        values = []
        for name in self.param_names:
            values.append(self.get(name))
        return values
            
    def get(self, name):
        """
        Gets the value of a specific parameter
        
        Parameters
        ----------
        self : this model
        name :  parameter that value is required
        
        returns :  value parameter value 'name'
        """
        if name[0] == 'h' or name == 's':
            return self.models[0].get(name)
        else:
            return self.models[int(name[-1])-1].get(name[:-1])
            
    def set(self, name, value):
        """
        Sets the value of a specific parameter
        
        Parameters
        ----------
        self : this model
        name :  parameter value that is to be changed
        
        """
        if name[0] == 'h' or name == 's':
            for i in range(self.nimg):
                self.models[i].update({name:value})
        else:
            self.models[int(name[-1])-1].update({name[0:-1]:value})
            
    def plot(self,data, zp=25., zpsys='ab', ncol=2):
        """
        Plots the data for all bands
        
        Parameters
        ----------
        self : this model
        data : photometric data
        
        returns : Visualization of the photometric data and this model
        in 16geu_nestfit.jpg
        """
        
        all_bands = list(set(data['band']))
        nbands = len(all_bands)
        nrow = int(np.ceil(nbands/float(ncol)))
        
        tmin,tmax = np.min(data['time']),np.max(data['time'])
        fig = plt.figure(figsize=(10,12))
        
        ground = ['P60g','P48R','P60i','P60r']
    


        #timegrid = np.linspace(model.mintime(), model.maxtime(),
        #                       int(model.maxtime() - model.mintime() + 1))
        #totmag=0.
        #mag_array,i = np.zeros(len(ground)),0
        tgrid = np.linspace(tmin,tmax,300)
        marker = ['.','o','s','d','^']
        color  = ['k','c','m','g','b']
        for n,band in enumerate(all_bands):
            mask = data['band'] == band
            
            ax = fig.add_subplot(nrow,ncol,n+1)
            ax.axhline(0.,color='k',ls=':')
            #plt.axhline(0.,color='k',ls=':')            

            sum_model_flux,plot_sum,diff,sigma = None,False,None,None
            for n in range(self.nimg+1):
                imask = data['imageid'] == n                
                m = mask*imask

                if n > 0:
                    model_flux = self.models[n-1].bandflux(band, tgrid, zp=zp, zpsys=zpsys)
                    #if band in ground:
                        #varmag = self.cleanmax(band,tgrid,zp,zpsys)*magn[n-1]*lensing[n-1]/self.physmax(band,tgrid,zp,zpsys)
                        #totmag += varmag
                        #mag_array[i] += varmag
                    if sum_model_flux is None:
                        sum_model_flux = model_flux
                    else:
                        sum_model_flux += model_flux
                
                if not np.any(m) : continue
                        
                d = photometric_data(data[m])
                d = d.normalized(zp=zp, zpsys=zpsys)
                
                    
                if not np.all(np.ediff1d(d.time) >= 0.0):
                    sortidx = np.argsort(d.time)
                    d = d[sortidx]
                else:
                    sortidx = None
                                
                ax.errorbar(d.time,d.flux,yerr=d.fluxerr,marker=marker[n],capsize=0.,ls='',color=color[n])
                
                if n > 0:   
                    myfile = open('tables/16geufit_'+band+'_nimg'+str(n)+'.dat','w')
                    myfile.write('date flux\n')
                    for l in range(len(tgrid)):
                        myfile.write(str(tgrid[l])+' '+str(model_flux[l])+'\n')
                    myfile.close()
                    ax.plot(tgrid,model_flux,color=color[n],ls='-')
                    var = interp1d(tgrid,model_flux)
                    if diff is None:
                        diff = (d.flux-var(d.time))/d.fluxerr
                        sigma = np.sum(((d.flux-var(d.time))/d.fluxerr)**2)
                    else:
                        diff += (d.flux-var(d.time))/d.fluxerr
                        sigma += np.sum(((d.flux-var(d.time))/d.fluxerr)**2)
                else:
                    plot_sum = True
                    
            if (sigma is not None) and not (self.nimg*len(d.time) == 1):
                sigma=sigma/(self.nimg*len(d.time)-1.)
            
            if plot_sum:
                myfile = open('tables/16geufit_'+band+'.dat','w')
                myfile.write('date flux\n')
                for l in range(len(tgrid)):
                    myfile.write(str(tgrid[l])+' '+str(model_flux[l])+'\n')
                myfile.close()
                var = interp1d(tgrid,sum_model_flux)
                diff = (d.flux-var(d.time))/d.fluxerr
                if not (len(d.time) == 1):
                    sigma = np.sum(((d.flux-var(d.time))/d.fluxerr)**2)/(len(d.time)-1.)
                else:
                    sigma = np.sum(((d.flux-var(d.time))/d.fluxerr)**2)/(len(d.time))
                ax.plot(tgrid,sum_model_flux,color=color[0],ls='-')
            
            #if band in ground: 
            #    i += 1
                
            ax.set_title(band)
            ax.set_xlim((tmin,tmax))
            ax.set_ylabel('flux')
            ca = fig.gca().get_xticks()
            ax.set_xticks(ca[1:-1])
            divider = make_axes_locatable(ax)
            sigma = np.sqrt(sigma)
            axpulls = divider.append_axes('bottom', size='30%', pad=0.0,
                                          sharex=ax)
            axpulls.plot(d.time, diff, 'x')
            axpulls.set_ylim(-4.*sigma,4.*sigma)
            axpulls.set_yticks([-3*sigma,-1*sigma,1*sigma,3*sigma])
            axpulls.set_yticklabels(['$-3\sigma$','$-1\sigma$','$+1\sigma$','$+3\sigma$'])
            axpulls.fill_between(tgrid,-1*sigma,1*sigma,facecolor='green',alpha=0.3)
            axpulls.set_xlabel('Julian Date')                      
            
        
        fig.tight_layout()
        fig.savefig('16geu_nestfit.png')
        plt.close()
        
        #mean =  totmag/len(ground)
        #sigmag = np.sqrt(np.sum((mag_array-mean)**2)/(len(ground)-1))
        
        

    def chisq(self,data):
            """
            Calculates the chisquare of this model and the data. imid (=imageid)
            is used to distinguish between ground (=0) where the sum of the fluxes 
            matter and space (=1,..,4) where we can use the individual fluxes
        
            Parameters
            -----------
            self : this model
            data : photometric data
        
            returns chisquare of data and model
            """
            chi2 = 0.
            for imid in range(self.nimg+1):
                mask = data['imageid'] == imid
                if not np.any(mask): continue
                
                d = photometric_data(data[mask])
                
                # if the image-id is zero it should be compared to the sum
                # of the flux from all images
                if imid > 0:
                    m = self.models[imid-1]
                    model_flux = m.bandflux(d.band, d.time, zp=d.zp, zpsys=d.zpsys)
                else:
                    model_flux = np.zeros(len(d.time))
                    for m in self.models:
                        mf = m.bandflux(d.band, d.time, zp=d.zp, zpsys=d.zpsys)
                        model_flux += mf
                    
                if not np.all(d.fluxerr > 0.):
                    print d.fluxerr
                    print "OOPS"
                    
                diff = d.flux - model_flux
                
                # is this going to work when we are masking?
                cov = np.diag(d.fluxerr**2) if d.fluxcov is None else d.fluxcov            
                invcov = np.linalg.inv(cov)            
                chi2 += np.dot(np.dot(diff, invcov), diff)
            #print 'chi='+str(chi2)            
            return chi2    

def nest_lc(data, model, vparam_names, bounds, guess_amplitude_bound=False,
            minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
            maxiter=None, maxcall=None, modelcov=False, rstate=None,
            verbose=False, warn=True, **kwargs):
    """Run nested sampling algorithm to estimate model parameters and evidence.
    Parameters
    ----------
    data : `~astropy.table.Table` or `~numpy.ndarray` or `dict`
        Table of photometric data. Must include certain columns.
        See the "Photometric Data" section of the documentation for
        required columns.
    model : `~sncosmo.Model`
        The model to fit.
    vparam_names : list
        Model parameters to vary in the fit.
    bounds : `dict`
        Bounded range for each parameter. Bounds must be given for
        each parameter, with the exception of ``t0``: by default, the
        minimum bound is such that the latest phase of the model lines
        up with the earliest data point and the maximum bound is such
        that the earliest phase of the model lines up with the latest
        data point.
    guess_amplitude_bound : bool, optional
        If true, bounds for the model's amplitude parameter are determined
        automatically based on the data and do not need to be included in
        `bounds`. The lower limit is set to zero and the upper limit is 10
        times the amplitude "guess" (which is based on the highest-flux
        data point in any band). Default is False.
    minsnr : float, optional
        Minimum signal-to-noise ratio of data points to use when guessing
        amplitude bound. Default is 5.
    priors : `dict`, optional
        Prior probability distribution function for each parameter. The keys
        should be parameter names and the values should be callables that
        accept a float. If a parameter is not in the dictionary, the prior
        defaults to a flat distribution between the bounds.
    ppfs : `dict`, optional
        Prior percent point function (inverse of the cumulative distribution
        function) for each parameter. If a parameter is in this dictionary,
        the ppf takes precedence over a prior pdf specified in ``priors``.
    npoints : int, optional
        Number of active samples to use. Increasing this value increases
        the accuracy (due to denser sampling) and also the time
        to solution.
    method : {'classic', 'single', 'multi'}, optional
        Method used to select new points. Choices are 'classic',
        single-ellipsoidal ('single'), multi-ellipsoidal ('multi'). Default
        is 'single'.
    maxiter : int, optional
        Maximum number of iterations. Iteration may stop earlier if
        termination condition is reached. Default is no limit.
    maxcall : int, optional
        Maximum number of likelihood evaluations. Iteration may stop earlier
        if termination condition is reached. Default is no limit.
    modelcov : bool, optional
        Include model covariance when calculating chisq. Default is False.
    rstate : `~numpy.random.RandomState`, optional
        RandomState instance. If not given, the global random state of the
        ``numpy.random`` module will be used.
    verbose : bool, optional
        Print running evidence sum on a single line.
    warn : bool, optional
        Issue warning when dropping bands outside the model range. Default is
        True.
        *New in version 1.5.0*
    Returns
    -------
    res : Result
        Attributes are:
        * ``niter``: total number of iterations
        * ``ncall``: total number of likelihood function calls
        * ``time``: time in seconds spent in iteration loop.
        * ``logz``: natural log of the Bayesian evidence Z.
        * ``logzerr``: estimate of uncertainty in logz (due to finite sampling)
        * ``h``: Bayesian information.
        * ``vparam_names``: list of parameter names varied.
        * ``samples``: 2-d `~numpy.ndarray`, shape is (nsamples, nparameters).
          Each row is the parameter values for a single sample. For example,
          ``samples[0, :]`` is the parameter values for the first sample.
        * ``logprior``: 1-d `~numpy.ndarray` (length=nsamples);
          log(prior volume) for each sample.
        * ``logl``: 1-d `~numpy.ndarray` (length=nsamples); log(likelihood)
          for each sample.
        * ``weights``: 1-d `~numpy.ndarray` (length=nsamples);
          Weight corresponding to each sample. The weight is proportional to
          the prior * likelihood for the sample.
        * ``parameters``: 1-d `~numpy.ndarray` of weighted-mean parameter
          values from samples (including fixed parameters). Order corresponds
          to ``model.param_names``.
        * ``covariance``: 2-d `~numpy.ndarray` of parameter covariance;
          indicies correspond to order of ``vparam_names``. Calculated from
          ``samples`` and ``weights``.
        * ``errors``: OrderedDict of varied parameter uncertainties.
          Corresponds to square root of diagonal entries in covariance matrix.
        * ``ndof``: Number of degrees of freedom (len(data) -
          len(vparam_names)).
        * ``bounds``: Dictionary of bounds on varied parameters (including
          any automatically determined bounds).
        * ``data_mask``: Boolean array the same length as data specifying
          whether each observation was used.
          *New in version 1.5.0.*
    estimated_model : `~sncosmo.Model`
        A copy of the model with parameters set to the values in
        ``res.parameters``.
    """

    try:
        import nestle
    except ImportError:
        raise ImportError("nest_lc() requires the nestle package.")

    # warnings
    if "nobj" in kwargs:
        warnings.warn("The nobj keyword is deprecated and will be removed in "
                      "sncosmo v2.0. Use `npoints` instead.")
        npoints = kwargs.pop("nobj")

    # experimental parameters
    tied = kwargs.get("tied", None)
    """
    # sort by time
    if not np.all(np.ediff1d(data.time) >= 0.0):
        sortidx = np.argsort(data.time)
        data = data[sortidx]
    else:
        sortidx = None
    """
    model = copy.copy(model)
    bounds = copy.copy(bounds)  # need to copy this b/c we modify it below

    # Order vparam_names the same way it is ordered in the model:
    vparam_names = [s for s in model.param_names if s in vparam_names]
    print vparam_names
    """
    # Drop data that the model doesn't cover.
    fitdata, data_mask = cut_bands(data, model,
                                   z_bounds=bounds.get('z', None),
                                   warn=warn)
    

    if guess_amplitude_bound:
        if model.param_names[2] not in vparam_names:
            raise ValueError("Amplitude bounds guessing enabled but "
                             "amplitude parameter {0!r} is not varied"
                             .format(model.param_names[2]))
        if model.param_names[2] in bounds:
            raise ValueError("cannot supply bounds for parameter {0!r}"
                             " when guess_amplitude_bound=True"
                             .format(model.param_names[2]))

        # If redshift is bounded, set model redshift to midpoint of bounds
        # when doing the guess.
        if 'z' in bounds:
            model.set(z=sum(bounds['z']) / 2.)
        _, amplitude = guess_t0_and_amplitude(fitdata, model, minsnr)
        bounds[model.param_names[2]] = (0., 10. * amplitude)

    # Find t0 bounds to use, if not explicitly given
    if 't0' in vparam_names and 't0' not in bounds:
        bounds['t0'] = t0_bounds(data, model)
    """
    if ppfs is None:
        ppfs = {}
    if tied is None:
        tied = {}

    # Convert bounds/priors combinations into ppfs
    if bounds is not None:
        for key, val in six.iteritems(bounds):
            if key in ppfs:
                continue  # ppfs take priority over bounds/priors
            a, b = val
            if priors is not None and key in priors:
                # solve ppf at discrete points and return interpolating
                # function
                x_samples = np.linspace(0., 1., 101)
                ppf_samples = ppf(priors[key], x_samples, a, b)
                f = Interp1D(0., 1., ppf_samples)
            else:
                f = Interp1D(0., 1., np.array([a, b]))
            ppfs[key] = f

    # NOTE: It is important that iparam_names is in the same order
    # every time, otherwise results will not be reproducible, even
    # with same random seed.  This is because iparam_names[i] is
    # matched to u[i] below and u will be in a reproducible order,
    # so iparam_names must also be.
    iparam_names = [key for key in vparam_names if key in ppfs]
    ppflist = [ppfs[key] for key in iparam_names]
    npdim = len(iparam_names)  # length of u
    ndim = len(vparam_names)  # length of v

    # Check that all param_names either have a direct prior or are tied.
    for name in vparam_names:
        if name in iparam_names:
            continue
        if name in tied:
            continue
        raise ValueError("Must supply ppf or bounds or tied for parameter '{}'"
                         .format(name))

    def prior_transform(u):
        d = {}
        for i in range(npdim):
            d[iparam_names[i]] = ppflist[i](u[i])
        v = np.empty(ndim, dtype=np.float)
        for i in range(ndim):
            key = vparam_names[i]
            if key in d:
                v[i] = d[key]
            else:
                v[i] = tied[key](d)
        return v

    # Indicies of the model parameters in vparam_names
    #idx = np.array([model.param_names.index(name) for name in vparam_names])

    def loglike(parameters):
        for i, name in enumerate(vparam_names):
            model.set(name,parameters[i])
        #model.parameters[idx] = parameters
        return -0.5 * model.chisq(data)

    t0 = time.time()
    res = nestle.sample(loglike, prior_transform, ndim, npdim=npdim,
                        npoints=npoints, method=method, maxiter=maxiter,
                        maxcall=maxcall, rstate=rstate,
                        callback=(nestle.print_progress if verbose else None))
    """
    # estimate parameters and covariance from samples
    vparameters, cov = nestle.mean_and_cov(res.samples, res.weights)
    # update model parameters to estimated ones.
    for i, name in enumerate(vparam_names):
        model.set(name,vparameters[i])
    
    # If we need to, unsort the mask so mask applies to input data
    if sortidx is not None: 
        unsort_idx = np.argsort(sortidx)  # indicies that will unsort array
        data_mask = data_mask[unsort_idx]
    
    # `res` is a nestle.Result object. Collect result into a sncosmo.Result
    # object for consistency, and add more fields.
    res = Result(niter=res.niter,
                 ncall=res.ncall,
                 logz=res.logz,
                 logzerr=res.logzerr,
                 h=res.h,
                 samples=res.samples,
                 weights=res.weights,
                 logvol=res.logvol,
                 logl=res.logl,
                 vparam_names=copy.copy(vparam_names),
                 ndof=len(data) - len(vparam_names),
                 bounds=bounds,
                 time=elapsed,
                 parameters=model._get_values(),
                 covariance=cov,
                 errors=OrderedDict(zip(vparam_names,
                                        np.sqrt(np.diagonal(cov)))),
                 param_dict=OrderedDict(zip(model.param_names,
                                            model._get_values())),
               )
    
    # Deprecated result fields.
    depmsg = ("The `param_names` attribute is deprecated in sncosmo v1.0 "
              "and will be removed in sncosmo v2.0."
              "Use `vparam_names` instead.")
    #res.__dict__['deprecated']['param_names'] = (res.vparam_names, depmsg)

    depmsg = ("The `logprior` attribute is deprecated in sncosmo v1.2 "
              "and will be changed in sncosmo v2.0."
              "Use `logvol` instead.")
    #res.__dict__['deprecated']['logprior'] = (res.logvol, depmsg)
    """
    return model