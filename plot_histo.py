# -*- coding: utf-8 -*-
# This program calls the Wambsganss inverse ray shooting software.
# It takes the output, visualizes the microlensing pattern, calculates the
# histogram and computes the probability of finding the given mircolensing
"""
Created on Tue Dec  5 13:07:17 2017

@author: kalender
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import time
from astropy.io import fits
from numpy.random import randint
from scipy.interpolate import interp1d

def plot(array,form,value):
    """
    Plots the histrgram for the microlensing pattern
    
    Parameters
    ----------
    array : the 2D array of the microlensing pattern
    form  : amplification or magnification on the axis, choice by user
    """
    n, bins, patches = plt.hist(array.flatten(),100,normed=1)
    interp = interp1d(bins[:-1],n)
    if form == 'amp':
        plt.xlabel('amplification')
    elif form == 'magn':
        plt.xlabel('magnification')
    plt.ylabel('counts/total')
    plt.plot([value,value],[0,1],'r')
    plt.ylim(0,max(n))
    plt.title('Microlensing Histogram for image IV, $\kappa = 0.45,\,\gamma = 0.46$')
    plt.savefig('microlens_histo_'+form+'.png')
    plt.close()
    try:
        print interp(value)
    except ValueError:
        print 0
    
def prob(m,pattern,l):
    """
    Computes the probability of finding the given mircolensing,
    selects random points in the pattern, draws circle and takes 2 random
    points on this circle. Compares their ratio to the given one.
    
    Parameters
    ----------
    m : array consisting of the ratio of the amplifications and error
    pattern : the 2D microlensing pattern
    l :  lxl is the size of the pattern
    runs : # of runs
    counts : counts the matches with the given microlensing
    """
    runs = 1000000
    counts = 0.
    r = 0.25*l
    for i in range(0,runs):
        c1,c2 = randint(int(l*0.25),int(l*0.75)),randint(int(l*0.25),int(l*0.75))
        phi1,phi2 = randint(0,100),randint(0,100)
        m1_r = pattern[int(c1+r*np.cos(phi1/50.*np.pi)),int(c2+r*np.sin(phi1/50.*np.pi))]
        m2_r = pattern[int(c1+r*np.cos(phi2/50.*np.pi)),int(c2+r*np.sin(phi2/50.*np.pi))]
        if m1_r >= m2_r:
            if abs(m1_r/m2_r-m[0])<=m[1]:
                counts += 1.
        else:
            if abs(m2_r/m1_r-m[0])<=m[1]:
                counts += 1.
    print 'The probability of finding this setting is: '+str(100.*counts/float(runs))+'%'
        
        
    
def start(name,form,value):
    """
    Calculates the magnification/amplification of the microlensing pattern. 
    Visualizes the pattern, calls probality and plot function
    
    Parameters
    ----------
    name :  name of the file with the microlensing pattern
    form : magnfication or amplification, choice by user. Ratio in magnification
    not yet calculated. Visualisation of pattern only in magnification
    """
    try:
        data = fits.open(name)
    except IOError:
        print 'Shooting ray could not be initialized with given input.'
        quit()
    pattern = data[0].data
    data.close()
    length = len(pattern)
    real_pattern=np.empty((length,length))
    pixmax = 0.
    pixmin = 10000.
    if form == 'amp':
        for i in range(0,length):
            for j in range(0,length):
                pixmax = max(pixmax,pattern[i][j])
                pixmin = min(pixmin,pattern[i][j])
                real_pattern[i,j] = 10**(0.4*(pattern[i][j]-1024.)/256.)
        pixmax = 10**(0.4*(pixmax-1024.)/256.)
        pixmin = 10**(0.4*(pixmin-1024.)/256.)
        #prob(np.array([3.77/0.17,5.8]),real_pattern,length) 
        plt.imshow(real_pattern,cmap='Spectral')
        plt.colorbar(label='amplification')
        plt.title('Mircolensing pattern in amplification')
        plt.xticks([0.,0.25*length,0.5*length,0.75*length,length],
                    ['-50','-25','0','25','50'])
        plt.yticks([0.,0.25*length,0.5*length,0.75*length,length],
                    ['-50','-25','0','25','50'])
        plt.xlabel('Einstein radii')
        plt.ylabel('Einstein radii')
        plt.savefig('pattern_'+form+'.png')
        plt.close()                    
        plot(real_pattern,form,value)
    elif form == 'magn':
        for i in range(0,length):
            for j in range(0,length):
                pixmax = max(pixmax,pattern[i][j])
                pixmin = min(pixmin,pattern[i][j])
                real_pattern[i,j] = (pattern[i][j]-1024.)/256.
        pixmax = (pixmax-1024.)/256.
        pixmin = (pixmin-1024.)/256.
        #prob(np.array([]]),real_pattern,length)
        plt.imshow(real_pattern,cmap='Spectral')
        plt.colorbar(label='magnification')
        plt.title('Mircolensing pattern in Magnitudes')
        plt.xticks([0.,0.25*length,0.5*length,0.75*length,length],
                    ['-2','-1','0','1','2'])
        plt.yticks([0.,0.25*length,0.5*length,0.75*length,length],
                    ['-2','-1','0','1','2'])
        plt.xlabel('Einstein radii')
        plt.ylabel('Einstein radii')
        plt.savefig('pattern_'+form+'.png')
        plt.close()                   
        plot(real_pattern,form,value)
    else:
        print 'form must be amp or magn'
        
    
    
    
if __name__ == '__main__':
    # Modify the 'input' file in this folder. Then run the script. Choose
    # 'amp' or 'magn' for the form. Make sure name of the file is correct.
    os.popen('./microlens')
    time.sleep(2)
    start('IRIS401.fits','magn',3.01)
    os.remove('IRIS401.fits')
    
    
    
