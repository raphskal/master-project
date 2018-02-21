# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 10:19:21 2018

@author: kalender
"""
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
from scipy.interpolate import interp1d
from scipy.integrate import quad

def samples(samples,params,names):
    table = []
    for name in names:
        table.append(samples[:,params.index(name)])
    table = np.array(table).T
    corner.corner(table, smooth=1, labels=names)
    plt.savefig('samples.png')
    plt.close()


def histogram(data):
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(data,bins=100, normed=1,zorder=1)
    myfile = open('plots/amp_histo.dat','w')
    myfile.write('bins  normalized counts\n')
    for i in range(0,len(bins)-1):
        myfile.write(str(bins[i+1])+' '+str(n[i])+'\n')
    myfile.close()
    interp,grid = interp1d(bins[1:],n),np.linspace(bins[1],bins[-1],num=1e+3)
    i,area = 1,0
    while area < 0.159:
        area = quad(interp,grid[0],grid[i])[0]
        i += 1
    a1,n_a1 = grid[i],interp(grid[i])
    i,area = len(grid)-2,0
    while area < 0.159:
        area = quad(interp,grid[i],grid[-1])[0]
        i-=1
    b1,n_b1 = grid[i],interp(grid[i])
    i,area = 1,0
    while area < 0.023:
        area = quad(interp,grid[0],grid[i])[0]
        i += 1
    a2,n_a2 = grid[i],interp(grid[i])
    i,area = len(grid)-2,0
    while area < 0.023:
        area = quad(interp,grid[i],grid[-1])[0]
        i -= 1
    b2,n_b2 = grid[i],interp(grid[i])
    i,area = 1,0
    while area < 0.002:
        area = quad(interp,grid[0],grid[i])[0]
        i += 1
    a3,n_a3 = grid[i],interp(grid[i])
    i,area = len(grid)-2,0
    while area < 0.002:
        area = quad(interp,grid[i],grid[-1])[0]
        i-=1
    b3,n_b3 = grid[i],interp(grid[i])
    ax.fill_between(grid,0.,interp(grid),where=grid<=a1,alpha=0.75,facecolor='green',zorder=2)
    plt.text(a1-5,n_a1,round(a1,1),fontsize=16)
    ax.fill_between(grid,0.,interp(grid),where=grid>=b1,alpha=0.75,facecolor='green',zorder=2)
    plt.text(b1-5,n_b1,round(b1,1),fontsize=16)
    ax.fill_between(grid,0.,interp(grid),where=grid<=a2,alpha=0.75,facecolor='yellow',zorder=3)
    plt.text(a2-5,n_a2,round(a2,1),fontsize=16)
    ax.fill_between(grid,0.,interp(grid),where=grid>=b2,alpha=0.75,facecolor='yellow',zorder=3)
    plt.text(b2-5,n_b2,round(b2,1),fontsize=16)
    ax.fill_between(grid,0.,interp(grid),where=grid<=a3,alpha=0.75,facecolor='red',zorder=4)
    plt.text(a3-5,n_a3,round(a3,1),fontsize=16)
    ax.fill_between(grid,0.,interp(grid),where=grid>=b3,alpha=0.75,facecolor='red',zorder=4)
    plt.text(b3-5,n_b3,round(b3,1),fontsize=16)
    ax.set_title('Histogram of total amplification')
    ax.set_xlabel('total amplification')
    ax.set_ylabel('counts/total')
    fig.savefig('totamp_histo.png')
    plt.close()