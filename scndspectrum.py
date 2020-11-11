#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 10:50:45 2020

@author: fzeuner

This code will plot the solar intensity spectrum in color 
and the second solar spectrum at grey scale.

"""
################### 
# -----------------
###################

import numpy as np
import pylab as p
from matplotlib import gridspec
from matplotlib.ticker import NullFormatter
import matplotlib as mp

def scndspectrum(file_name,save_file_name,starto,endo,lines,xsize,ysize,zoom):


    # Style
    
    params = {'backend': 'pdf',
              'axes.linewidth': zoom*0.5,
              'axes.edgecolor': 'grey',
              'font.family': 'serif',
              'lines.dash_joinstyle': 'round',
              'lines.dash_capstyle': 'round',
              'lines.solid_joinstyle': 'round',
              'lines.solid_capstyle': 'round',
              'lines.markersize': 0.5*zoom,
              'lines.linewidth': 0.6*zoom,
              }
    p.rcParams.update(params)
    p.style.use('dark_background')
    
    ################### 
    # Definitions
    ###################
    
    # This color scale comes from 
    # https://stackoverflow.com/questions/53658296/mandelbrot-set-color-spectrum-suggestions
    
    def spectral_color(l): #// RGB <0,1> <- lambda l <400,700> [nm]
    
        t=0.0 
        r=0.0
        g=0.0
        b=0.0
        if ((l>=400.0) and (l<410.0)):
            t=(l-400.0)/(410.0-400.0)
            r= +(0.33*t)-(0.20*t*t)
        elif ((l>=410.0) and (l<475.0)):
            t=(l-410.0)/(475.0-410.0)
            r=0.14         -(0.13*t*t)
        elif ((l>=545.0) and (l<595.0)):
           t=(l-545.0)/(595.0-545.0)
           r=    +(1.98*t)-(     t*t)
        elif ((l>=595.0) and (l<650.0)):
          t=(l-595.0)/(650.0-595.0)
          r=0.98+(0.06*t)-(0.40*t*t)
        elif ((l>=650.0) and (l<700.0)):
            t=(l-650.0)/(700.0-650.0)
            r=0.65-(0.84*t)+(0.20*t*t)
        if ((l>=415.0) and (l<475.0)):
            t=(l-415.0)/(475.0-415.0)
            g=+(0.80*t*t)
        elif ((l>=475.0) and (l<590.0)):
            t=(l-475.0)/(590.0-475.0)
            g=0.8 +(0.76*t)-(0.80*t*t)
        elif ((l>=585.0) and (l<639.0)):
            t=(l-585.0)/(639.0-585.0)
            g=0.84-(0.84*t)           
        if ((l>=400.0) and (l<475.0)):
            t=(l-400.0)/(475.0-400.0)
            b=+(2.20*t)-(1.50*t*t)
        elif ((l>=475.0) and (l<560.0)):
            t=(l-475.0)/(560.0-475.0)
            b=0.7 -(t)+(0.30*t*t)
        return np.array([r,g,b])
    
    
    
    
    cdict = {'red':   ((0.0, 0.0, 0.0), (0.25, 0.40, 0.3), (0.5, 0.75, 0.75),(0.68,1.0,1.0),(1.0,0.0,1.0)),
             'green': ((0.0, 0.0, 0.0), (0.25, 0.40, 0.3), (0.5, 0.75, 0.75),(0.68,1,1.0),(1.0,1.0,1.0)),
             'blue':  ((0.0, 0.0, 0.0), (0.25, 0.40, 0.3), (0.5, 0.75, 0.75),(0.68,1.0,1.0),(1.0,1.0,1.0))}
    
    p_cmap = mp.colors.LinearSegmentedColormap('p_colormap',cdict,256)
    
    
    ################### 
    # Data reading and preparation
    ###################
    
    w1,i1,q1,qc = p.genfromtxt(file_name,
    	usecols=(0,1,2,4), unpack=True, comments='#', skip_header=20)
    q1=abs(q1)
    
    s_idx=np.where(min(abs(w1-10*starto)) == abs(w1-10*starto)) #in Angstrom
    e_idx=np.where(min(abs(w1-10*endo)) == abs(w1-10*endo)) #in Angstrom
    
    w1=w1[s_idx[0][0]-1:]
    i1=i1[s_idx[0][0]-1:]
    q1=q1[s_idx[0][0]-1:]
    dwl=abs((w1[0]-w1[1])/10.) #wavelength step [nm]
    
    
    step=int((endo-starto)/dwl/lines)
    end=starto*10+lines*step*10.*dwl #in Angstrom
    e_idx=np.where(min(abs(w1-end)) == abs(w1-end)) #in Angstrom
    w1=w1[:e_idx[0][0]-1]
    i1=i1[:e_idx[0][0]-1]
    q1=q1[:e_idx[0][0]-1]
    qc=qc[:e_idx[0][0]-1]
    
    wl_chunk=(w1.reshape(lines,step))/10. # in nm
    intensity_chunk=i1.reshape(lines,step)
    q_chunk=abs(q1.reshape(lines,step))
    qc_chunk=abs(qc.reshape(lines,step))
    
    qc_chunk/=qc_chunk.max()
    q_chunk/=q_chunk.max()
    
    
    q_final=q_chunk+qc_chunk
    q_final/=q_final.max()
    
    
    spectrum=np.ones([step,4]) #number for color and intensity
    
    #### PLOTTING
    
    f=p.figure(num=0) 
    f.set_size_inches([xsize*zoom,ysize*zoom],forward=True) 
    
    gs = gridspec.GridSpec(lines,2)
    gs.update(left=0.01, right=0.99,top=0.99,bottom=0.01)  
    
    p.clf()
    p.cla() 
    
    # =============================================================================
    # intensity spectrum
    
    for i in range(lines):
     ax=p.subplot(gs[i,0]) 
     ax.get_xaxis().set_visible(False)
     ax.get_yaxis().set_visible(False)
    
     spectrum_plot=1.*spectrum
     for l in range(np.size(wl_chunk[i,:])):
        spectrum_plot[l,0:3]=spectral_color(wl_chunk[i,l])
        spectrum_plot[l,3]=intensity_chunk[i,l]
     spectrum_plot=np.array([spectrum_plot,spectrum_plot])
     im=ax.imshow(spectrum_plot.repeat(200,axis=0))
    
    # =============================================================================
    # second solar spectrum
    
    for i in range(lines):
     ax=p.subplot(gs[i,1]) 
     ax.get_xaxis().set_visible(False)
     ax.get_yaxis().set_visible(False)
    
    
     spectrum_plot=1.*spectrum
     
         
     for l in range(np.size(wl_chunk[i,:])):
    
        spectrum_plot[l,:]=p_cmap(q_final[i,l])
        
     spectrum_plot=np.array([spectrum_plot,spectrum_plot])
     im=ax.imshow(spectrum_plot.repeat(200,axis=0))
     
        
    f.savefig(save_file_name, dpi=300)
