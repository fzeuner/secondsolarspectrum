# -*- coding: utf-8 -*-

"""
Created on 22.10.2020

@author: fzeuner

This code will plot the solar intensity spectrum in color 
and the second solar spectrum at grey scale.

"""

from scndspectrum import scndspectrum


################### 
# Change parameters
###################

# second solar spectrum
# Get data from http://www.irsol.usi.ch/data-archive/second-solar-spectrum-ss2-atlas/
file_name='/home/franziskaz/Downloads/SSSatlas.txt'
save_file_name='/home/franziskaz/2ndspectrum_poster.png' # figure save name

starto=415. # wavelength start [nm]
endo=480 # wavelength end [nm] - will be handled flexible to serve number of lines
lines=19 # number of lines in the plot

# use below parameters for full spectrum - but you will need some memory!
# starto=400. # wavelength start [nm]
# endo=650 # wavelength end [nm] - will be handled flexible to serve number of lines
# lines=42 # number of lines in the plot

# size of the figure [inches]
xsize=7
ysize=4.5

zoom=1 # zoom factor 

scndspectrum(file_name,save_file_name,starto,endo,lines,xsize,ysize,zoom)
