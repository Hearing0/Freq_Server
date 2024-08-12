#!/usr/bin/env python3

import numpy as np
from numpy import sin,cos

import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib


params = {'text.usetex': True,
          'text.latex.preamble': [r'\usepackage{cmbright} \usepackage{amsmath} \usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'],
          'axes.linewidth':1,
          'axes.facecolor':(1.0, 1.0, 1.0),
          'savefig.facecolor':(1.0, 1.0, 1.0)
          }

plt.rcParams.update(params)

col1='cadetblue'
col2='tomato'
PSUblue=(38/255,62/255,126/255)

fig=plt.figure(figsize=[6,4])
ax=fig.add_subplot(1,1,1)

parser = argparse.ArgumentParser(description='bin_radar_to_map argument parser')
parser.add_argument('num_avg', type=int)

args = parser.parse_args()
num_avg=args.num_avg

freq=100
npts=1000


pts=np.linspace(0,freq*np.pi,num=npts)    
spectrum=np.complex(npts)

for j in range(num_avg):

    sig=np.cos(pts)+np.random.normal(0.0,5.0,size=npts)
    spectrum+=np.fft.fft(sig)


ax.plot(pts,np.abs(spectrum)/num_avg)


plt.show()