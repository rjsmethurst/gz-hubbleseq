""" This python code will take the ../data/GZ_hubbleseq_sample_featured_faceon_spiralarms.fits
file and priduce a contour plot for the locations of the arm winding and bulge size and save
it as bulge_armwinding.pdf"""

import numpy as np 
from matplotlib import pyplot as plt 
from astropy.table import Table 
import os
from hist2d import hist2d
from scipy.stats import binned_statistic as bs

# Here redefine the matplotlib defaults to make nice plots
os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'

plt.rc('figure', facecolor='none', edgecolor='none', autolayout=True)
plt.rc('path', simplify=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes', labelsize='large', facecolor='none', linewidth=0.7, color_cycle = ['k', 'r', 'g', 'b', 'c', 'm', 'y'])
plt.rc('xtick', labelsize='medium')
plt.rc('ytick', labelsize='medium')
plt.rc('lines', markersize=4, linewidth=1, markeredgewidth=0.2)
plt.rc('legend', numpoints=1, frameon=False, handletextpad=0.3, scatterpoints=1, handlelength=2, handleheight=0.1)
plt.rc('savefig', facecolor='none', edgecolor='none', frameon='False')

params =   {'font.size' : 16,
            'xtick.major.size': 8,
            'ytick.major.size': 8,
            'xtick.minor.size': 3,
            'ytick.minor.size': 3,
            }
plt.rcParams.update(params) 

def stdem(x):
  ''' Function to calculate standard error on the median i.e. stdev(x)/sqrt(N)'''
  return np.std(x)/np.sqrt(float(len(x)))

# Load data from fits table

data = Table.read("../data/GZ_hubbleseq_sample_featured_faceon_spiralarms.fits")

# Calculate bulge size and arm winding as per equations (1) and (2) respectively 
bulgesize = 0.0*data['t05_bulge_prominence_a10_no_bulge_debiased'] + 0.2*data['t05_bulge_prominence_a11_just_noticeable_debiased'] + 0.8*data['t05_bulge_prominence_a12_obvious_debiased'] + 1.0*data['t05_bulge_prominence_a13_dominant_debiased']

armwind = 0.0*data['t10_arms_winding_a30_loose_debiased'] + 0.5*data['t10_arms_winding_a29_medium_debiased'] + 1.0*data['t10_arms_winding_a28_tight_debiased']

# Make the pretty python plot with all the data

medianb, binedgesb, binnumberb = bs(bulgesize, armwind, 'median', bins=10, range=(-0.05,1.05))
semb, binedgesb, binnumberb = bs(bulgesize, armwind, stdem, bins=10, range=(-0.05,1.05))

mediana, binedgesa, binnumbera = bs(armwind, bulgesize, 'median', bins=10, range=(-0.05,1.05))
sema, binedgesa, binnumbera = bs(armwind, bulgesize, stdem, bins=10, range=(-0.05,1.05))


binmiddleb = binedgesb[:-1] + np.diff(binedgesb)/2.
binmiddlea = binedgesa[:-1] + np.diff(binedgesa)/2.

x = np.linspace(-0.05, 1.05, 100)

plt.figure(figsize=(6.5,6))
ax = plt.subplot(111)
# ax.contourf(Xbins, Ybins, H.T, origin='lower', cmap=plt.cm.binary, alpha=0.2)
# ax.contour(Xbins, Ybins, H.T, origin='lower', colors='k')
csl = hist2d(bulgesize, armwind, bins=21, range=((-0.05,1.05),(-0.05,1.05)), smooth=0.5,
           ax=ax, plot_datapoints=True, plot_density=True,
           plot_contours=True, fill_contours=True, data_kwargs={'color':"k", "alpha":0.5})
ax.plot(x, x, color='k', linestyle=':')
ax.errorbar(binmiddleb, medianb, yerr=semb, ecolor='r', marker='x', c='r')
ax.errorbar(mediana, binmiddlea, xerr=sema, ecolor='b', marker='x', c='b')
ax.clabel(csl, inline=1, fontsize=12, fmt='%i')
ax.set_xlabel(r'$B_{avg}$')
ax.set_ylabel(r'$w_{avg}$')
# ax.set_xlim(0,1)
# ax.set_ylim(0,1)
ax.minorticks_on()
ax.tick_params(axis='both', which='both', direction='in', top='on', right='on')
plt.tight_layout()
plt.savefig('../bulge_armwinding_rolling_median.pdf')



# Make the pretty python plot with data split by p_bar < 0.2 and p_bar > 0.5

nobar = np.where(data['t03_bar_a06_bar_debiased']<0.2)
# There are 2166 galaxies satisfying this condition 

median_nobarb, binedges_nobarb, binnumber_nobarb = bs(bulgesize[nobar], armwind[nobar], 'median', bins=10, range=(-0.05,1.05))
sem_nobarb, binedges_nobarb, binnumber_nobarb = bs(bulgesize[nobar], armwind[nobar], stdem, bins=10, range=(-0.05,1.05))

median_nobara, binedges_nobara, binnumber_nobara = bs(armwind[nobar], bulgesize[nobar],  'median', bins=10, range=(-0.05,1.05))
sem_nobara, binedges_nobara, binnumber_nobara = bs(armwind[nobar], bulgesize[nobar], stdem, bins=10, range=(-0.05,1.05))

binmiddle_nobarb = binedges_nobarb[:-1] + np.diff(binedges_nobarb)/2.
binmiddle_nobara = binedges_nobara[:-1] + np.diff(binedges_nobara)/2.


bar = np.where(data['t03_bar_a06_bar_debiased']>0.5)
# There are 1364 galaxies satisfying this condition 
median_barb, binedges_barb, binnumber_barb = bs(bulgesize[bar], armwind[bar], 'median', bins=10, range=(-0.05,1.05))
sem_barb, binedges_barb, binnumber_barb = bs(bulgesize[bar], armwind[bar], stdem, bins=10, range=(-0.05,1.05))

median_bara, binedges_bara, binnumber_bara = bs(armwind[bar], bulgesize[bar],  'median', bins=10, range=(-0.05,1.05))
sem_bara, binedges_bara, binnumber_bara = bs(armwind[bar], bulgesize[bar], stdem, bins=10, range=(-0.05,1.05))

binmiddle_barb = binedges_barb[:-1] + np.diff(binedges_barb)/2.
binmiddle_bara = binedges_bara[:-1] + np.diff(binedges_bara)/2.



plt.figure(figsize=(14,6))
ax1 = plt.subplot(121)
csl = hist2d(bulgesize[nobar], armwind[nobar], bins=21, range=((-0.05,1.05),(-0.05,1.05)), smooth=0.7,
           ax=ax1, plot_datapoints=True, plot_density=True,
           plot_contours=True, fill_contours=True, data_kwargs={'color':"k", "alpha":0.5})
ax1.plot(x, x, color='k', linestyle=':')
ax1.errorbar(binmiddle_nobarb, median_nobarb, yerr=sem_nobarb, ecolor='r', marker='x', c='r')
ax1.errorbar(median_nobara, binmiddle_nobara, xerr=sem_nobara, ecolor='b', marker='x', c='b')
ax1.clabel(csl, inline=1, fontsize=12, fmt='%i')
ax1.set_xlabel(r'$B_{avg}$')
ax1.set_ylabel(r'$w_{avg}$')
# ax.set_xlim(0,1)
# ax.set_ylim(0,1)
ax1.minorticks_on()
ax1.text(0.05, 0.9, r'$\rm{p}_{\rm{bar}} < 0.2$', transform=ax1.transAxes, fontsize=14)
ax1.tick_params(axis='both', which='both', direction='in', top='on', right='on')

ax2 = plt.subplot(122)
# ax.contourf(Xbins, Ybins, H.T, origin='lower', cmap=plt.cm.binary, alpha=0.2)
# ax.contour(Xbins, Ybins, H.T, origin='lower', colors='k')
csl = hist2d(bulgesize[bar], armwind[bar], bins=21, range=((-0.05,1.05),(-0.05,1.05)), smooth=0.7,
           ax=ax2, plot_datapoints=True, plot_density=True,
           plot_contours=True, fill_contours=True, data_kwargs={'color':"k", "alpha":0.5})
ax2.plot(x, x, color='k', linestyle=':')
ax2.errorbar(binmiddle_barb, median_barb, yerr=sem_barb, ecolor='r', marker='x', c='r')
ax2.errorbar(median_bara, binmiddle_bara, xerr=sem_bara, ecolor='b', marker='x', c='b')
ax2.clabel(csl, inline=1, fontsize=12, fmt='%i')

ax2.set_xlabel(r'$B_{avg}$')
ax2.set_ylabel(r'$w_{avg}$')
# ax.set_xlim(0,1)
# ax.set_ylim(0,1)
ax2.minorticks_on()
ax2.text(0.05, 0.9, r'$\rm{p}_{\rm{bar}} > 0.5$', transform=ax2.transAxes, fontsize=14)
ax2.tick_params(axis='both', which='both', direction='in', top='on', right='on')

plt.tight_layout()
plt.savefig('../bulge_armwinding_split_bar_rolling_median.pdf')

