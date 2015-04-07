#!/usr/bin/env python
#-*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Cairo')
import matplotlib.pyplot as plt
import numpy as np
from SkyNet import *
import sys
from scipy.interpolate import interp1d
import os

deltaLogT = 0.02

def make_plot(in_path, out_png_path, out_txt_path):
  out = NetworkOutput.ReadFromFile(in_path)
  original_times = np.array(out.Times())

  #initial_time = original_times[0]
  initial_time = np.max([1.0e-3, original_times[0]])
  final_time = original_times[-1]

  num = (np.log(final_time) - np.log(initial_time)) / deltaLogT
  num = int(num + 1.5)

  actual_deltaLogT = (np.log(final_time) - np.log(initial_time)) / (num - 1)
  times = np.exp(np.log(initial_time) + np.arange(num) * actual_deltaLogT)
  times[0] = initial_time

  original_temp = np.array(out.TemperatureVsTime()) * 1.0E9
  original_rho = np.array(out.DensityVsTime())
  original_edot = np.array(out.HeatingRateVsTime())
  original_edot[original_edot <= 0.0] = None

  interp = interp1d(original_times, original_temp, bounds_error=False, fill_value=None)
  temp = interp(times)
  interp = interp1d(original_times, original_rho, bounds_error=False, fill_value=None)
  rho = interp(times)
  interp = interp1d(original_times, original_edot, bounds_error=False, fill_value=None)
  edot = interp(times)

  font = {'family' : 'sans-serif',
          'weight' : 'normal',
          'size'   : font_size}
  #text = {'hinting_factor': 8.0}

  matplotlib.rc('font', **font)
  #matplotlib.rc('text', **text)

  t0 = np.min(times)
  tidx0 = np.min(np.arange(len(times)))

  fig, ax1 = plt.subplots(figsize=(4.5, 2.75), dpi=100)
  ax2 = ax1.twinx()

  ax1.set_xlim(t0, 1E9)
  ax1.set_xscale('log')
  ax2.set_xlim(t0, 1E9)
  ax2.set_xscale('log')

  ax1.set_ylim(-28, 12)
  ax1.set_yscale('linear')
  ax2.set_ylim(3, 22)
  ax2.set_yscale('linear')

  l1 = ax1.plot(times, 8*np.log10(temp)-70, 'r-', label='Temperature [K]')
  l2 = ax1.plot(times, np.log10(rho), 'b-', label='Density [g / cm$^3\!$]')
  l3 = ax2.plot(times, np.log10(edot), 'k-', label='Heating rate [erg / s / g]')
  lines = l1 + l2 + l3
  labels = [l.get_label() for l in lines]

  #ax1.set_xlabel('Time [s]')
  ax1.text(0.5, -0.16, 'Time [s]',
    horizontalalignment="center",transform=ax1.transAxes)

  #ax1.set_ylabel('Temperature [K] / Density [g cm$^{-3}$]')
  ax1.text(-0.16, 0.43, u"8 log10 Temp. – 70  /  log10 Density",
    verticalalignment="center",rotation='vertical',transform=ax1.transAxes)

  #ax2.set_ylabel('dE/dt [?]')
  ax2.text(1.10, 0.5, 'log10 Heating rate',
    verticalalignment="center",rotation='vertical',transform=ax2.transAxes)

  ax1.legend(lines, labels, loc='upper right', fontsize=font_size, labelspacing=0.3,frameon=False)
  plt.minorticks_off()
  fig.subplots_adjust(top = 0.97, left = 0.125, right = 0.89, bottom = 0.15)

  ax1.set_yticks([-28, -24, -20, -16, -12, -8, -4, 0, 4, 8, 12])
  ax1.set_yticks([], minor=True)

  ax2.set_yticks([22, 20, 18, 16, 14, 12, 10, 8, 6, 4])
  ax2.set_yticks([], minor=True)

  ax1.set_xticks([1E-3, 1E-1, 1E1, 1E3, 1E5, 1E7, 1E9])
  ax1.set_xticks([1E-2, 1E0, 1E2, 1E4, 1E6, 1E8], minor=True)

  ax2.set_xticks([1E-3, 1E-1, 1E1, 1E3, 1E5, 1E7, 1E9])
  ax2.set_xticks([1E-2, 1E0, 1E2, 1E4, 1E6, 1E8], minor=True)

  # Get the x and y data and transform it into pixel coordinates
  x1, y1 = l1[0].get_data()
  #y1[np.isnan(y1)] = -1.0
  x2, y2 = l2[0].get_data()
  x3, y3 = l3[0].get_data()
  pixel_coords1 = ax1.transData.transform(np.vstack([x1,y1]).T)
  pixel_coords2 = ax1.transData.transform(np.vstack([x2,y2]).T)
  pixel_coords3 = ax2.transData.transform(np.vstack([x3,y3]).T)

  #pixel_coords1[:,0] = pixel_coords2[:,1] # x in 1 is nan when temperature is small
  coords = np.append(pixel_coords1, pixel_coords2[:,1:2], axis=1)
  coords = np.append(coords, pixel_coords3[:,1:2], axis=1)

  plt.savefig(out_png_path)
  np.savetxt(out_txt_path, coords,
  header='x                      temperature              density                  dE/dt')

font_size = 10;
h5path = sys.argv[1]
dirname = os.path.dirname(os.path.abspath(h5path))

make_plot(h5path, dirname + '/plot_background.png', dirname + '/plot_coords.txt')
