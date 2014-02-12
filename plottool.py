#!/usr/bin/env python

from pylab import *

lw = 2.0

# for two-axes technique, see
#   http://matplotlib.sourceforge.net/plot_directive/mpl_examples/api/two_scales.py
#   http://matplotlib.sourceforge.net/plot_directive/mpl_examples/api/fahrenheit_celcius_scales.py

def plottwo(fig, x, y1, y2, label1, label2, y1min=-inf, y1max=inf, y2min=-inf, y2max=inf):
  import matplotlib.pyplot as plt

  ax1 = fig.add_subplot(111) # left axis
  ax2 = ax1.twinx()          # right axis

  ax1.plot(x,y1,'k',linewidth=lw)
  ax2.plot(x,y2,'k--',linewidth=lw)

  ax1.set_xlabel("x  (km)")

  if y1min == -inf:
    y1min = y1.min()
  if y1max == inf:
    y1max = y1.max()
  ax1.set_ylim(y1min, y1max)

  if y2min == -inf:
    y2min = y2.min()
  if y2max == inf:
    y2max = y2.max()
  ax2.set_ylim(y2min, y2max)

  # make the y-axis label and tick labels match the plot colors
  ax1.set_ylabel(label1, color='k')
  for tl in ax1.get_yticklabels():
      tl.set_color('k')
  ax2.set_ylabel(label2, color='k')
  for tl in ax2.get_yticklabels():
      tl.set_color('k')

  return ax1, ax2

