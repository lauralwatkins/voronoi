#!/usr/bin/env python
# -----------------------------------------------------------------------------
# VORONOI.REASSIGN_BAD_BINS
# Laura L Watkins [lauralwatkins@gmail.com]
# - converted from IDL code by Michele Cappellari (bin2d_reassign_bad_bins)
# -----------------------------------------------------------------------------

from numpy import *


def reassign_bad_bins(x, y, signal, noise, targetsn, clas):
    
    """
    Reassign bad bins -- steps vi-vii of eq 5.1 of Cappellari & Copin (2003)
    
    INPUTS
      x        : x-coordinates of pixels
      y        : y-coordinates of pixels
      signal   : signal in pixels
      noise    : noise in pixels
      targetsn : target signal/noise required
      clas     : bin number for each pixel
    """
    
    
    # get number of pixels in each bin (clas=0 are unassigned pixels)
    area, lim = histogram(clas, bins=clas.ptp(), range=(0.5,clas.max()+0.5))
    cent = (lim[:-1]+lim[1:])/2.
    
    # indices of good bins
    good = where(area > 0)[0]
    
    # centroids of good bins
    xnode = zeros(good.size)
    ynode = zeros(good.size)
    for j in range(good.size):
        p = where(clas == cent[good[j]])
        xnode[j] = x[p].mean()
        ynode[j] = y[p].mean()
    
    # reassign pixels of bins with S/N < targetSN to closest good bin
    bad = where(clas == 0)[0]
    for j in range(bad.size):
        index = ((x[bad[j]] - xnode)**2 + (y[bad[j]] - ynode)**2).argmin()
        clas[bad[j]] = good[index] + 1
    
    # recompute all centroids of the reassigned bins
    # these will be used as starting points for the CVT
    area, lim = histogram(clas, bins=clas.ptp()+1, range=(0.5,clas.max()+0.5))
    cent = (lim[:-1]+lim[1:])/2.
    
    good = where(area > 0)[0]
    for j in range(good.size):
        p = where(clas == cent[good[j]])
        xnode[j] = x[p].mean()
        ynode[j] = y[p].mean()
    
    return xnode, ynode
