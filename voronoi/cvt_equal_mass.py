#!/usr/bin/env python
# -----------------------------------------------------------------------------
# VORONOI.CVT_EQUAL_MASS
# - converted from IDL code by Michele Cappellari (bin2d_cvt_equal_mass)
# -----------------------------------------------------------------------------

import numpy as np
from .weighted_centroid import weighted_centroid


def cvt_equal_mass(x, y, signal, noise, xnode, ynode, quiet=True, wvt=False):
    
    """
    Modified Lloyd algorithm -- section 4.1 of Cappellari & Copin (2003).
    When the keyword wvt is set, the routine includes the modification
    proposed by Diehl & Statler (2006).
    
    INPUTS
      x      : x-coordinates of pixels
      y      : y-coordinates of pixels
      signal : signal in pixels
      noise  : noise in pixels
      xnode  : x-coordinates of bins
      ynode  : y-coordinates of bins
    
    OPTIONS
      quiet  : suppress output [default True]
      wvt    : use modification of Diehl & Statler (2006) [default False]
    """
    
    
    clas = np.zeros(len(signal))   # see beginning of section 4.1 of CC03
    if wvt: dens = np.ones(len(signal))
    else: dens = signal**2/noise**2
    scale = 1                   # start with the same scale length for all bins
    sn = np.zeros(len(xnode))
    
    iters = 1
    diff = 1
    while diff!=0:
        
        xold = xnode.copy()
        yold = ynode.copy()
        
        # computes (weighted) voronoi tessellation of the pixels grid
        for j in range(len(signal)):
            index = (((x[j]-xnode)/scale)**2+((y[j]-ynode)/scale)**2).argmin()
            clas[j] = index
        
        # Computes centroids of the bins, weighted by dens^2.
        # Exponent 2 on the density produces equal-mass Voronoi bins.
        # The geometric centroids are computed if /wvt keyword is set.
        
        area, lim = np.histogram(clas, bins=int(clas.ptp())+1,
            range=(clas.min()-0.5, clas.max()+0.5))
        cent = (lim[:-1]+lim[1:])/2
        
        nonzero = np.where(area>0)[0]       # check for zero-size voronoi bins
        for j in range(len(nonzero)):
            k = nonzero[j]                  # only loop over nonzero bins
            pix = clas==cent[k]
            xnode[k], ynode[k] = weighted_centroid(x[pix],y[pix],dens[pix]**2)
            sn[k] = signal[pix].sum()/np.sqrt((noise[pix]**2).sum())
        
        if wvt: scale = np.sqrt(area/sn)    # eq 4 of Diehl & Statler (2006)
        diff = ((xnode-xold)**2 + (ynode-yold)**2).sum()
        iters = iters + 1
        
        if not quiet:
            print("  iteration: {:}  difference: {:}".format(iters, diff))
    
    # only return the generators of the nonzero voronoi bins
    xnode = xnode[nonzero]
    ynode = ynode[nonzero]
    
    return scale, iters
