#!/usr/bin/env python
# -----------------------------------------------------------------------------
# VORONOI.BIN_ROUNDESS
# Laura L Watkins [lauralwatkins@gmail.com]
# - converted from IDL code by Michele Cappellari (bin2d_roundness)
# -----------------------------------------------------------------------------

from numpy import *


def bin_roundness(x, y, pixelsize):
    
    """
    Computes roundness of a bin -- eq 5 of Cappellari & Copin (2003)
    
    INPUTS
      x         : x-coordinates of pixels in bin
      y         : y-coordinates of pixels in bin
      pixelsize : size of pixels
    """
    
    
    equivalentradius = sqrt(x.size/pi)*pixelsize
    
    # geometric centroid
    xbar = x.mean()
    ybar = y.mean()
    
    maxdistance = sqrt((x-xbar)**2 + (y-ybar)**2).max()
    roundness = maxdistance/equivalentradius - 1.
    
    return roundness
