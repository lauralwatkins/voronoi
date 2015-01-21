#!/usr/bin/env python
# -----------------------------------------------------------------------------
# VORONOI.WEIGHTED_CENTROID
# Laura L Watkins [lauralwatkins@gmail.com]
# - converted from IDL code by Michele Cappellari (bin2d_weighted_centroid)
# -----------------------------------------------------------------------------


def weighted_centroid(x, y, density):
    
    """
    Computes weighted centroid of a bin -- eq 4 of Cappellari & Copin (2003).
    
    INPUTS
      x       : x-coordinate of pixels in bin
      y       : y-coordinate of pixels in bin
      density : pixel weights
    """
    
    
    mass = density.sum()
    xbar = (x*density).sum()/mass
    ybar = (y*density).sum()/mass
    
    return xbar, ybar
