#!/usr/bin/env python
# -----------------------------------------------------------------------------
# VORONOI.ACCRETION
# Laura L Watkins [lauralwatkins@gmail.com]
# - converted from IDL code by Michele Cappellari (bin2d_accretion)
# -----------------------------------------------------------------------------

from numpy import *
from .bin_roundness import bin_roundness


def accretion(x, y, signal, noise, targetsn, pixelsize=False, quiet=False):
    
    """
    Initial binning -- steps i-v of eq 5.1 of Cappellari & Copin (2003)
    
    INPUTS:
      x        : x coordinates of pixels to bin
      y        : y coordinates of pixels to bin
      signal   : signal associated with each pixel
      noise    : noise (1-sigma error) associated with each pixel
      targetsn : desired signal-to-noise ration in final 2d-binned data
    
    OPTIONS:
      pixelsize : pixel scale of the input data
      quiet     : if set, suppress printed outputs
    """
    
    
    n = x.size
    clas = zeros(x.size, dtype="<i8")   # bin number of each pixel
    good = zeros(x.size, dtype="<i8")   # =1 if bin accepted as good
    
    # for each point, find distance to all other points and select minimum
    # (robust but slow way of determining the pixel size of unbinned data)
    if not pixelsize:
        dx = 1.e30
        for j in range(x.size-1):
            d = (x[j] - x[j+1:])**2 + (y[j] - y[j+1:])**2
            dx = min(d.min(), dx)
        pixelsize = sqrt(dx)
    
    # start from the pixel with highest S/N
    sn = (signal/noise).max()
    currentbin = (signal/noise).argmax()
    
    # rough estimate of the expected final bin number
    # This value is only used to have a feeling of the expected
    # remaining computation time when binning very big dataset.
    wh = where(signal/noise<targetsn)
    npass = size(where(signal/noise >= targetsn))
    maxnum = int(round( (signal[wh]**2/noise[wh]**2).sum()/targetsn**2 ))+npass
    
    # first bin assigned CLAS = 1 -- with N pixels, get at most N bins
    for ind in range(1, n+1):
        
        if not quiet:
            print("  bin: {:} / {:}".format(ind, maxnum))
        
        # to start the current bin is only one pixel
        clas[currentbin] = ind
        
        # centroid of bin
        xbar = x[currentbin]
        ybar = y[currentbin]
        
        while True:
            
            # stop if all pixels are binned
            unbinned = where(clas == 0)[0]
            if unbinned.size == 0: break
            
            # find unbinned pixel closest to centroid of current bin
            dist = (x[unbinned]-xbar)**2 + (y[unbinned]-ybar)**2
            mindist = dist.min()
            k = dist.argmin()
            
            # find the distance from the closest pixel to the current bin
            mindist = ((x[currentbin]-x[unbinned[k]])**2 \
                + (y[currentbin]-y[unbinned[k]])**2).min()
            
            # estimate roundness of bin with candidate pixel added
            nextbin = append(currentbin, unbinned[k])
            roundness = bin_roundness(x[nextbin], y[nextbin], pixelsize)
            
            # compute sn of bin with candidate pixel added
            snold = sn
            sn = signal[nextbin].sum()/sqrt((noise[nextbin]**2).sum())
            
            # Test whether the CANDIDATE pixel is connected to the
            # current bin, whether the POSSIBLE new bin is round enough
            # and whether the resulting S/N would get closer to targetsn
            if sqrt(mindist) > 1.2*pixelsize or roundness > 0.3 \
                or abs(sn-targetsn) > abs(snold-targetsn):
                if (snold > 0.8*targetsn):
                    good[currentbin] = 1
                break
            
            # if all the above tests are negative then accept the CANDIDATE
            # pixel, add it to the current bin, and continue accreting pixels
            clas[unbinned[k]] = ind
            currentbin = nextbin
            
            # update the centroid of the current bin
            xbar = x[currentbin].mean()
            ybar = y[currentbin].mean()
        
        # get the centroid of all the binned pixels
        binned = where(clas != 0)[0]
        unbinned = where(clas == 0)[0]
        
        # stop if all pixels are binned
        if unbinned.size == 0: break
        xbar = x[binned].mean()
        ybar = y[binned].mean()
        
        # find the closest unbinned pixel to the centroid of all
        # the binned pixels, and start a new bin from that pixel
        k = ((x[unbinned]-xbar)**2 + (y[unbinned]-ybar)**2).argmin()
        currentbin = unbinned[k]    # the bin is initially made of one pixel
        sn = signal[currentbin] / noise[currentbin]
    
    # set to zero all bins that did not reach the target S/N
    clas = clas*good
    
    return clas
