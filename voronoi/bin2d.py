#!/usr/bin/env python
# -----------------------------------------------------------------------------
# VORONOI.BIN2D
# Laura L Watkins [lauralwatkins@gmail.com]
# - converted from IDL code by Michele Cappellari (voronoi_2d_binning)
# -----------------------------------------------------------------------------

from numpy import *
from matplotlib.pyplot import *
from .weighted_centroid import *
from .bin_roundness import *
from .accretion import *
from .reassign_bad_bins import *
from .cvt_equal_mass import *
from .bin_quantities import *


def bin2d(x, y, signal, noise, targetsn, cvt=True, wvt=False, quiet=True,
    graphs=True):
    
    """
    This is the main program that has to be called from external programs.
    It simply calls in sequence the different steps of the algorithms
    and optionally plots the results at the end of the calculation.
    
    INPUTS
      x        : x-coordinates of pixels
      y        : y-coordinates of pixels
      signal   : signal in pixels
      noise    : noise in pixels
      targetsn : target S/N required
    
    OPTIONS
      cvt      : use Modified-Lloyd algorithm [default True]
      wvt      : use additional modification by Diehl & Statler [default False]
      quiet    : supress output [default True]
      graphs   : show results graphically [default True]
    """
    
    
    npix = x.size
    if y.size != x.size or signal.size != x.size or noise.size != x.size:
        print("ERROR: input vectors (x, y, signal, noise) must have same size")
        return
    if any(noise < 0):
        print("ERROR: noise cannot be negative")
        return
    
    # prevent division by zero for pixels with signal=0 and
    # noise=sqrt(signal)=0 as can happen with X-ray data
    noise[noise==0] = noise[noise>0].min() * 1e-9
    
    # Perform basic tests to catch common input errors
    if signal.sum()/sqrt((noise**2).sum()) < targetsn:
        print("Not enough S/N in the whole set of pixels. " \
            + "Many pixels may have noise but virtually no signal. " \
            + "They should not be included in the set to bin, " \
            + "or the pixels should be optimally weighted." \
            + "See Cappellari & Copin (2003, Sec.2.1) and README file.")
        return
    if (signal/noise).min() > targetsn:
        print("EXCEPTION: all pixels have enough S/N -- binning not needed")
        return
    
    if not quiet: print("Bin-accretion...")
    clas = accretion(x, y, signal, noise, targetsn, quiet=quiet)
    if not quiet: print("{:} initial bins\n".format(clas.max()))
    
    if not quiet: print("Reassign bad bins...")
    xnode, ynode = reassign_bad_bins(x, y, signal, noise, targetsn, clas)
    if not quiet: print("{:} good bins\n".format(xnode.size))
    
    if cvt:
        if not quiet: print("Modified Lloyd algorithm...")
        scale, iters = cvt_equal_mass(x, y, signal, noise, xnode, ynode,
            quiet=quiet, wvt=wvt)
        if not quiet: print("  iterations: {:}".format(iters-1))
    else:
        scale = 1.
    
    if not quiet: print("Recompute bin properties...")
    clas, xbar, ybar, sn, area = bin_quantities(x, y, signal, noise, xnode,
        ynode, scale)
    unb = where(area == 1)[0]
    binned = where(area != 1)[0]
    if not quiet: print("Unbinned pixels: {:} / {:}".format(unb.size, npix))
    fracscat = ((sn[binned]-targetsn)/targetsn*100.).std()
    if not quiet: print("Fractional S/N scatter (%):", fracscat)
    
    if graphs:
        
        # set up plotting
        rc("font", family="serif")
        rc("text", usetex=True)
        rc("xtick", labelsize="8")
        rc("ytick", labelsize="8")
        rc("axes", labelsize="10")
        rc("legend", fontsize="9")
        
        # pixel map
        fig = figure(figsize=(4,3))
        fig.subplots_adjust(left=0.13, bottom=0.13, top=0.97, right=0.98)
        rnd = random.rand(xnode.size).argsort()      # randomize bin colors
        scatter(x, y, lw=0, c=rnd[clas])
        plot(xnode, ynode, "k+", ms=2)
        xlim(x.min()-x.ptp()*0.05, x.max()+x.ptp()*0.05)
        ylim(y.min()-y.ptp()*0.05, y.max()+y.ptp()*0.05)
        xlabel("coordinate 1")
        ylabel("coordinate 2")
        show()
        
        # signal-to-noise profile
        fig = figure(figsize=(4,3))
        fig.subplots_adjust(left=0.12, bottom=0.13, top=0.97, right=0.97)
        rad = sqrt(xbar**2 + ybar**2)     # use centroids, NOT generators
        rmin = max(0., rad.min()-rad.ptp()*0.05)
        rmax = rad.max()+rad.ptp()*0.05
        plot([rmin, rmax], ones(2)*targetsn, c="k", lw=2, alpha=0.8)
        scatter(rad[binned], sn[binned], lw=0, c="b", alpha=0.8)
        if unb.size > 0: scatter(rad[unb], sn[unb], lw=0, c="r", alpha=0.8)
        xlim(rmin, rmax)
        ylim(0., sn.max()*1.05)
        xlabel(r"$R_{\rm bin}$")
        ylabel(r"$SN_{\rm bin}$")
        show()
    
    
    return clas, xnode, ynode, sn, area, scale
 