#!/usr/bin/env python
# -----------------------------------------------------------------------------
# VORONOI.EXAMPLE
# Laura L Watkins [lauralwatkins@gmail.com]
# - converted from IDL code by Michele Cappellari (voronoi_2d_binning_example)
# -----------------------------------------------------------------------------

from astropy import table
from bin2d import bin2d


def example():
    
    """
    Example program to show how binning routine can be implemented.
    The test data is contained in example_input.dat:
      columns 1-4 contain respectively the x, y coordinates of each SAURON
      lens and the corresponding Signal and Noise.
    The resulting output is contained in example_output.dat:
      columns 1-3 contain the number of the bin to which the pixel has been
      assigned and the positions of the pixel.
    """
    
    
    # read in test data
    data = table.Table.read("example_input.dat",
        format="ascii.commented_header")
    
    # target signal-to-noise required
    targetsn = 50.
    
    # perform the actual computation
    pix_bin, bin_x, bin_y, bin_sn, bin_npix, scale = bin2d(data["x"],
        data["y"], data["signal"], data["noise"], targetsn)
    
    # write results to file
    result = table.Table()
    result["x"] = data["x"]
    result["y"] = data["y"]
    result["bin"] = pix_bin
    result.write("example_output.dat", format="ascii.fixed_width",
        bookend=False, delimiter=None)
    
    return
