VORONOI
=======

> **AUTHORS**  
Laura L Watkins (STScI), <lauralwatkins@gmail.com>

> **PUBLICATIONS**  
[Watkins et al. 2015, ApJ, submitted][watkins2015]

-------------------------------------------------------------------------------

CONTENTS
--------

* license and referencing
* code description
* requirements

-------------------------------------------------------------------------------

LICENSE AND REFERENCING
-----------------------

This code is released under a BSD 2-clause license.

If you use this code for your research, please cite:  
[Watkins et al. 2015, ApJ submitted][watkins2015]

Please don't forget to also cite [Cappellari & Copin (2003)][cappellari2003] who wrote the original IDL code upon which this Python version was based.

-------------------------------------------------------------------------------

CODE DESCRIPTION
----------------

This code performs a Voronoi 2D-binning algorithm on an input set of pixels. This is a Python version of [the IDL implementation](http://www-astro.physics.ox.ac.uk/~mxc/software/#binning) by Michele Cappellari. A Python version of the code is also available at that URL, though it was not at the time that this code was written.

example.py offers an example for how the code should be run; example\_in.dat and example\_out.dat are the input and output files for the example.

The code takes five arguments:
> x      : the x-coordinates of the input pixels  
  y      : the y-coordinates of the input pixels  
  signal : the signal in each of the input pixels  
  noise  : the noise in each of the input pixels  
  target : the target signal/noise required in each bin

And returns six:
> pix\_bin : bin number into which each pixel was placed  
  bin\_x   : the x-coordinate of the output bins  
  bin\_y   : the y-coordinate of the output bins  
  bin\_sn  : the signal-to-noise level in each of the output bins  
  bin\_pix : the number of original pixels in each of the output bins  
  scale    : the scale


-------------------------------------------------------------------------------


REQUIREMENTS
----------------------------------------

This code uses the standard python libraries numpy and matplotlib. The example also makes use of the astropy libraries for file i/o.


[cappellari2003]: http://adsabs.harvard.edu/abs/2003MNRAS.342..345C
[watkins2015]: http://adsabs.harvard.edu/abs/
