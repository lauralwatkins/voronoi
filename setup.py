import os
from setuptools import setup, Extension, find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
	name='voronoi',
	version='0.0',
	description='This code performs a Voronoi 2D-binning algorithm on an input set of pixels. This is a Python version of the IDL implementation by Michele Cappellari. A Python version of the code is also available at that URL, though it was not at the time that this code was written.',
	author='Laura L Watkins',
	author_email='lauralwatkins@gmail.com',
    package_dir = {
        'voronoi': 'voronoi',
        },
    packages=["voronoi"],
)