#!/usr/bin/env python

# import os
# import string
# 
# for f in os.listdir( os.path.dirname( __file__ ) ):
#     name, ext = string.split( f, "." )
#     if ( name != "__init__" and ext == "py" ):
#         tmp = __import__( name, locals(), globals(), [name], -1 )
#         vars()[name] = vars(tmp)[name]
# 
# del f, name, ext, os, string

from . import bin2d
