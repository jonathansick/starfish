#!/bin/bash
#
#     This file is part of StarFISH
#     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
#     
# -=> StarFISH is free software; you can redistribute it 
# -=> and/or modify it under the terms of the GNU General Public 
# -=> License (GPL) as published by the Free Software Foundation;
# -=> either version 2 of the License, or (at your option) any 
# -=> later version.
#
# Usage: ./mkgif.sh <synthCMD_prefix>
#
# mkgif creates images of the synthetic CMD files output by 
# synth.  Just copy this script and pxl2ppm.awk into your 
# synthCMD directory.  You will probably have to modify it to 
# suit your needs:
#
# The original script assumes 3 CMDs per isochrone, and that 
# the synthCMD files have suffixes ".ub", ".bv", and ".vi".
#
# adjust the contrast by modifying the 'max' parameter, which 
# controls the pixel value that will be fully saturated in the 
# image.
#
# adjust the xd and yd parameters to reflect the geometry of 
# your synthCMDs.
#
# Note: the script assumes you have imageMagick installed; if 
# you don't have it, you can just comment out the "convert" 
# lines and the "rm -f *.ppm" line.  This will leave you with 
# the ppm images (note that these will be upside-down).
#
test $# -ne 1 && echo "Usage: ./mkgif.sh <synthcmd_prefix>" && exit

gawk -f pxl2ppm.awk -v xd=70 -v yd=200 -v clr=0 -v max=0.00008 -v out=$1.ub.ppm $1.ub
gawk -f pxl2ppm.awk -v xd=50 -v yd=200 -v clr=0 -v max=0.00008 -v out=$1.bv.ppm $1.bv
gawk -f pxl2ppm.awk -v xd=60 -v yd=200 -v clr=0 -v max=0.00008 -v out=$1.vi.ppm $1.vi

convert $1.ub.ppm -flip $1.ub.png
convert $1.bv.ppm -flip $1.bv.png
convert $1.vi.ppm -flip $1.vi.png

rm -f *.ppm
