#!/bin/bash
#
# mkimages.sh:  Generate CMD images from a chifile.
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
#     Usage: ./mkgif.sh <output file prefix> <chifile>
#
#     This script must be customized to suit your needs.  You have 
#     to specify the CMD suffixes and the size of each CMD image.
#     You may also wish to adjust the pixel value ranges.
#
#     The script assumes you have ImageMagick installed; if you do 
#     not, comment out the indicated code block.
#
############################################################
## Usage message
##
test $# -ne 2 && echo "Usage: ./mkimages.sh <output_prefix> <chifile>" && exit
############################################################
## CUSTOMIZATION:  Adjust these settings to suit your needs:
####################################
## CMD indexes  
## (add or remove indexes to match your number of CMDs)
icmd=( 0 1 2 )
####################################
## CMD suffixes  
## (add or remove values to match your number of CMDs)
sfx=( "ub"  "bv"  "vi" )
####################################
## CMD suffixes (all in one string, delimited by '.')
cmdstring="ub.bv.vi"
####################################
## CMD pixel dimensions  
## (add or remove values to match your number of CMDs)
x=(  70   50   60 )
y=( 200  200  200 )
####################################
## image type names
typ1="mod"  typ2="dat"  typ3="chi"  typ4="dif"
####################################
## Pixel value ranges for model, data, chi, and diff images
max1=300    max2=300    max3=50     max4=200
##
## End of Customization section
############################################################

for i in ${icmd[*]}
do
  ((j=$i+1))
  stem1[$i]=$1.$typ1.${sfx[$i]}
  stem2[$i]=$1.$typ2.${sfx[$i]}
  stem3[$i]=$1.$typ3.${sfx[$i]}
  stem4[$i]=$1.$typ4.${sfx[$i]}

  gindex[$i]=grid$j.index
done

for i in ${icmd[*]}
do
  ((j=$i+1))
  # generate the *.pxl files from the chifile using the gridN.index files
  gawk -f grid.awk   -v out=${stem1[$i]}.pxl -v grid=${gindex[$i]} -v cmd=$j -v icol=3 -v nx=${x[$i]} -v ny=${y[$i]} $2
  gawk -f grid.awk   -v out=${stem2[$i]}.pxl -v grid=${gindex[$i]} -v cmd=$j -v icol=4 -v nx=${x[$i]} -v ny=${y[$i]} $2
  gawk -f grid.awk   -v out=${stem3[$i]}.pxl -v grid=${gindex[$i]} -v cmd=$j -v icol=5 -v nx=${x[$i]} -v ny=${y[$i]} $2
  gawk -f mkdiff.awk -v out=${stem4[$i]}.pxl -v grid=${gindex[$i]} -v cmd=$j -v p0=100 -v nx=${x[$i]} -v ny=${y[$i]} $2

  # convert the PXL files to PPM images
  gawk -f pxl2ppm.awk -v xd=${x[$i]} -v yd=${y[$i]} -v clr=0 -v max=$max1 -v out=${stem1[$i]}.ppm ${stem1[$i]}.pxl
  gawk -f pxl2ppm.awk -v xd=${x[$i]} -v yd=${y[$i]} -v clr=0 -v max=$max2 -v out=${stem2[$i]}.ppm ${stem2[$i]}.pxl
  gawk -f pxl2ppm.awk -v xd=${x[$i]} -v yd=${y[$i]} -v clr=0 -v max=$max3 -v out=${stem3[$i]}.ppm ${stem3[$i]}.pxl
  gawk -f pxl2ppm.awk -v xd=${x[$i]} -v yd=${y[$i]} -v clr=0 -v max=$max4 -v out=${stem4[$i]}.ppm ${stem4[$i]}.pxl
done
rm -f *.pxl

# Comment out the rest of the script if you do not have ImageMagick installed 
# (specifically, the convert program):
for i in ${icmd[*]}
do
  mkdir -p images
  convert ${stem1[$i]}.ppm images/${stem1[$i]}.png
  convert ${stem2[$i]}.ppm images/${stem2[$i]}.png
  convert ${stem3[$i]}.ppm images/${stem3[$i]}.png
  convert ${stem4[$i]}.ppm images/${stem4[$i]}.png
done
rm -f *.ppm

# Generate an HTML page displaying the images in a table
rm -f $1.html
gawk -f mkhtml.awk -v pre=$1 -v cmds=$cmdstring

# END
