#     This file is part of StarFISH
#     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
#     
# -=> StarFISH is free software; you can redistribute it 
# -=> and/or modify it under the terms of the GNU General Public 
# -=> License (GPL) as published by the Free Software Foundation;
# -=> either version 2 of the License, or (at your option) any 
# -=> later version.
#
# Sample data photometry script.  StarFISH expects a separate 
# photometry file for each CMD.  Often, you will have one single 
# photometry table.  Use something like this script to generate 
# individual CMD files.
#
($3>0 && $5>0) {
  printf("%7.3f %7.3f %9.6f %9.7f\n", $3-$5, $5, $1, $2) >> pre ".ub";
}
($5>0 && $7>0) {
  printf("%7.3f %7.3f %9.6f %9.7f\n", $5-$7, $7, $1, $2) >> pre ".bv";
}
($7>0 && $9>0) {
  printf("%7.3f %7.3f %9.6f %9.7f\n", $7-$9, $9, $1, $2) >> pre ".vi";
}
