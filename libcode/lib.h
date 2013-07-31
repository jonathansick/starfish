c *** Include file for mklib
c
c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
      integer NISO   ! number of isochrones in library [250]
      integer NPTS   ! maximum number of interpolated isoc points [20000]
      integer MAXDIM ! maximum dimension of fitted polynomial [20]
      integer MMAG   ! maximum number of magnitudes present in each isoc [8]
      real TOL       ! tolerance parameter for SVDCMP routine [1.0e-24]

      parameter(NISO=2000,NPTS=100000,MAXDIM=20,MMAG=8)
      parameter(TOL=1.0e-24)
