c *** Include file for synth
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
      integer MMAG   ! max. number of magnitudes [8]
      integer MCMD   ! max. number of CMDs [7]
      integer NISOC  ! max. number of isochrones in library [500]
      integer NISOS  ! max. number of stars per isochrone [20000]
      integer NAV    ! max. number of extinction (Av) measurements [5000]
      integer MX     ! max. number of Hess pixels in X-direction [200]
      integer MY     ! max. number of Hess pixels in Y-direction [600]
      integer BINX   ! max. number of crowding bins in X [10]
      integer BINY   ! max. number of crowding bins in Y [52]
      integer MAXSTR ! max. number of artificial stars [1000000]
      integer MINPOP ! min. number of fake stars for non-empty crowd bin [10]
      integer MAXPOP ! max. number of fake stars in a cowd bin [20000]
      integer NBINS  ! number of delta-mag bins in the error histograms [50]
      real PI        ! acos(-1.0)
      
      parameter(MMAG=8,MCMD=7,NISOC=2000,NISOS=100000,NAV=5000)
      parameter(MX=200,MY=600,BINX=10,BINY=52)
      parameter(MAXSTR=1000000,MINPOP=10,MAXPOP=20000,NBINS=50)
      parameter(PI=3.141592654)
