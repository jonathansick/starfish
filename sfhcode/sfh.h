c *** Include file for sfh
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
      integer MP    ! max. number of independent SFR amplitudes [55]
      integer NDATA ! max. number of input data stars [500000]
      integer MCMD  ! max. number of CMDs [7]
      integer MBOX  ! max. number of boxes per CMD for chi**2 [3000]

      parameter(MP=55,NDATA=500000,MCMD=7,MBOX=3000)
