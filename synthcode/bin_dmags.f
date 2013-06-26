      subroutine bin_dmags(icmd,ix,iy,dm)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Select photometric error offsets for the four 
c     crowding bins surrounding the current model star.
c     The photometric error offsets for each bracketing 
c     crowding bin are stored in dm(2,2,2).  The first
c     index is the dimension (1,2)==(color,magnitude).  
c     The other two indices identify the crowding bin by 
c     column and row.  e.g., dm(2,1,1) is the magnitude 
c     offset for the bin at position (1,1).
c
      include 'synth.h'

      integer icmd,ix,iy,jx,jy,j,kk,kk0(2,2)
      double precision dr
      real rr,dm(2,2,2)

      ! Choose a random number for indexing the cumulative 
      ! error functions.
      ! random returns a number less than 1.00, but casting 
      ! it as a real can round it up to 1.00.  Redraw the 
      ! number in this case.
 10   call random( dr, 0.0D0, 1.0D0 )
      rr = real(dr)
      if ( rr.eq.1.0 ) goto 10
      
      do jx=1,2     ! loop over the four 
         do jy=1,2  ! bracketing crowd bins
            iix = ix + (jx-1)  ! iix is either ix or ix+1
            iiy = iy + (jy-1)  ! iiy is either iy or iy+1

            ! select a delta-mag value from the cumulative fcn
            ! if each of the 4 bracketing crowding bins, using 
            ! the same random number rr to index each cumfcn.
            call dmags( icmd, iix, iiy, rr, dm(1,jx,jy), dm(2,jx,jy) )
         end do
      end do

      return
      end
