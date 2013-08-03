      subroutine direction(ex)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     Return the NP direction cosines in vector ex.
C     Numbers are assigned randomly to ex, and then these numbers 
C     are converted to direction cosines by dividing each element 
C     by the magnitude of the vector (SUM(ex(i)**2)).
C
      include 'sfh.h'

      integer j,np
      double precision ex(MP),r
      common /cnp/ np
      
 10   r = 0.0D0
      do j=1,np
         if (ex(j).eq.1.0D0) then ! only randomize significant amps
            call random(ex(j),-1.0D0,1.0D0)
            r = r + ex(j)**2
         endif
      end do

      r = sqrt(r)

      if (r.gt.0.0D0) then
         do j=1,np
            ex(j) = ex(j)/r
         end do
      else
         ! pause 'Warning: Null direction picked.  Selecting default' 
 
         ex(1) = 1.0D0
         do j=2,np
            ex(j) = 0.0D0
         end do
      endif

      return 
      end

