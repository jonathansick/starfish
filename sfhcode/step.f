      subroutine step(p,pmin,size,flag)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     step: Take a step in N-dimensional parameter space.
C
C     flag=1: From position p(MP), take a step along unit vector pmin(MP).
C             size determines the magnitude of the step.
C     flag=0: From position p(MP), take a fractional step toward position
C             pmin(MP).  The fraction is given by size.
C     
      include 'sfh.h'

      integer j,nfree,flag
      double precision p(MP),pmin(MP)
      real size

      common /cnp/ np

      if (flag.eq.0) then
c        Take a fractional step from p to pmin
         do j=1,np
            p(j) = (1.0D0-dble(size))*p(j) + dble(size)*pmin(j)
         end do

      else if (flag.eq.1) then
c        Take a step in a random direction.  The direction cosines 
c        are contained in pmin.
         do j=1,np
            p(j) = p(j) + dble(size)*pmin(j)
            if (p(j).lt.0.0D0) p(j) = 0.0D0
         end do
      endif

      return
      end
