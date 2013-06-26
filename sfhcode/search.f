      subroutine search(p,y,n,stepsize)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     search: look for lower fitstat values by stepping along random 
C     directions.
C
C     Choose a random parameter space direction, and take a small step.  
C     Compute fitstat for the new location.  If it is smaller than the last
C     location's fitstat, step again.  Repeat until the next step produces
C     a larger fitstat than the previous step.  Repeat this procedure for 
C     n random directions, and return the parameter space location 
C     with the lowest fitstat value in p(MP)
C
C     This procedure complements the fitstat minimization
C     algorithm, and is run first to get a "ballpark" solution that will
C     converge quickly.  See the User Manual for more information.
C
      include 'sfh.h'

      integer i,j,n,imask(MP),jper,jnext
      real stepsize
      double precision p(MP),p0(MP),pold(MP),pdir(MP)
      double precision y,y0,yold,psum

      integer np
      common /cnp/ np

      integer iverb
      common /verb/ iverb

      y0 = y
      psum = 0.0D0

      do j=1,np
         p0(j)=p(j)
         psum = psum + p0(j)
      end do

      do j=1,np
         imask(j) = 1
         if (p0(j)/psum.gt.0.01) imask(j) = 0 ! unmask significant amps
      end do

      jper = 10
      jnext = n/10

      do i=1,n
         if (iverb.ge.3) then
            if (i.ge.jnext) then
               write(*,*) '  ', jper, '% complete'
               jper = jper + 10
               jnext = jnext + n/10
            endif
         endif

         do j=1,np
            pdir(j) = 1.0D0  ! default amp state is significant
            if (imask(j).eq.1) then
c               call random(xx,0.0,1.0)
               pdir(j) = 0.0D0  ! flag amp as insignificant
            endif
         end do

         call direction(pdir)
		 
 10      yold = y0
         do j=1,np
            pold(j) = p0(j)
         end do

         call step(p0,pdir,stepsize,1)
         call fitstat(p0,y0,0)

         if (y0.lt.yold) goto 10

c     Do we have a new minimum?
         if (yold.lt.y) then
            y = yold
            do j=1,np
               p(j) = pold(j)
            end do
         endif
      end do

      return
      end
