      subroutine errs(x,y,yerr,n,stepsize,ppos,pneg)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     errs: given a best-fit set of amplitudes (x(MP)), compute the confidence 
C     interval on each amplitude by exploring the local fitstat topology.
C     Take small steps away from the minimum, recomputing fitstat at each 
C     location until y > yerr, where yerr is the fitstat value of the
C     desired confidence interval (default is 68%, see synth.dat input file).
C
      include 'sfh.h'

      integer i,j,n,istep,jper,jnext
      double precision ppos(MP),pneg(MP),pdir(MP)
      double precision x(MP),xtry(MP),xold(MP),xtest,xsum
      double precision y,ytry,yold,yerr,dy
      double precision updown(2)
      real stepsize

      integer np
      common /cnp/ np
      integer iverb
      common /verb/ iverb

      updown(1) = -1.0D0
      updown(2) =  1.0D0

      if (iverb.ge.2) write(*,*) ' Computing independent error ' //
     x     'of each amplitude... '
 
c     First, compute the uncorrelated confidence interval on each amplitude:
      do i=1,2   ! do both plus and minus directions
         do j=1,np
            istep = 0
            do k=1,np
               pdir(k) = 0.0D0
               if (k.eq.j) pdir(j) = updown(i)  ! +/- 1.0
               xtry(k) = x(k)
            end do
         
            ytry = y
            
c        step in this dir. until (ytry.gt.yerr)
 2          yold = ytry
            do k=1,np
               xold(k) = xtry(k)
            end do
            
            call step(xtry,pdir,stepsize,1)
            call fitstat(xtry,ytry,0)
            if (xtry(j).eq.0.0) ytry = yerr 

            istep = istep + 1
            if (ytry.lt.yerr) goto 2
         
c        xold was the last location below yerr.  Do a linear interpolation 
c        to estimate p(j) where y=yerr.
            dy = (yerr - yold)/(ytry - yold)
            xtest = xold(j) + dy*(xtry(j) - xold(j))
            if (xtest.gt.ppos(j)) ppos(j) = xtest
            if (xtest.lt.pneg(j)) pneg(j) = xtest
         end do
      end do

      if (iverb.ge.2) 
     x        write(*,*) ' Computing pairwise correlated errors...'

c     next, find correlated errors between adjacent isochrone pairs.
c     Give isochrone j a random diection cosine, and isochrone j+1 its 
c     negative complement.
      do i=1,100
         if ( iverb.ge.3 ) write(*,*) '  ', i, '% complete'
         do j=1,np-1
            istep = 0

            do k=1,np-1
               xtry(k) = x(k)

               if (k.eq.j) then
                  call random(pdir(k),-1.0D0,1.0D0) 
                  pdir(k+1) = 1.0D0 - abs(pdir(k))
                  if (pdir(k).gt.0.0D0) pdir(k+1) = -1.0D0*pdir(k+1)
               else
                  pdir(k) = 0.0D0
               endif
            end do
            
            ytry = y
            
c        step in this dir. until (ytry.gt.yerr)
 3          yold = ytry
            do k=1,np
               xold(k) = xtry(k)
            end do
            
            call step(xtry,pdir,stepsize,1)
            call fitstat(xtry,ytry,0)
            if (xtry(j).le.0.0) ytry = yerr 

            istep = istep + 1
			if (istep.gt.100) write(*,*) istep,ytry,yerr
            if (ytry.lt.yerr) goto 3
         
c        xold was the last location below yerr.  Do a linear interpolation 
c        to estimate p(j) where y=yerr.
            dy = (yerr - yold)/(ytry - yold)
            xtest = xold(j) + dy*(xtry(j) - xold(j))
            if (xtest.gt.ppos(j)) ppos(j) = xtest
            if (xtest.lt.pneg(j)) pneg(j) = xtest
         end do
      end do


c      if (iverb.ge.3) write(*,*) 'Checking active isochrones...'
c
cc     next, find correlated errors between all "significant" isochrones
c      xsum = 0.0D0
c      do j=1,np
c         xsum = xsum + x(j)
c      end do
c
c      jper = 10
c      jnext = n/10
c
c      do i=1,n
c         if (iverb.ge.3) then
c            if (i.ge.jnext) then
c               write(*,*) jper, '% complete'
c               jper = jper + 10
c               jnext = jnext + n/10
c            endif
c         endif
c            
c         do j=1,np ! randomize significant amplitudes
c            xtry(j) = x(j)
c
c            pdir(j) = 0.0D0
c            if (x(j)/xsum.gt.0.01) pdir(j) = 1.0D0
c         end do
c         call direction(pdir)
c
c         istep = 0
c         ytry = y
c
cc        step in this dir. until (ytry.gt.yerr)
c 4       yold = ytry
c         do j=1,np
c            xold(j) = xtry(j)
c         end do
c
c         call step(xtry,pdir,1000.0,1)
c         call fitstat(xtry,ytry,0)
c         istep = istep + 1
c         if (ytry.lt.yerr) goto 4
c
cc        xold was the last location below yerr.  Do a linear interpolation 
cc        to estimate p(j) where y=yerr.
c         write(*,*) istep
c         do j=1,np
c            dy = (yerr - yold)/(ytry - yold)
c            xtest = xold(j) + dy*(xtry(j) - xold(j))
c            if (xtest.gt.ppos(j)) ppos(j) = xtest
c            if (xtest.lt.pneg(j)) pneg(j) = xtest
c         end do
c      end do
c

      if ( iverb.ge.2 ) write(*,*) 
     x     ' Computing correlated errors among all amps... '

c     finally, compute general correlated confidence intervals by exploring 
c     random parameter directions.
      jper = 10
      jnext = n/10

      do i=1,n
         if (iverb.ge.3) then
            if (i.ge.jnext) then
               write(*,*) jper, '% complete'
               jper = jper + 10
               jnext = jnext + n/10
            endif
         endif
         do j=1,np ! randomize every amplitude in making new pdir vector
            xtry(j) = x(j)
            pdir(j) = 1.0D0
         end do
         call direction(pdir)

         istep = 0
         ytry = y

c        step in this dir. until (ytry.gt.yerr)
 10      yold = ytry
         do j=1,np
            xold(j) = xtry(j)
         end do

         call step(xtry,pdir,stepsize,1)
         call fitstat(xtry,ytry,0)
         istep = istep + 1
         if (ytry.lt.yerr) goto 10

c        xold was the last location below yerr.  Do a linear interpolation 
c        to estimate p(j) where y=yerr.
         do j=1,np
            dy = (yerr - yold)/(ytry - yold)
            xtest = xold(j) + dy*(xtry(j) - xold(j))
            if (xtest.gt.ppos(j)) then
               ppos(j) = xtest
               write(*,*) 'max: ',j,istep,ppos(j),pdir(j)
            endif

            if (xtest.lt.pneg(j)) then
               pneg(j) = xtest
               write(*,*) 'min: ',j,istep,pneg(j),pdir(j)
            endif
         end do
      end do

      return
      end







