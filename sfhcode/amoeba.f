      subroutine amoeba(p,y,iter,ftol)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     amoeba:  find the minimum of a function in MP dimensions by
c     following the local gradient.  
c
      
      include 'sfh.h'

      integer iter
      double precision p(MP+1,MP),y(MP+1)
      double precision rtol,ftol,ytry,ysave,amotry

      integer i,j,ihi,ilo,inhi,m,n,pnum(MP)
      double precision sum,swap,psum(MP),py(MP),psave(MP)

      integer np
      common /cnp/ np

      integer iverb
      common /verb/ iverb

      iter = 0
      
      do i=1,np
         py(i) = 0.D0
      end do

 1    do 12 n=1,np
         sum = 0.D0
         do 11 m = 1,np+1
            sum = sum + p(m,n)
 11      continue
         psum(n) = sum
 12   continue

 2    ilo = 1

      if (y(1).gt.y(2)) then
         ihi = 1
         inhi = 2
      else 
         ihi = 2
         inhi = 1
      endif

      do 13 i=1,np+1
         if (y(i).le.y(ilo)) ilo = i
         if (y(i).gt.y(ihi)) then
            inhi = ihi
            ihi = i
         else if (y(i).gt.y(inhi)) then
            if (i.ne.ihi) inhi = i
         endif
 13   continue

      rtol = 2.D0*abs(y(ihi) - y(ilo))/(abs(y(ihi)) + abs(y(ilo)))
      if (iverb.ge.5) write(*,*) rtol,ftol,y(ilo)

 20   if (rtol.lt.ftol) then
         swap = y(1)
         y(1) = y(ilo)
         y(ilo) = swap

         do 14 n=1,np
            swap = p(1,n)
            p(1,n) = p(ilo,n)
            p(ilo,n) = swap
 14      continue

         return
      endif

      iter = iter + 2

      ytry = amotry(p,y,psum,ihi,-1.0)

      if (ytry.le.y(ilo)) then
         if ( iverb.ge.3 ) 
     x        call dblemessage( "  current fitstat: ", ytry )

         ytry = amotry(p,y,psum,ihi,2.0)

      else if (ytry.ge.y(inhi)) then
         ysave = y(ihi)
         ytry = amotry(p,y,psum,ihi,0.5)
         if (ytry.ge.ysave) then
            do 16 i=1,np+1
               if (i.ne.ilo) then
                  do 15 j=1,np
                     psum(j) = 0.5*(p(i,j)+p(ilo,j))
                     if (psum(j).lt.0.D0) psum(j) = 0.005D0
                     p(i,j) = psum(j)
 15               continue
                  call fitstat(psum,y(i),0)
               endif
 16         continue
            iter = iter + np
            goto 1
         endif
      else
         iter = iter - 1
      endif
      
      goto 2

 6    format(f5.2,3g16.9,i4)
 7    format(f12.4,2g14.4,i4,f10.4)

      end

      function amotry(p,y,psum,ihi,fac)

      include 'sfh.h'

      integer ihi,j,np
      double precision p(MP+1,MP),psum(MP),sum
      double precision y(MP+1),ytry,amotry,ptry(MP)
      real fac,fac1,fac2

      common /cnp/ np
      fac1 = (1.0 - fac)/np
      fac2 = fac1 - fac

      do 11 j=1,np
         ptry(j) = psum(j)*fac1 - p(ihi,j)*fac2
         if (ptry(j).lt.0.D0) ptry(j) = 0.005D0
 11   continue
      
      call fitstat(ptry,ytry,0)
      
      if (ytry.lt.y(ihi)) then
         y(ihi) = ytry

         do 12 j=1,np
            psum(j) = psum(j) - p(ihi,j) + ptry(j)
            p(ihi,j) = ptry(j)
 12      continue
      endif

      amotry = ytry
      
      return
      end

