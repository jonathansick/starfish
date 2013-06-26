      subroutine dfstat( p, xi )

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     Given a point in multidimensional parameter space p,
C     return the multidimensional gradient of the fit-statistic 
C     function at that location (xi).

      include 'sfh.h'

      double precision p(MP), xi(MP), dfdm(MCMD,MBOX)
      
      integer np
      common /cnp/ np

      real dat(MCMD,2,NDATA),chem(MP),age(MP)
      integer fstat,gtype
      integer nstars(MCMD),pnum(MP),mask(MCMD,MBOX)
      character*40 cmdfile,chifile
      common /cstat/ dat,cmdfile,chifile,nstars,mask,
     x               chem,age,pnum,fstat,gtype

      double precision psynth(MCMD,MP,MBOX)
      double precision pmod(MCMD,MBOX),pdat(MCMD,MBOX)
      common /cfit/ psynth, pmod, pdat

      character*8 suffix(MCMD)
      integer npix, ncmd, nbox(MCMD)
      real xlo(MCMD), xhi(MCMD), ylo(MCMD), yhi(MCMD), dpix
      common /cmds/ suffix, xlo, xhi, ylo, yhi, nbox, dpix, npix, ncmd

C     Initialize xi
      do i=1,np
         xi(i) = 0.0D0
      end do

C     dfdm is the derivative of the fit statistic w.r.t.
C     the number of model stars in each CMD box.  Calculate 
C     this one time for each box.
      do icmd=1,ncmd
         do ibox=1,nbox(icmd)
            
            if ( fstat .eq. 2 ) then !Poisson
               if ( pmod(icmd,ibox).eq.0.0D0 .and. 
     x              pdat(icmd,ibox).eq.0.0D0 ) then
                  dfdm(icmd,ibox) = 0.0D0
               else if ( pmod(icmd,ibox).eq.0.0D0 ) then
c                 if pmod=0, use pmod=1.0, which is "close to zero"
                  dfdm(icmd,ibox) = 1.0D0 - 
     x                 pdat(icmd,ibox)/1.0D0
               else
                  dfdm(icmd,ibox) = 1.0D0 - 
     x                 pdat(icmd,ibox)/pmod(icmd,ibox)
               endif
            else   !chi-squared or Lorentzian
               if ( pmod(icmd,ibox).eq.0.0D0 .and. 
     x              pdat(icmd,ibox).eq.0.0D0 ) then
                  dfdm(icmd,ibox) = 0.0D0
               else if ( pdat(icmd,ibox).eq.0.0D0 ) then
c                 if pdat=0, use pdat=0.01, which is "close to zero"
                  dfdm(icmd,ibox) = pmod(icmd,ibox)/0.01D0 - 1.0D0
               else
                  dfdm(icmd,ibox) = 
     x                 pmod(icmd,ibox)/pdat(icmd,ibox) - 1.0D0
               endif
            endif
         end do 
      end do
      
      do i=1,np
         do icmd=1,ncmd
            do ibox=1,nbox(icmd)
               xi(i) = xi(i) + dfdm(icmd,ibox)*psynth(icmd,i,ibox)
            end do
         end do

         if ( fstat.ne.1 ) then  !Poisson or Chi-squared
            xi(i) = 2.0D0*xi(i)
         endif
      end do

      return
      end

