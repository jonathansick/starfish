      subroutine outpix(nstars, frac, pix, pixsum, spsum)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Output the pixel array to three synthCMD files.
c     
      include 'synth.h'

      integer icmd, nstars
      integer ii,jj
      real frac
      double precision p,padd,pix(MCMD,MX,MY),pixsum(MCMD), spsum(MCMD)

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      integer nxpix(MCMD),nypix(MCMD)
      common /c_pix/ nxpix, nypix

      do icmd=1,ncmd
         pixsum(icmd) = 0.0
         do jj=nypix(icmd),1,-1
            do ii=1,nxpix(icmd)
               p = pix(icmd,ii,jj)*dble(frac)/dble(nstars)
               write(20+icmd,6) p

               pixsum(icmd) = pixsum(icmd) + p
               pix(icmd,ii,jj) = 0.0D0
            end do
         end do
c     calculate error on the pixsum.  Assume it's dominated by poisson noise
c     of padd (the sum of pix) and nstars
         padd = pixsum(icmd)*dble(nstars)/dble(frac)
         spsum(icmd) = padd**-1 + dble(nstars)**-1
      end do

 6    format(g14.7)

      return
      end
