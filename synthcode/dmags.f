      subroutine dmags( icmd, ix, iy, rr, dx, dy )

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Select photometric error offsets from the crowding 
c     bin ix, iy in CMD icmd.
c
      include 'synth.h'

      integer icmd, ix, iy
      double precision dr
      real dx, dy, rr
      real x1, x2, f

      real emin,emax,epix
      common /c_chist/ emin,emax,epix

      real pdrop(MCMD,BINX,BINY)
      real cumerr(MCMD,BINX,BINY,NBINS*NBINS)
      common /c_cum/ pdrop,cumerr

      integer kk,kk0,nbin,iflag
      real xb(NBINS)
      save iflag,nbin,xb
      
      ! Only do these initializations on the first call
      if (iflag.eq.0) then
         iflag = 1
         nbin = nint((emax-emin)/epix)
         do i=1,nbin
            xb(i) = emin + epix*(i-0.5)
         end do
      endif
     
      ! find the index kk at which the cumfcn crosses the number rr
      dk = int(nbin*nbin/2)
      kk = dk

 10   c = cumerr(icmd,ix,iy,kk)

      dk = int(dk/2)
      if ( dk.le.0 ) goto 20
      if ( c.gt.rr ) then
         kk = kk - dk
      else
         kk = kk + dk
      endif
      goto 10

      write(*,*) 'Error: point not found in cumulative err function!'
      write(*,*) icmd,ix,iy,rr,cumerr(icmd,ix,iy,nbin*nbin)
      stop

      ! The exact value rr is crossed somewhere between 
      ! kk-1 and kk.  Interpolate to find the exact point.

      ! Convert from the 1-D cum.fcn. index (kk) to 
      ! the 2-D delta-mag pixel indices (kx,ky)
 20   ky = int( kk/nbin ) + 1
      kx = kk - nbin*( ky - 1 )

      ! determine the values of the bracketing cumfcn bins (x1 and x2)
      if ( kk.eq.1 ) then  
         x1 = 0.0
      else 
         x1 = cumerr(icmd,ix,iy,(kk-1))
      endif
      x2 = cumerr(icmd,ix,iy,kk)

      ! fractional distance of random number 
      ! between the bracketing cumfcn bins
      f = (rr - x1)/(x2 - x1)

      if ( kx.eq.1 ) then ! we can't interpolate because kk-1 and kk 
                          ! are on different rows 
         if ( f.lt.0.5 ) then ! need to extrapolate from kk-2 and kk-1
            kk = kk - 1
            ky = int( kk/nbin ) + 1
            kx = kk - nbin*( ky-1 )
            f = f + 1.0
         else           ! need to extrapolate back from kk and kk+1
            kk = kk + 1
            ky = int( kk/nbin ) + 1
            kx = kk - nbin*( ky-1 )
            f = f - 1.0
         endif
      endif

      ! the dX value is just interpolated from 
      ! the bracketing bin values
      dx = xb( kx-1 ) + f*( xb( kx ) - xb( kx-1 ) )
      ! the dY value is the row's bin value, plus random scatter.
      call random( dr, -0.50D0, 0.50D0 )
      dy = epix*real(dr) + xb( ky )

      return
      end

