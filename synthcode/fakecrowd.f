      subroutine fakecrowd( x, y, dropflag )

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Add photometric errors according to an analytic error model

      include 'synth.h'

      integer i,j,dropflag(MCMD)
      double precision xtest,ytest,ptest
      real x(MCMD),y(MCMD)
      real mag0(MMAG),dmag(MMAG)
      real sigma, sigma0, pdrop

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      double precision xran(105)
      common /rnd/ xran

      common /sig0/ sigma0

      ! The minimum photometric error for bright stars
      sigma0 = 0.02

      ! need magnitudes from the CMD axes; need to modify by hand (sorry!)
      mag0(1) = x(1) + y(1)  ! U mag
      mag0(2) = y(1)         ! B mag
      mag0(3) = y(2)         ! V mag
      mag0(4) = y(3)         ! I mag

      do i=1,nmag
         ! Determine dropout fraction for this mag
         call getdropfrac( pdrop, mag0(i) )
         call random( ptest, 0.0D0, 1.0D0 )

         if ( real(ptest).ge.pdrop ) then ! not a dropout

            ! draw random number from a unit Gaussian distribution
 11         call random(ytest, 0.0D0, 1.0D0)
            call random(xtest, -5.0D0, 5.0D0)
            if ( (exp(-0.5*xtest**2)).lt.ytest ) goto 11
               
            ! determine phot. err. sigma for this star
            call getsigma( sigma, mag0(i) )
             
            dmag(i) = real(xtest)*sigma ! scale random number by sigma
         else 
            dmag(i) = 9.999
         endif
      end do
         
      do i=1,ncmd
         if ( dmag(i).lt.9.0 .and. dmag(i+1).lt.9.0 ) then
            dropflag(i) = 0
            x(i) = x(i) + dmag(i) - dmag(i+1)
            y(i) = y(i) + dmag(i+1)
         else
            dropflag(i) = 1
         endif
      end do

      return
      end


      subroutine getsigma( sigma, mag )

c     Determine mag-dependent sigma.  sigma0 is the minimum error 
c     for bright mags.  sigma increases exponentially; the parameters
c     are set by hand to produce a reasonable curve.

      real amp, mag, mag0, sigma, sigma0

      common /sig0/ sigma0

c     amp characterizes the scale of the scatter
      amp = 0.06

c     mag0 is the mag where sigma increases to sigma0 + amp.
      mag0 = 22.0

      sigma = sigma0 + amp * exp( 0.5*(mag - mag0) )

      return
      end


      subroutine getdropfrac( pdrop, mag )

c     determine mag-dependent dropout rate.  The function and its params
c     are hand-selected to get a reasonable shape.
c     
c     We use an arctan function centered on the 50% completeness mag.
c     This function looks quite similar to real completeness functions.
c
      real mag, pdrop, mag0, s

c     mag0 is the mag at which the completeness falls to 50%
      mag0 = 22.0

c     s characterizes how steep the dropout "shelf" is in the 
c     neighborhood of mag0
      s = 0.325

      if ( mag.lt.mag0 ) then
         pdrop = 0.5 + s * atan( 2.0*(mag-mag0) )
      else
         pdrop = 0.5 + s * atan( 4.0*(mag-mag0) )
      endif

      if ( pdrop.gt.0.95 ) pdrop = 1.0
      if ( pdrop.lt.0.0 ) pdrop = 0.0

      return
      end
