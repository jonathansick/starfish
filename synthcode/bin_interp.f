      subroutine bin_interp(icmd,ix,iy,dxbin,dybin,dm,pdrop,dx,dy,p0)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Bilinear interpolation to determine the dropout rate (p0) 
c     and photometric error offsets (dx,dy) at the position of the 
c     current star in the current CMD (icmd) given the dropout rates
c     (pdrop) and stochastically-selected delta-mag values from the 
c     four bracketing  crowd bins (dm), and the relative position of 
c     the star between the bins (dxbin,dybin).
c
      include 'synth.h'

      integer icmd,ix,iy
      real dxbin,dybin
      real dm(2,2,2),pdrop(MCMD,BINX,BINY)
      real p0,dpx,dpy
      real dx,dy

      
      p0 = pdrop(icmd,ix,iy) +
     x     dxbin * ( pdrop(icmd,ix+1,iy) - pdrop(icmd,ix,iy) ) +
     x     dybin * ( pdrop(icmd,ix,iy+1) - pdrop(icmd,ix,iy) ) +
     x     dxbin*dybin * ( pdrop(icmd,ix,iy) + pdrop(icmd,ix+1,iy+1) -
     x                     pdrop(icmd,ix+1,iy) - pdrop(icmd,ix,iy+1) )

c     Some of the bracketting bins may have 100% dropout; we should 
c     not interpolate from these bins.  As a workaround, I will
c     duplicate the dm value from a neighboring bin.
      !! ix,iy
      if ( pdrop(icmd,ix,iy).eq.1.0 ) then
         if ( pdrop(icmd,ix+1,iy).ne.1.0 ) then
            dm(1,1,1) = dm(1,2,1)
            dm(2,1,1) = dm(2,2,1)
         else if ( pdrop(icmd,ix,iy+1).ne.1.0 ) then
            dm(1,1,1) = dm(1,1,2)
            dm(2,1,1) = dm(2,1,2)
         else if ( pdrop(icmd,ix+1,iy+1).ne.1.0 ) then
            dm(1,1,1) = dm(1,2,2)
            dm(2,1,1) = dm(2,2,2)
         endif
      endif

      !! ix+1,iy
      if ( pdrop(icmd,ix+1,iy).eq.1.0 ) then
         if ( pdrop(icmd,ix,iy).ne.1.0 ) then
            dm(1,2,1) = dm(1,1,1)
            dm(2,2,1) = dm(2,1,1)
         else if ( pdrop(icmd,ix+1,iy+1).ne.1.0 ) then
            dm(1,2,1) = dm(1,2,2)
            dm(2,2,1) = dm(2,2,2)
         else if ( pdrop(icmd,ix,iy+1).ne.1.0 ) then
            dm(1,2,1) = dm(1,1,2)
            dm(2,2,1) = dm(2,1,2)
         endif
      endif

      !! ix,iy+1
      if ( pdrop(icmd,ix,iy+1).eq.1.0 ) then
         if ( pdrop(icmd,ix+1,iy+1).ne.1.0 ) then
            dm(1,1,2) = dm(1,2,2)
            dm(2,1,2) = dm(2,2,2)
         else if ( pdrop(icmd,ix,iy).ne.1.0 ) then
            dm(1,1,2) = dm(1,1,1)
            dm(2,1,2) = dm(2,1,1)
         else if ( pdrop(icmd,ix+1,iy).ne.1.0 ) then
            dm(1,1,2) = dm(1,2,1)
            dm(2,1,2) = dm(2,2,1)
         endif
      endif

      !! ix+1,iy+1
      if ( pdrop(icmd,ix+1,iy+1).eq.1.0 ) then
         if ( pdrop(icmd,ix,iy+1).ne.1.0 ) then
            dm(1,2,2) = dm(1,1,2)
            dm(2,2,2) = dm(2,1,2)
         else if ( pdrop(icmd,ix+1,iy).ne.1.0 ) then
            dm(1,2,2) = dm(1,2,1)
            dm(2,2,2) = dm(2,2,1)
         else if ( pdrop(icmd,ix,iy).ne.1.0 ) then
            dm(1,2,2) = dm(1,1,1)
            dm(2,2,2) = dm(2,1,1)
         endif
      endif

      dx = dm(1,1,1) +
     x     dxbin * ( dm(1,2,1) - dm(1,1,1) ) +
     x     dybin * ( dm(1,1,2) - dm(1,1,1) ) +
     x     dxbin*dybin * ( dm(1,1,1) + dm(1,2,2) -
     x                     dm(1,2,1) - dm(1,1,2) )

      dy = dm(2,1,1) +
     x     dxbin * ( dm(2,2,1) - dm(2,1,1) ) +
     x     dybin * ( dm(2,1,2) - dm(2,1,1) ) +
     x     dxbin*dybin * ( dm(2,1,1) + dm(2,2,2) -
     x                     dm(2,2,1) - dm(2,1,2) )

      if (p0.lt.0.0) p0 = 0.0
      if (p0.gt.1.0) p0 = 1.0

      return
      end
