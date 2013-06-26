      subroutine cumfcn(dat,ndata,xmin,xmax,x,y,nbin,dbin)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     histogram the set of points dat(ndata) from xmin to xmax, with binsize
c     dbin.  Construct the cumulative distribution from the histogram,
c     return the cumulative distribution in y(nbin).
c
      include 'synth.h'

      integer i,nbin,ndata,n(NBINS)
      real dat(ndata),xmin,xmax,dbin
      real x(NBINS),y(NBINS)

      do i=1,nbin
         x(i) = xmin + dbin*(i - 0.5)
         y(i) = 0.0
         n(i) = 0
      end do

c     Lop through the data.  Identify each point's bin, and 
c     increment the bin's population. 
      do i = 1,ndata
         ibin = int((dat(i) - xmin)/dbin) + 1
         n(ibin) = n(ibin) + 1
      end do

      y(1) = 0.0
      do i=2,nbin
         y(i) = y(i-1) + real(n(i))/real(ndata)
      end do

      do i=1,nbin
         y(i) = y(i)/y(nbin)
      end do

      if (y(nbin).lt.1.0) then
         write(*,*) y(nbin),ndata
         pause
      endif

      return
      end
      
