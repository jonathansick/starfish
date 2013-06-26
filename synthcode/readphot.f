      subroutine readphot(libfile, n)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     read in the photometry and OP values for the current isochrone 
c     (libfile).
c
      include 'synth.h'

      character*40 libfile
      integer n,j
      real f

      integer vindex,imin
      real mag(MMAG,NISOS),vfaint,fbinary
      double precision mass(NISOS),prob(NISOS),pcum(NISOS)
      common /c_phot/ vindex,imin,vfaint,fbinary,mag,mass,prob,pcum

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      call oldfile(12,libfile)
      
      read(12,*) (mag(imag,1),imag=1,nmag),prob(1),mass(1)
      pcum(1) = 0.0
      imin = 1

      do j=2,NISOS
         read(12,*,end=10) (mag(imag,j),imag=1,nmag),
     x        prob(j),mass(j)
         pcum(j) = pcum(j-1) + prob(j)
         if (imin.eq.1.and.mag(vindex,j).lt.vfaint) imin = j
      end do
      
 10   close(12)
      n = j - 1

      return
      end
