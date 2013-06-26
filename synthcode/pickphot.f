      subroutine pickphot(starmag, np)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Draw stellar photometry randomly from the current isochrone file, 
c     according to its OP distribution.  A binary fraction is implemented: 
c     For some fraction of the stars, another set of photometry is drawn and 
c     flux-added to the original.  The combined magnitudes are returned as a 
c     single detected object.
c
      include 'synth.h'

      integer k,k1,k2,ibinary,np,ix
      real f1,f2,mstar,q
      real mag1(MMAG),mag2(MMAG),starmag(MMAG)
      
      integer vindex,imin0
      real mag(MMAG,NISOS),vfaint,fbinary
      double precision mass(NISOS),prob(NISOS),pcum(NISOS),prand,p
      common /c_phot/ vindex,imin0,vfaint,fbinary,mag,mass,prob,pcum

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      ibinary = 0
      imin = imin0

 11   if (ibinary.eq.1) imin = 1  ! the second star can be very faint
      call random(prand,pcum(imin),pcum(np))

c     Find the isochrone points which bracket this probability:
 12   do k=imin,np
         if (pcum(k).ge.prand) then
            k2 = k
            k1 = k-1
            goto 14
         endif
      end do
      
      write(*,*) 'Pickphot: Random draw outside prob. interval.'
      pause

c     interpolate to find the intrinsic photometry.
 14   if (ibinary.eq.0) then
         f2 = (prand - pcum(k1))/(pcum(k2) - pcum(k1))
         f1 = 1 - f2

         do imag=1,nmag
            starmag(imag) = f1*mag(imag,k1) + f2*mag(imag,k2)
         end do
         mstar = f1*mass(k1) + f2*mass(k2)

         call random(p,0.0D0,1.0D0)
         if (real(p).lt.fbinary) then
            ibinary = 1
            goto 11
         endif
      else
         f2 = real((prand - pcum(k1))/(pcum(k2) - pcum(k1)))
         f1 = 1 - f2
         
         do imag=1,nmag
            mag2(imag) = f1*mag(imag,k1) + f2*mag(imag,k2)
            mag1(imag) = starmag(imag)

c           Flux-add the two components:
            starmag(imag) = -2.5*alog10(10**(-0.4*mag1(imag)) 
     x                                + 10**(-0.4*mag2(imag)))
         end do
         mstar = mstar + f1*mass(k1) + f2*mass(k2)

      endif
      
      return
      end

