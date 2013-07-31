c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
      subroutine read_av(avfile,ntot,av)
c     Fill the extinction vector from the input file.
c
      include 'synth.h'

      character*40 avfile
      integer i,n,ntot
      real av(NAV)

      open(unit=4,file=avfile,status="old")
      n = 1

 20   read(4,*,end=30) av(n)
      n = n + 1
      goto 20

 30   close(4)
      ntot = n - 1

      return
      end


      subroutine redden(mag,age)
c     Adjust stellar photometry for interstellar extinction.
c
      include 'synth.h'

      integer i,redflag
      real mag(MMAG),age
      real av,ebv,q
      double precision p,xav

      integer nhav,ncav
      real hav(NAV),cav(NAV)
      common /c_av1/ nhav,ncav,hav,cav
      
      redflag = 0
      if (age.lt.7.0) then
         redflag = 1
      else 
         q = (age - 7.0)/2.0
         call random(p,0.0D0,1.0D0)
         if (real(p).gt.q) redflag = 1
      endif
      
      if (redflag.eq.1) then
         call random(xav,0.0D0,dble(nhav))
         i = int(xav) + 1
         call redden2(mag,hav(i))
      else
         call random(xav,0.0D0,dble(ncav))
         i = int(xav) + 1
         call redden2(mag,cav(i))
      endif
      
      return 
      end

      subroutine redden2(mag,av)

      include 'synth.h'

      real mag(MMAG),av

      real red(MMAG)
      common /c_av2/ red

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      do imag=1,nmag
         mag(imag) = mag(imag) + av*red(imag)
      end do

      return
      end
