      program geteep

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c
c -=> StarFISH is free software; you can redistribute it
c -=> and/or modify it under the terms of the GNU General Public
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any
c -=> later version.
c
c     Determine EEP-like points along an isochrone empirically.
c     An Equivalent Evolutionary Point (EEP) is a flag indicating
c     a particular evolutionary state in an isochrone, such as
c     the main sequence turn-off or the tip of the red giant branch.
c
c     Such points are very useful for interpolating between isochrones
c     with adjacent values in age or metallicity, to determine an 
c     isochrone for any intermediate age or metallicity value.  However,
c     while EEPs used to be available for Padua isochrones in the form
c     of a table, this is no longer the case.  We must find a way to 
c     determine EEPs for ourselves.
c     
c     This program attempts to find empirical EEPs (eEEPs) by
c     identifying "nodes" along the isochrones, where a magnitude 
c     or color has an extremum.
c
c     Currently, we identify eEEPs in one CMD and use them to 
c     interpolate all of the photometry.

      parameter(NN=400,MMAG=8,NINTERVAL=4, DMAG=0.2)
      character*40 isofile, eepfile
      character*16 fstr(MMAG)
      real mass(NN),mag(MMAG,NN)
      real color(NN),dx(NN),dy(NN), strength(NN)
      integer i,j,n,nmag,lastnode,inode(NN)

      fstr(1) = "(f13.8,f8.3,i4)"
      fstr(2) = "(f13.8,2f8.3,i4)"
      fstr(3) = "(f13.8,3f8.3,i4)"
      fstr(4) = "(f13.8,4f8.3,i4)"
      fstr(5) = "(f13.8,5f8.3,i4)"
      fstr(6) = "(f13.8,6f8.3,i4)"
      fstr(7) = "(f13.8,7f8.3,i4)"
      fstr(8) = "(f13.8,8f8.3,i4)"

c     read in the input file (specified on the command line as STDIN):
      call fetchint( 5, nmag )

      if ( nmag.gt.MMAG ) then
         write(*,*) "Sorry, nmag is too large."
         write(*,*) nmag, MMAG
         stop
      endif
      call fetchint( 5, j0 )  ! the primary CMD mag

 10   call fetch2char40( 5, isofile, eepfile )
      if ( isofile.eq.'' ) goto 100  ! end-of-file
      do i=1,NN
         inode(i) = 0
         strength(i) = 0.0
      end do

      dx(1) = 0.0
      dy(1) = 0.0

      lastnode = 1
      i = 1
      open(unit=11,file=isofile,type="old")
 11   read(11,*,end=20) mass(i),(mag(j,i),j=1,nmag)
      color(i) = mag(j0-1,i) - mag(j0,i)

c     dx,dy are the change in color,mag from the previous point
      if ( i.gt.1 ) then
         dx(i) = color(i)  - color(i-1)
         dy(i) = mag(j0,i) - mag(j0,i-1)
      endif

c     check for an extremum in color or magnitude.  extremum defined 
c     as change in sign of dx or dy over its previous value
      if ( i.gt.2 ) then
         if ( dy(i).gt.0.0 .and. dy(i-1).lt.0.0 ) then ! i-1 was min
            inode(i-1) = 3
            strength(i-1) = dy(i) - dy(i-1)
         else if ( dy(i).lt.0.0 .and. dy(i-1).gt.0.0 ) then ! i-1 was max
            inode(i-1) = 4
            strength(i-1) = dy(i-1) - dy(i)
         else                   ! if ( dy(i).eq.0.0 ) then
            dy(i) = dy(i-1)     !keep last non-zero dy
         endif
         
         if ( dx(i).gt.0.0 .and. dx(i-1).lt.0.0 ) then ! i-1 was min
            inode(i-1) = 1
            strength(i-1) = 4*(dx(i) - dx(i-1))
         else if ( dx(i).lt.0.0 .and. dx(i-1).gt.0.0 ) then ! i-1 was max
            inode(i-1) = 2
            strength(i-1) = 4*(dx(i-1) - dx(i))
         else                   ! if ( dx(i).eq.0.0 ) then
            dx(i) = dx(i-1)     !keep last non-zero dx
         endif
            
c     Only keep nodes which are separated by more than DMAG 
c     in the CMD.  If the present node is within DMAG of the last
c     node, discard the weaker of the two.
         if ( inode(i-1).gt.0 ) then
            r = sqrt((color(i-1) - color(lastnode))**2 + 
     x           0.25*(mag(j0,i-1) - mag(j0,lastnode))**2 )
            
            if ( r .lt. DMAG ) then
               if ( strength(i-1) .gt. strength(lastnode) ) then
                  inode(lastnode) = 0
                  lastnode = i-1
               else 
                  inode(i-1) = 0
               endif
            else 
               lastnode = i-1
            endif
         endif

      endif

      i = i + 1
      goto 11

 20   close(11)
      n = i - 1

c     Write out eep points
      open(unit=12,file=eepfile,type="unknown")
      do i=1,n
c         if ( inode(i).gt.0 ) then
            write(12,fstr(nmag)) mass(i),(mag(j,i),j=1,nmag),inode(i)
c         endif
      end do
      close(12)

      goto 10

 100  stop
      end
