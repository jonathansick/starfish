      program synth

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Construct a set of synthetic CMD Hess diagrams from the isochrone
c     library constructed by mklib.  Each synthetic CMD represents the 
c     predicted probability distribution of stars from a given epoch, 
c     under given seeing and reddening conditions.
c     The seeing conditions are determined by artificial star tests;
c     The reddening conditions are determined by examining the data (using
c     either reddening-free parameters on the OB population or 
c     blackbody fitting to the hot and cool stellar populations).
c     The CMDs are constructed by randomly placing points
c     in the CMDs, according to the OP's found by mklib.  Then, reddening 
c     and photometric error offsets are applied.
c
      include 'synth.h'

      character*40 isofile,infile,libfile1,libfile2
      character*40 hotfile,coldfile,lockfile,outfile(MCMD)
      character*8 postfix(MCMD)
      character*40 outstem
      integer i,redflag,ib,oldb,idrop(MCMD)
      integer jnext,jpercent
      integer nstars,nadd,np
      integer interp_errs, err_method
      real x(MCMD),y(MCMD),age,dpix
      real xmax(MCMD),ycmax(MCMD),ypmax(MCMD)
      real gamma,mass1,mass2,mrange,fsum,frac
      real starmag(MMAG), vmsto
      double precision pix(MCMD,MX,MY), pixsum(MCMD), spsum(MCMD)

      integer vindex,imin
      real mag(MMAG,NISOS),vfaint,fbinary
      double precision mass(NISOS),prob(NISOS),pcum(NISOS)
      common /c_phot/ vindex,imin,vfaint,fbinary,mag,mass,prob,pcum

      integer nxpix(MCMD),nypix(MCMD)
      common /c_pix/ nxpix, nypix

      integer nhav,ncav
      real hav(NAV),cav(NAV)
      common /c_av1/ nhav,ncav,hav,cav

      real red(MMAG)
      common /c_av2/ red

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      character*40 crowd1,crowd2
      common /c_cfile/ crowd1,crowd2
      
      integer nx(MCMD),ny(MCMD)
      real dbinx,dbiny,xmin(MCMD),ymin(MCMD)
      common /c_crowd/ dbinx,dbiny,nx,ny,xmin,ymin
      
      real emin,emax,epix
      common /c_chist/ emin,emax,epix

      character*8 xeqn(MCMD),yeqn(MCMD)
      common /c_eqn/ xeqn,yeqn

      integer iverb
      common /c_control/ iverb

      real pdrop(MCMD,BINX,BINY)
      real cumerr(MCMD,BINX,BINY,NBINS*NBINS)
      common /c_cum/ pdrop,cumerr

c *** Filenames ***
      call fetchchar40( 5, isofile )
      call fetchchar40( 5, lockfile )
      call fetchchar40( 5, hotfile )
      call fetchchar40( 5, coldfile )
      call fetchchar40( 5, crowd1 )
      call fetchchar40( 5, crowd2 )

c *** CMD limits ***
      call fetchint( 5, nmag )
      call fetchint( 5, ncmd )
      call fetchint( 5, vindex )
      call fetchreal( 5, dpix )

      do i=1,ncmd
         call fetchchar8( 5, xeqn(i) )
         call fetchchar8( 5, yeqn(i) )
         call fetchreal( 5, xmin(i) )
         call fetchreal( 5, xmax(i) )
         call fetchreal( 5, ymin(i) )
         call fetchreal( 5, ycmax(i) )
         call fetchreal( 5, ypmax(i) )
         call fetchchar8( 5, postfix(i) )
      end do

c *** Crowding parameters ***
      call fetchreal( 5, dbinx )
      call fetchreal( 5, dbiny )
      call fetchreal( 5, emin )
      call fetchreal( 5, emax )
      call fetchreal( 5, epix )

c *** Reddening parameters ***
      do i=1,nmag
         call fetchreal( 5, red(i) )
      end do

c *** Miscellaneous ***
      call fetchint( 5, iverb )
      call fetchint( 5, interp_errs )
      call fetchint( 5, err_method )
      call fetchint( 5, nscale )
      call fetchint( 5, iseed )
      call fetchreal( 5, mass1 )
      call fetchreal( 5, mass2 )
      call fetchreal( 5, gamma )
      call fetchreal( 5, vfaint )
      call fetchreal( 5, fbinary )

      niso = 1
      nhav = 0
      ncav = 0
      oldb = 0
      nstars = 0
      mrange = mass2**gamma - mass1**gamma
      fsum = 0.0
 
      ! Compute the number of crowding bins and number of pixels 
      ! in each dimension of each CMD.
      do icmd=1,ncmd
         nx(icmd) = nint((xmax(icmd)-xmin(icmd))/dbinx)
         ny(icmd) = nint((ycmax(icmd)-ymin(icmd))/dbiny)
         nxpix(icmd) = nint((xmax(icmd)-xmin(icmd))/dpix) 
         nypix(icmd) = nint((ypmax(icmd)-ymin(icmd))/dpix) 
         pixsum(icmd) = 0.0

         do jj=1,nypix(icmd)
            do ii=1,nxpix(icmd)
               pix(icmd,ii,jj) = 0.0D0
            end do
         end do
      end do

      ! Initialize random number generator
      call sdprnd(iseed)

      ! Fill reddening vectors
      call read_av(hotfile,nhav,hav)
      call read_av(coldfile,ncav,cav)

      call oldfile(10,isofile)
      call oldfile(11,lockfile)

      ! Loop through isochrones
      do 30 i=1,NISOC
         call readisodatline( 10, age, infile, libfile1, vmsto )
         if ( age.lt.0.0 ) goto 99  ! end-of-file reached

         call readisolockline( 11, ib, libfile2, outstem )
         if ( ib.lt.0 ) goto 99  ! end-of-file reached

         if (libfile1.ne.libfile2) then
            write(*,*) ' Error: mismatch between isofile and lockfile '
            write(*,*) libfile1
            write(*,*) libfile2
            pause
         endif

         if (ib.ne.oldb) then
            ! We have finished with the current synthCMD.
            ! output the pixels and prepare for the next synthCMD
            if (oldb.gt.0) then  !only output if not the first iso group
               frac = fsum/niso
               call outpix(nstars, frac, pix, pixsum,spsum)

               nstars = 0
               fsum = 0.0

               do icmd=1,ncmd
                  close(20+i)
               end do
            endif

            do icmd=1,ncmd
               call strcat(outstem, postfix(icmd), outfile(icmd))
               call newfile(20+icmd,outfile(icmd))
               
               pixsum(icmd) = 0.0
            end do

            niso = 1
            oldb = ib
         else
            niso = niso + 1
         endif

         if (iverb.gt.0) write(*,*) 'Reading ',libfile2

         ! Read in the UBVI photometry, the masses and the relative 
         ! occupation probabilities for points in the current isochrone.
         call readphot(libfile2, np)

         fsum = fsum + (mass(np)**gamma - mass(imin)**gamma)/mrange
         nn = nscale

         jnext = nn/10
         jpercent = 10
         do 20 j=1,nn
            if (j.gt.jnext) then
               if ( iverb.gt.0) write(*,*) jpercent, '% complete'
               jnext = jnext + nn/10
               jpercent = jpercent + 10
            endif

            ! choose stellar photometry randomly from the isochrone
            call pickphot(starmag, np)
            nstars = nstars + 1

            ! Choose a random extinction value from an extinction vector.
            call redden(starmag, age)

            ! Combine magnitudes according to xeqn(icmd), yeqn(icmd) 
            ! to form CMD axes
            call axes(starmag,x,y)
            
            ! Add photometric scatter 
            if ( err_method .eq. 2 ) then
               call scatter( x, y, idrop )
            else if ( err_method .eq. 1 ) then
               call fakecrowd( x, y, idrop )
            else 
               call crowd( x, y, idrop )
            endif

            ! Increment pixel values
            do icmd=1,ncmd
               ix = int((x(icmd) - xmin(icmd))/dpix) + 1
               iy = int((y(icmd) - ymin(icmd))/dpix) + 1

               if (idrop(icmd).eq.0) then
                  if (ix.ge.1.and.ix.le.nxpix(icmd).and.
     x                 iy.ge.1.and.iy.le.nypix(icmd)) then
                     pix(icmd,ix,iy) = pix(icmd,ix,iy) + 1.0D0
                  endif
               endif
            end do
            
 20      continue
 30   continue
      
      ! Write out final pixel file
 99   frac = fsum/niso
      call outpix(nstars, frac, pix, pixsum, spsum)

      do icmd=1,ncmd
         close(20+icmd)
      end do

      stop
      end
