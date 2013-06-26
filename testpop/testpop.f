      program testpop

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Constructs artificial stellar photometry based on an input SFH,
c     isochrones, and data-derived photometric conditions (reddening and 
c     crowding errors).

      include '../synthcode/synth.h'

      character*40 agefile,isofile,libfile,infile
      character*40 hotfile,coldfile,lockfile,outfile(MCMD)
      character*40 stem1,stem
      character*8  pre,postfix(MCMD)
      integer lockflag,sfrflag,i,ib,oldb,idrop(MCMD)
      integer jnext,jpercent,iseed
      integer nstars,np
      integer err_method, interp_errs
      real age,t,z,t2,z2,dtime
      real x(MCMD),y(MCMD),x0(MCMD),y0(MCMD)
      real xmax(MCMD),ycmax(MCMD),ypmax(MCMD)
      real gamma,mass1,mass2,mrange,frac
      real dmod,avfactor,flum(NISOS)
      real starmag(MMAG), dpix

      integer vindex,imin
      real mag(MMAG,NISOS),vfaint,fbinary
      double precision mass(NISOS),prob(NISOS),pcum(NISOS)
      common /c_phot/ vindex,imin,vfaint,fbinary,mag,mass,prob,pcum

      integer nxpix(MCMD),nypix(MCMD)
      common /c_pix/ nxpix, nypix

      integer namg,ncmd
      common /c_nmag/ nmag,ncmd

      integer nhav,ncav
      real hav(NAV),cav(NAV)
      common /c_av1/ nhav,ncav,hav,cav

      real red(MMAG)
      common /c_av2/ red

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

      call oldfile(1,"synth.dat")
c *** Filenames ***
      call fetchchar40( 1, agefile )
      call fetchchar40( 1, lockfile )
      call fetchchar40( 1, hotfile )
      call fetchchar40( 1, coldfile )
      call fetchchar40( 1, crowd1 )
      call fetchchar40( 1, crowd2 )

      agefile  = "../" // agefile
      lockfile = "../" // lockfile
      hotfile  = "../" // hotfile
      coldfile = "../" // coldfile
      crowd1   = "../" // crowd1
      crowd2   = "../" // crowd2

c *** CMD limits ***
      call fetchint( 1, nmag )
      call fetchint( 1, ncmd )
      call fetchint( 1, vindex )
      call fetchreal( 1, dpix )

      do i=1,ncmd
         call fetchchar8( 1, xeqn(i) )
         call fetchchar8( 1, yeqn(i) )
         call fetchreal( 1, xmin(i) )
         call fetchreal( 1, xmax(i) )
         call fetchreal( 1, ymin(i) )
         call fetchreal( 1, ycmax(i) )
         call fetchreal( 1, ypmax(i) )
         call fetchchar8( 1, postfix(i) )
      end do

c *** Crowding parameters ***
      call fetchreal( 1, dbinx )
      call fetchreal( 1, dbiny )
      call fetchreal( 1, emin )
      call fetchreal( 1, emax )
      call fetchreal( 1, epix )

c *** Reddening parameters ***
      do i=1,nmag
         call fetchreal( 1, red(i) )
      end do

c *** Miscellaneous ***
      call fetchint( 1, iverb )
      call fetchint( 1, interp_errs )
      call fetchint( 1, err_method )
      call fetchint( 1, nscale )
      call fetchint( 1, iseed )
      call fetchreal( 1, mass1 )
      call fetchreal( 1, mass2 )
      call fetchreal( 1, gamma )
      call fetchreal( 1, vfaint )
      call fetchreal( 1, fbinary )
      close(1)

      nhav = 0
      ncav = 0
      oldb = 0
      nstars = 0
      mrange = mass2**gamma - mass1**gamma

      ! Compute the number of crowding bins and number of pixels 
      ! in each dimension of each CMD.
      do icmd=1,NCMD
         nx(icmd) = int((xmax(icmd)-xmin(icmd))/dbinx)
         ny(icmd) = int((ycmax(icmd)-ymin(icmd))/dbiny)
         nxpix(icmd) = int((xmax(icmd)-xmin(icmd))/dpix) + 1
         nypix(icmd) = int((ypmax(icmd)-ymin(icmd))/dpix) + 1
      end do

      ! Initialize random number generator
      call sdprnd(iseed)
      
      call fetchchar8( 5, pre )
      call fetchint( 5, lockflag )
      call fetchint( 5, nscale )
      call fetchreal( 5, dmod )
      call fetchreal( 5, avfactor )
      call fetchreal( 5, gamma )
      call fetchreal( 5, fbinary )
      call fetchint( 5, sfrflag )

      call oldfile(10,agefile)
      if (lockflag.eq.1) then
         call oldfile(11,lockfile)
      endif

      ! Fill reddening vectors
      call read_av(hotfile,nhav,hav)
      call read_av(coldfile,ncav,cav)

      ! Tweak reddening distributions:
      do i=1,nhav
         hav(i) = hav(i)*avfactor
      end do

      do i=1,ncav
         cav(i) = cav(i)*avfactor
      end do

      do icmd=1,ncmd
         call strcat(pre,postfix(icmd),outfile(icmd))
         call unkfile(20+icmd, outfile(icmd))
      end do
      

      ! dtfile needed for calculating nn from SFRs
      if ( sfrflag.eq.1 ) then
         call oldfile(15,"input/dtime.dat")
      endif 

      ! Loop through isochrones
      do 30 i=1,NISOC
         if (lockflag.eq.1) then
            call readisodatline( 10, age, infile, libfile, vmsto )
            if ( age.lt.0.0 ) goto 99 ! end-of-file reached

c            do ipos=len(libfile),1,-1
c               if ( libfile(ipos:ipos).eq.'/' ) goto 10
c            end do
c 10         ipos = ipos + 1

            call readisolockline( 11, ib, stem1, stem )
            if ( ib.lt.0 ) goto 99 ! end-of-file reached

            if (ib.ne.oldb) then
               ! Done with current isochrone; initialize the next one.
               oldb = ib
               call readtestpopline( 5, a0, z, t, niso )
               if ( sfrflag.eq.1 ) read(15,*)   z2,t2,dtime

               if (stem1.ne.libfile) then
                  write(*,*) 'Error: mismatched libfiles'
                  write(*,*) libfile
                  write(*,*) stem1
                  pause
               endif

               if ( sfrflag.eq.1 .and. (z2.ne.z.or.t2.ne.t) ) then
                  write(*,*) "Error: mismatched dtime file"
                  write(*,*) z,t
                  write(*,*) z2,t2
               endif

               if (ib.gt.1) then
                  nstars = 0
                  nobs = 0
               endif
            endif
         else
            call readisodatline( 10, age, infile, libfile, vmsto )
            call readtestpopline( 5, a0, z, t, niso )

            if ( sfrflag.eq.1 ) read(15,*)   z2,t2,dtime

            libfile = '../' // libfile

            if (age.ne.t) then
               write(*,*) 'Error: mismatched libfiles'
               write(*,*) age, t
               pause
            endif

            if ( sfrflag.eq.1 .and. (z2.ne.z.or.t2.ne.t) ) then
               write(*,*) "Error: mismatched dtime file"
               write(*,*) z,t
               write(*,*) z2,t2
            endif
            
         endif

         if (a0.eq.0.0) goto 30
         if (iverb.gt.0) then
            write(*,*) 'Reading ',libfile
         endif

         libfile = "../" // libfile
         
         ! Read in the UBVI photometry, the masses and the relative 
         ! occupation probabilities for points in the current isochrone.
         call readphot(libfile, np)

         ! Recompute probabilities:
         flum(1) = mass(1)**gamma/gamma
         pcum(1) = 0.0

         do j=2,np
            flum(j) = mass(j)**gamma/gamma
            prob(j) = flum(j) - flum(j-1)
            pcum(j) = pcum(j-1) + prob(j)
         end do
         
         do j=1,np
            pcum(j) = pcum(j)/pcum(np)
         end do
         
         nstars = 0
         nobs = 0

         frac = (mass(np)**gamma - mass(imin)**gamma)/mrange

         ! Determine number of stars to create for this isochrone
         nn = nscale*a0*frac
         if (sfrflag.eq.1)  nn = nn*dtime/1.0e9
         if (lockflag.eq.1) nn = nn/niso

         write(*,*) imin,np,mass(imin),mass(np)
         do 20 j=1,nn

            ! choose stellar photometry randomly from the isochrone
            call pickphot(starmag, np)
            nstars = nstars + 1
            do imag=1,nmag
               starmag(imag) = starmag(imag) + dmod  !! tweak distance modulus
            end do

            ! Choose a random extinction value from an extinction vector.
            call redden(starmag, age)

            ! Combine magnitudes according to xeqn(icmd), yeqn(icmd) 
            ! to form CMD axes
            call axes(starmag,x,y)
            
            ! Record intrinsic photometry
            do icmd=1,ncmd
               x0(icmd) = x(icmd)
               y0(icmd) = y(icmd)
            end do

            ! Add photometric scatter
            if ( err_method .eq. 1 ) then
               call fakecrowd( x, y, idrop )
            else
               call crowd( x, y, idrop )
            endif

            ! Output star to files
            do icmd=1,ncmd
               if (idrop(icmd).eq.0) then
                  if (x(icmd).ge.xmin(icmd).and.
     x                x(icmd).le.xmax(icmd).and.
     x                y(icmd).ge.ymin(icmd).and.
     x                y(icmd).le.ypmax(icmd)) then
                     write(20+icmd,7) x(icmd),y(icmd),
     x                    x(icmd)-x0(icmd),y(icmd)-y0(icmd)
                     if (icmd.eq.2) nobs = nobs + 1
                  endif
               endif
            end do

 20      continue
         write(*,*) nstars,nobs,nobs*frac/nstars
 30   continue
      
 99   close(10)
      if (lockflag.eq.1) close(11)
      if (sfrflag.eq.1)  close(15)

      do icmd=1,ncmd
         close(20+icmd)
      end do

 7    format(4f7.3)

      stop
      end
