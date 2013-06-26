      program mklib
c *** mklib:  Build a library of theoretical isochrones.  
c
c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     In the following discussion, the term "epoch" refers to the age 
c     and metallicity of a stellar population.
c
c     The input (Padua) isochrones give the photometry of a range of  
c     test masses at a given epoch.  However, the mass 
c     resolution is too coarse.  I need to simulate the photometry
c     resulting from a continuum of masses formed at a given epoch.  
c     The artifact of discrete mass values results in very lumpy 
c     Hess diagrams (CMDs populated by an isochrone, and according to some 
c     initial mass function).  In order to avoid this problem, I need 
c     to interpolate between neighboring isochrone points, such that 
c     the photometric distance between neighboring points will be 
c     much less than the characteristic photometric errors.  The Hess 
c     diagram is then indistinguishable from one drawn from a true continuum
c     of masses, if the points are error-blurred (as they will be).
c
c     For clarity, I will call the photometric locations of the original 
c     isochrone masses "isopoints", and the phot. locations of the 
c     interpolated isochrone points "interpoints".  The 
c     interpolation between neighboring isopoints is linear in photometric 
c     space.  However, the true isochrone is not a connected set 
c     of line segments, it is a smooth curve.  Thus, a simple linear 
c     interpolation of isopoint masses results in a slight mismatch between
c     the photometry of an interpoint and its mass (a better way to think 
c     of it is that the photometry of the interpoints is slightly wrong for
c     their mass).  The result of this mismatch is that the mass interval
c     per unit photometric interval is not a very smooth function of 
c     photometric position.  This is a problem, since the occupation 
c     probability at any point is determined by applying the IMF to the 
c     mass interval represented by the interpoint.  In other words, the 
c     Hess diagrams are still lumpy.  I kludge a solution to this problem 
c     by fitting a smooth curve through the \Delta_mass - phot. pos'n data,
c     and forcing the isochrones to have masses that conform to the curve 
c     fit.  I never change an interpoint's mass by more than a few
c     hundredths of a solar mass in making this correction.  
c
c     The interpolation is done only for masses below the main sequence
c     turn-off.  Beyond the turnoff, the mass resolution of the isopoints
c     becomes much finer.  Also, the \Delta_mass - phot. pos'n relation
c     should not necssarily be smooth in the post-main sequence phases.
c     I do not want to smooth over any real "lumps" in the Hess diagram.
c
      include "lib.h"

      character*120 linetc
      character*40 pfile,isofile,infile,libfile
      character*8 fz
      character*14 fstr(8)
      integer i,j,ii,jj,k,nintpts(NPTS),nmod(NPTS),iq(NPTS)
      integer dim,dim1,vindex
      real gamma,dmag,age,fiterr
      real msum,mass0,mcorrection
      real mag(MMAG,NPTS),mag2(MMAG,NPTS),dcos(MMAG)
      real mass(NPTS),mass2(NPTS),dm(NPTS),dm2(NPTS)
      real prob(NPTS),flum(NPTS),prob2(NPTS),flum2(NPTS)
      real s(NPTS),s2(NPTS),xd(NPTS),yd(NPTS),wd(NPTS),yfit(NPTS)
      real a(3*NPTS+3*MAXDIM+3)

      integer iverb
      common /c_verb/ iverb
      
      ! output format depends on number of magnitudes:
      fstr(1) = '(1f8.4,2f13.8)'
      fstr(2) = '(2f8.4,2f13.8)'
      fstr(3) = '(3f8.4,2f13.8)'
      fstr(4) = '(4f8.4,2f13.8)'
      fstr(5) = '(5f8.4,2f13.8)'
      fstr(6) = '(6f8.4,2f13.8)'
      fstr(7) = '(7f8.4,2f13.8)'
      fstr(8) = '(8f8.4,2f13.8)'

      ! Read in parameters from STDIN (unit=5)
      call fetchchar40( 5, isofile )
      call fetchreal( 5, vmax )
      call fetchreal( 5, dmag )
      call fetchreal( 5, dmod )
      call fetchreal( 5, gamma )
      call fetchint( 5, nmag )
      call fetchint( 5, vindex )
      call fetchint( 5, iverb )
      
      ! Read input file, output file, and MSTO magnitude from isofile
      call oldfile(1,isofile)
      do 9 jjj=1,NISO
         call readisodatline( 1, age, infile, libfile, vmsto )

         if ( infile.eq.'' ) goto 98  ! finished.
         if (iverb.gt.0) write(*,*) libfile

         vmsto = vmsto + dmod

         ii = 1
         call oldfile(2,infile)

         do i=1,NPTS
            read(2,*,end=11) mass(ii),(mag(imag,ii),imag=1,nmag)

            ! apply distance modulus
            do imag=1,nmag
               mag(imag,ii) = mag(imag,ii) + dmod
            end do

            ! keep only points brighter than vmax:
            if (mag(vindex,ii).le.vmax) ii = ii + 1

         end do
 11      close(2)
         nisopt = ii - 1


         ! s(i) is the cumulative phot. path length along the isopoints.
         ! s2(i) is the same, but for the interpoints.
         s(1) = 0
         s2(1) = 0

         ! i is the current interpoint index.  
         ! j is the current upper isopoint.
         i = 2
         j = 2

c     Begin locating interpoints (mag2) for the current isochrone.  
c     Interpoints will be placed along straight line segments every 
c     dmag between isopoints.  The first interpoint is the first isopoint.  
c     r is the photometric distance between lower and upper isopoints
c     r1 is the distance already covered by the previous interval.
         r = 0.0
         r1 = 0.0
         do imag=1,nmag
            mag2(imag,1) = mag(imag,1)
            r = r + (mag(imag,j) - mag(imag,j-1))**2
         end do
         mass2(1) = mass(1)

         r = sqrt(r)
         s(j) = s(j-1) + r

         ! Next, determine the direction cosines for each phot. dimension.
         do imag=1,nmag
            dcos(imag) = (mag(imag,j) - mag(imag,j-1))/r
         end do
         
         ! ninterp is the number of interpoints that can fit in 
         ! the current isopoint interval
 100     ninterp = int((r - r1)/dmag)

         ! Determine the interpoint photometry.  For now, don't worry 
         ! about assigning the mass to each interpoint.
         do k=1,ninterp
            do imag=1,nmag
               mag2(imag,i) = mag2(imag,i-1) + dcos(imag)*dmag
            end do
            s2(i) = s2(i-1) + dmag
            i = i + 1
         end do

c     The next interpoint needs to be dmag away from the current one, 
c     but that would overshoot the upper isopoint.  To deal with this, 
c     fit a line to the next interval, and find a point along *that* line
c     such that the total phot. distance from the previous interpoint
c     through the upper isopoint, to the next interpoint is dmag.
         j = j + 1
         if (j.gt.nisopt) goto 102

c     r2 is the photometric distance from the previous interpoint to the 
c     upper isopoint (it's less than dmag).  r1 is the phot. dist. that 
c     must be applied to the next isopoint interval (r1= dmag - r2).
c     r is the phot. dist. of the next isopoint interval.
         r2 = 0.0
         do imag=1,nmag
            r2 = r2 + (mag(imag,j-1)-mag2(imag,i-1))**2
         end do
         r2 = sqrt(r2)
         r1 = dmag - r2

 101     r = 0.0
         do imag=1,nmag
            r = r + (mag(imag,j)-mag(imag,j-1))**2
         end do
         r = sqrt(r)

         ! If r1 is greater than the current isopoint interval, 
         ! then try the next isopoint.  Repeat as necessary.
         if (r1.gt.r) then
            s(j) = s(j-1) + r
            r1 = r1 - r
            j = j + 1
            if (j.gt.nisopt) goto 102
            goto 101
         endif
         
         s(j) = s(j-1) + r
         s2(i) = s2(i-1) + dmag

         do imag=1,nmag
            dcos(imag) = (mag(imag,j) - mag(imag,j-1))/r
            mag2(imag,i) = mag(imag,j-1) + dcos(imag)*r1
         end do

         i = i + 1
         goto 100

c     Now I have the photometric resolution that I need.  Next, I need to
c     determine the masses of the interpoints.  Get a first guess 
c     at the interpoint masses by a simple linear interpolation from the 
c     surrounding isopoint masses.
 102     call newfile(3,libfile)
         nn = i - 1
         jj = 1
         prob2(1) = 0.0

         if (iverb.gt.3) write(20,*) age,1,s(1),mag(vindex,1),mass(1)
         do 110 ii=2,nisopt
            if (iverb.gt.3) write(20,*) age,ii,s(ii),
     x                                  mag(vindex,ii),mass(ii)

 103        if (s2(jj).gt.s(ii)) goto 110

            if (jj.gt.nn) goto 120

            rr = s(ii) - s(ii-1)
            xa = (s2(jj) - s(ii-1))/rr
            xb = (s(ii) - s2(jj))/rr
            mass2(jj) = mass(ii)*xa + mass(ii-1)*xb

            if (iverb.gt.3) 
     x           write(21,*) age,jj,s2(jj),mag2(vindex,jj),mass2(jj)
            if (jj.gt.1) dm2(jj) = mass2(jj) - mass2(jj-1)

            flum2(jj) = (mass2(jj)**gamma)/gamma
            if (jj.gt.1) prob2(jj) = (flum2(jj) - flum2(jj-1))

            jj = jj + 1
            goto 103
 110     continue

 120     dm2(1) = dm2(2)
         ninterp = jj - 1

c     These interpolated masses are not smoothly distributed photometrically,
c     due to the slight error introduced by connecting isopoints
c     with straight line segments.  This error is amplified when computing 
c     the occupation probability, which involves dividing by the mass 
c     interval between adjacent interpoints; the interval is the difference 
c     between the masses of adjacent interpoints, two very similar numbers.
c     The solution is to force the delta masses to be a smooth function
c     with photometric position by fitting a low-order polynomial through 
c     the first-guess interpolated masses.

c     First, define a new origin for the phot. distances vs. 
c     delta-masses plot to simplify polynomial fitting.
         dm0 = dm2(1)       !! zeropoint of delta masses
         x0  = mag2(vindex,1)    !! zeropoint of phot. position
         if (age.lt.8.5) x0 = 23.0
         do i=1,ninterp
            if(mag2(vindex,i).lt.vmsto) goto 150
            xd(i) = mag2(vindex,i) - x0
            yd(i) = dm2(i) - dm0
            wd(i) = 1.0
         end do

 150     nd = i - 1
         mass0 = mass2(nd)
         chimin = 1.0e6

c     fitpoly returns the params and chisq of the best-fit polynomial 
c     of order dim (using linear SVD, se NR).
c     Loop through dim=1 to 20 and take the lowest chisq as the best
c     fitting polynomial
c         call fitpoly(xd,yd,sd,nd,a,ma,u,v,ww,chisq)
         fiterr=-1.0
         call polfit(nd,xd,yd,wd,MAXDIM,dim,fiterr,yfit,ierr,a)
         if (iverb.gt.3) write(*,*) age, ": ", dim

c     Need to correct the smoothed masses to ensure that the mass at 
c     vmsto (where the smoothing ends) is equal to the original isochrone 
c     mass at this point. 
         msum = mass2(1)
         do i=2,nd
            dm2(i) = yfit(i) + dm0
            msum = msum + dm2(i)
         end do
            
         mcorrection = (mass0 - msum)/(nd)

         ! Write delta-mass, polynomial dat to a file (for inspection)
         do i=2,nd
            xd(i) = xd(i) + x0
            yd(i) = yd(i) + dm0
            dm2(i) = dm2(i) + mcorrection
            if (iverb.gt.3) write(13,*) age,xd(i),yd(i),dm2(i)
         end do

         ! Output interpolated isochrone
         write(3,fstr(nmag)) (mag2(imag,1),imag=1,nmag),
     x        prob2(1),mass2(1)

         do i=2,ninterp
            mass2(i) = mass2(i-1) + dm2(i)
            flum2(i) = (mass2(i)**gamma)/gamma

            if (mass2(i).eq.mass2(i-1)) then
               prob2(i) = 0.0
            else 
               prob2(i) = (flum2(i) - flum2(i-1))
            endif

            write(3,fstr(nmag)) (mag2(imag,i),imag=1,nmag),
     x           prob2(i),mass2(i)
         end do
            
         close(3)

 9    continue

      close(1)
      close(2)

 4    format(f7.2)
 6    format(a16,i8,f7.3)
 7    format(f12.8)

 98   stop
      end


