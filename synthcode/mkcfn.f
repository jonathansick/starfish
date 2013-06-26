      subroutine mkcfn()

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Construct the *.cfn file, a lookup table of cumulative distributions 
c     for all crowding bins in the CMDs.
c
      include 'synth.h'

      integer i,ii,j,num
      integer nstars(BINX,BINY),ndrop(BINX,BINY)
      real xcoo,ycoo,w(MMAG,MAXSTR),dw(MMAG,MAXSTR)
      real mag(MMAG),dmag(MMAG)
      real x(MCMD),y(MCMD),xoff(MCMD),yoff(MCMD)
      real dx(2,BINX,BINY,MAXPOP),dx2(2,BINX,BINY,MAXPOP)
      real dat(MAXPOP),xc(NBINS),yc(NBINS)

      integer nbin
      real fcum
      integer epix(BINX,BINY,NBINS,NBINS)
      
      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      character*40 crowd1,crowd2
      common /c_cfile/ crowd1,crowd2

      integer nx(MCMD),ny(MCMD)
      real dbinx,dbiny,xmin(MCMD),ymin(MCMD)
      common /c_crowd/ dbinx,dbiny,nx,ny,xmin,ymin

      real emin,emax,dbin1
      common /c_chist/ emin,emax,dbin1

      real pdrop(MCMD,BINX,BINY)
      real cumerr(MCMD,BINX,BINY,NBINS*NBINS)
      common /c_cum/ pdrop,cumerr

      nbin = nint((emax-emin)/dbin1)
      i = 1

c     Read input file.
      call oldfile(1,crowd1)
 10   read(1,*,end=11) xcoo,ycoo,(w(imag,i),dw(imag,i),imag=1,nmag)

c     If the star was detected in band j, convert w(j) to 
c     the intrinsic mag by subtracting dw(j)
      do j=1,nmag
         if (dw(j,i).lt.9.9) w(j,i) = w(j,i) - dw(j,i)
      end do

      i = i + 1
      goto 10
      
 11   close(1)
      ii = i - 1

      call unkfile(3,crowd2)

      do 20 icmd=1,ncmd
         do ix=1,nx(icmd)
            do iy=1,ny(icmd)
               nstars(ix,iy) = 0
               ndrop(ix,iy) = 0
               do jx=1,nbin
                  do jy=1,nbin
                     epix(ix,iy,jx,jy) = 0
                  end do
               end do
            end do
         end do

c     Loop over all stars in crowding table.
         do i = 1,ii
            do imag=1,nmag
               mag(imag) = w(imag,i)
               dmag(imag) = dw(imag,i)
            end do

c     convert mag and dmag into color/mag CMD positions (x/y)
            call axes(mag,x,y)
            call delta_axes(dmag,xoff,yoff)
            
c     Determine to which crowding bin the current star belongs.
            ix = int((x(icmd) - xmin(icmd))/dbinx) + 1
            if (ix.lt.1) ix = 1
            if (ix.gt.nx(icmd)) ix = nx(icmd)
            
            iy = int((y(icmd) - ymin(icmd))/dbiny) + 1
            if (iy.lt.1) iy = 1
            if (iy.gt.ny(icmd)) iy = ny(icmd)
            
c     Determine the indices of the "error-image" pixel
c     to which the current star belongs
            jx = int((xoff(icmd) - emin)/dbin1) + 1
            jy = int((yoff(icmd) - emin)/dbin1) + 1

c     Add the star to the list of stars in its crowding bin.
c     If the star is a non-detection, increment ndrop(ix,iy)
c     Otherwise, increment nstars(ix,iy) and increment
c     the pixel value of the "error-image" for this crowd bin
            if (nstars(ix,iy).lt.MAXPOP) then
               if (xoff(icmd).gt.emax.or.xoff(icmd).lt.emin
     x              .or.yoff(icmd).gt.emax.or.yoff(icmd).lt.emin) then

                  ndrop(ix,iy) = ndrop(ix,iy) + 1
               else
                  nstars(ix,iy) = nstars(ix,iy) + 1
                  epix(ix,iy,jx,jy) = epix(ix,iy,jx,jy) + 1
               endif
            endif
         end do


ccc      By default, we do not call bincombine.  
ccc      See bincombine.f for details
ccc         call bincombine(icmd,nstars,ndrop,dx)

c     Loop through each crowding bin.  Compute the dropout fraction
c     and the cumulative 2-D delta-mag distribution
         do ix=1,nx(icmd)
            do iy=1,ny(icmd)
               ntot = nstars(ix,iy) + ndrop(ix,iy)

c     Set the default dropout fraction to 1.0/0.0 for faint/bright bins.
               pdrop(icmd,ix,iy) = 0.0
               if (iy.ge.int(0.5*ny(icmd))) then
                  pdrop(icmd,ix,iy) = 1.0
               endif
               
               if (ntot.gt.MINPOP)                
     x              pdrop(icmd,ix,iy)=real(ndrop(ix,iy))/real(ntot)
               
               ! The 2-D cumulative delta-mag distribution
               ! if the bin is unpopulated, use a null cumerr fcn
               if ( nstars(ix,iy) .eq. 0 ) then
                  do k=1,nbin*nbin
                     if ( k .lt. 0.5*nbin*(nbin + 1) ) then
                        cumerr(icmd,ix,iy,k) = 0.0
                     else
                        cumerr(icmd,ix,iy,k) = 1.0
                     endif
                  end do
               else                     
                  fcum = 0.0
                  k = 1
                  do jy=1,nbin
                     do jx=1,nbin
                        fcum = fcum + 
     x                       real(epix(ix,iy,jx,jy))/real(nstars(ix,iy))
                        cumerr(icmd,ix,iy,k) = fcum
                        k = k + 1
                     end do
                     
                     ! Make sure final pixel is equal to 1.0
                     if ( jy.eq.nbin .and. 
     x                    cumerr(icmd,ix,iy,nbin*nbin).ne.1.0 )
     x                    cumerr(icmd,ix,iy,nbin*nbin) = 1.0
                  end do
               endif

               ! Write to the output file
               do jy=1,nbin
                  write(3,3) (cumerr(icmd,ix,iy,
     x                 (jy-1)*nbin+jx),jx=1,nbin)
               end do

               write(3,4) pdrop(icmd,ix,iy),nstars(ix,iy)
            end do
         end do
 20   continue

      close(3)

 3    format(120f7.4)
 4    format(f11.8,2x,i6)

      return
      end
                  
