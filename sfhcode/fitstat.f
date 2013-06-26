      subroutine fitstat(x,y,iwrt)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     fitstat: Determine value of fitting statistic for current model.  
C     The fitstat is either chi-squared, Lorentian stat, or Poisson stat. 
C
C     The first time fitstat is called, it populates the data and synthetic 
C     CMD grid boxes. 
C
      include 'sfh.h'

      double precision x(MP),y,c,c2,pix,pjunk
      integer i,j,ii,jj,iflag,numx(MCMD),numy(MCMD),ix,iy
      integer ibox,nfull(MCMD),iwrt,ngridx(MCMD)
      integer ifull(MCMD,MBOX)
      character*40 pixbase
      character*40 pixfile

      save iflag,ifull,nfull, numx, numy, ngridx

      real dat(MCMD,2,NDATA),chem(MP),age(MP)
      integer fstat,gtype
      integer nstars(MCMD),pnum(MP),mask(MCMD,MBOX)
      character*40 cmdfile,chifile
      common /cstat/ dat,cmdfile,chifile,nstars,mask,
     x               chem,age,pnum,fstat,gtype

      double precision psynth(MCMD,MP,MBOX)
      double precision pmod(MCMD,MBOX),pdat(MCMD,MBOX)
      common /cfit/ psynth,pmod,pdat

      integer niso,nheld,ihold(MP)
      double precision pheld(MP)
      common /hold/ pheld,ihold,niso,nheld

      integer np
      common /cnp/ np

      integer iverb
      common /verb/ iverb

      character*8 suffix(MCMD)
      integer npix, ncmd, nbox(MCMD)
      real xlo(MCMD), xhi(MCMD), ylo(MCMD), yhi(MCMD), dpix
      common /cmds/ suffix, xlo, xhi, ylo, yhi, nbox, dpix, npix, ncmd

c     ***                                               ***
c     *** Do the following on the first call to fitstat ***
c     ***                                               ***
      if (iflag.ne.1) then
         iflag = 1              
         
         do icmd=1,ncmd
            numx(icmd) = int((xhi(icmd) - xlo(icmd))/dpix) + 1
            numy(icmd) = int((yhi(icmd) - ylo(icmd))/dpix) + 1
            ngridx(icmd) = numx(icmd)/npix 
         end do

         do i=1,np
            pnum(i)  = 0
         end do

         do icmd=1,ncmd
            nfull(icmd) = nbox(icmd)

            do i=1,nbox(icmd)
               pdat(icmd,i)  = 0.0D0
               ifull(icmd,i) = 0
            end do
         end do


c     *** Populate the CMD grids with model stars ***
         call oldfile(7,cmdfile)
         do i=1,niso            ! loop over all isochrones
            call readcmdline( 7, chem(i), age(i), pixbase )
            if ( chem(i).lt.0.0 ) goto 10

            if (iverb.ge.4) write(*,*) pixbase

            pnum(i) = pnum(i) + 1
            do icmd=1,ncmd
               call strcat(pixbase, suffix(icmd), pixfile)
               call oldfile(8,pixfile)

               do jj=1,numy(icmd)             ! loop over 
                  do ii=1,numx(icmd)    ! pixels
                     read(8,8) pix

                     call grid(gtype,icmd,npix,ngridx,ii,jj,ibox) 
                     if (mask(icmd,ibox).eq.1) ibox = 0
                     if (ibox.gt.0) then
                        psynth(icmd,i,ibox) = 
     x                       psynth(icmd,i,ibox) + pix
                     endif
                  end do                                         
               end do
               close(8)
            end do
         end do
 10      close(7)

c        Renormalize the basis functions.
         do i=1,np           !loop over indep. isoc.
            do icmd=1,ncmd
               sum = 0.0D0
               do j=1,nbox(icmd)
                  sum = sum + psynth(icmd,i,j)
               end do
               
               do j=1,nbox(icmd)
                  if (sum.gt.0.0) then
                     pjunk = psynth(icmd,i,j)/sum
                  else
                     pjunk = 0.0
                  endif

                  if (pjunk.gt.1.0e-9) ifull(icmd,j)=1 
               end do                                  
            end do                                
         end do      


c     *** Populate the CMD grids with data stars ***
         do icmd=1,ncmd
            do 21 i=1,nstars(icmd)    ! loop over data stars
               ix = int((dat(icmd,1,i)-xlo(icmd))/dpix) + 1
               iy = numy(icmd) - int((dat(icmd,2,i)-ylo(icmd))/dpix)
               if (ix.lt.1.or.ix.gt.numx(icmd) .or. 
     x             iy.lt.1.or.iy.gt.numy(icmd)) goto 21 
            
               call grid(gtype,icmd,npix,ngridx,ix,iy,ibox)
               if (mask(icmd,ibox).eq.1) ibox = 0
               if (ibox.gt.0) then
                  pdat(icmd,ibox) = pdat(icmd,ibox) + 1.0D0
                  ifull(icmd,ibox) = 2 ! ifull: 0=empty, 1=model, 2=data
               endif
 21         continue

c     decrement nfull when a box is empty in both data+model
            do j=1,nbox(icmd)
               if (ifull(icmd,j).eq.0) nfull(icmd) = nfull(icmd) - 1 
            end do              
         end do
      endif                    
c     *** End first-run if block ***

      c2 = 0.0D0
      if (iwrt.eq.1) call unkfile(16,chifile)

      do 22 icmd=1,ncmd
         do 23 ibox=1,nbox(icmd)

c           Reset pmod
            pmod(icmd,ibox) = 0.0D0

            if (ifull(icmd,ibox).ne.0) then

c              pmod is the predicted #stars in the current gridbox, 
c              according to the current model SFH.  Each amplitude is
c              either free (one of x(np)), or held fixed (one of 
c              pheld(nheld).
               jfree = 1
               do j=1,niso
                  if ( ihold(j) .eq. 0 ) then  ! free amplitude
                     pmod(icmd,ibox) = pmod(icmd,ibox) + 
     x                    x(jfree)*psynth(icmd,jfree,ibox)
                     jfree = jfree + 1
                  else  ! fixed amplitude
                     pmod(icmd,ibox) = pmod(icmd,ibox) + 
     x                    pheld(ihold(j))*psynth(icmd,ihold(j),ibox)
                  endif
               end do

c           Determine contribution of current gridbox to the fit statistic:
               if ( fstat.eq.0 ) then  !chi-squared
                  if (pdat(icmd,ibox).eq.0.0D0) then
c                    In this case, we ad-hoc insert pdat=1.0, 
c                    which is assumed to be "close to zero"
                     c = (pdat(icmd,ibox) - pmod(icmd,ibox))**2
                  else 
                     c = (pdat(icmd,ibox) - pmod(icmd,ibox))**2/
     x                    pdat(icmd,ibox)
                  endif

               else if ( fstat.eq.1 ) then  !Lorentzian 
c                 Equivalent to 1 + 0.5*chi^2
                  if (pdat(icmd,ibox).eq.0.0D0) then
                     c = (pdat(icmd,ibox) - pmod(icmd,ibox))**2
                  else 
                     c = (pdat(icmd,ibox) - pmod(icmd,ibox))**2/
     x                    pdat(icmd,ibox)
                  endif
                  c = 1.0 + 0.5*c

               else  !Poisson
                  if (pdat(icmd,ibox).eq.0.0D0 .and. 
     x                pmod(icmd,ibox).eq.0.0D0 ) then
                     c = 0.0D0
                  else if (pdat(icmd,ibox).eq.0.0D0 .and. 
     x                pmod(icmd,ibox).ne.0.0D0 ) then
                     c = 2.0 * pmod(icmd,ibox)
                  else if (pmod(icmd,ibox).eq.0.0D0 ) then
c                    In this case, we ad-hoc insert pmod=0.01, 
c                    which is assumed to be "close to zero"
                     c = 2.0 * ( 0.01 - pdat(icmd,ibox) + 
     x                    pdat(icmd,ibox) * log(pdat(icmd,ibox)/0.01 ))
                  else
                     c = 2.0 * ( pmod(icmd,ibox) - pdat(icmd,ibox) +
     x                    pdat(icmd,ibox) * 
     x                    log( pdat(icmd,ibox)/pmod(icmd,ibox) ))
                  endif
               endif
            else
               pmod(icmd,ibox) = 0.0D0
               c = 0.0D0
            endif
            
            c2 = c2 + c
            
            if (iwrt.eq.1) 
     x           write(16,*) icmd,ibox,pmod(icmd,ibox),pdat(icmd,ibox),c
 23      continue
 22   continue

      if (iwrt.eq.1) close(16)

      y = c2                  ! The result: fitting stat for current model.

 6    format(2a8,a24)
 8    format(g14.7)
 9    format(i4,31g14.7,g16.9)

      return
      end

