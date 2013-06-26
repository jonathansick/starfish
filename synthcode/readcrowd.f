      subroutine readcrowd()

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Read in the crowding bins' cumulative distributions from the *.cfn 
c     lookup table.  If the table does not exist, construct it with 
c     mkcfn.
c
      include 'synth.h'

      integer oflag,nbin
      integer ix,iy,icmd,ifilter
      real rowdata(NBINS)

      character*40 crowd1,crowd2
      common /c_cfile/ crowd1,crowd2

      integer nx(MCMD),ny(MCMD)
      real dbinx,dbiny,xmin(MCMD),ymin(MCMD)
      common /c_crowd/ dbinx,dbiny,nx,ny,xmin,ymin

      real emin,emax,dbin1
      common /c_chist/ emin,emax,dbin1

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      integer verb
      common /c_control/ verb

      real pdrop(MCMD,BINX,BINY)
      real cumerr(MCMD,BINX,BINY,NBINS*NBINS)
      common /c_cum/ pdrop,cumerr

      nbin = nint((emax-emin)/dbin1)

      open(unit=3,file=crowd2,status='old',iostat=oflag)

      if (oflag.eq.0) then
         if (verb.ge.2) then
            write(*,*) '   Reading existing crowding file...'
         endif

         do icmd=1,ncmd
            do ix=1,nx(icmd)
               do iy=1,ny(icmd)
                  k = 1
                  do jy=1,nbin
                     read(3,*) (rowdata(jx),jx=1,nbin)
                     do jx=1,nbin
                        cumerr(icmd,ix,iy,k) = rowdata(jx)
                        k = k + 1
                     end do
                  end do
                  read(3,*) pdrop(icmd,ix,iy), nobj
               end do
            end do
         end do
         close(3)

         return
      else
         close(3)

         if (verb.ge.2) write(*,*) '   Building crowding file...'
         call mkcfn()
      endif

      return
      end
