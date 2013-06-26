      subroutine scatter( x, y, dropflag )

      include 'synth.h'

      integer iflag, n(MCMD, BINX, BINY)
      integer dropflag(MCMD)
      real w(MMAG), dw(MMAG)
      real x(MCMD), y(MCMD), dx(MCMD), dy(MCMD)
      real dx0(MCMD, BINX, BINY, MAXPOP)
      real dy0(MCMD, BINX, BINY, MAXPOP)
      double precision p

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      character*40 crowd1,crowd2
      common /c_cfile/ crowd1,crowd2
      
      integer nx(MCMD),ny(MCMD)
      real dbinx,dbiny,xmin(MCMD),ymin(MCMD)
      common /c_crowd/ dbinx,dbiny,nx,ny,xmin,ymin

      save iflag, n, dx0, dy0

C *** Run this block on first call only:
      if ( iflag.eq.0 ) then
         iflag = 1

         do icmd=1,ncmd
            do ix=1,nx(icmd)
               do iy=1,ny(icmd)
                  n(icmd,ix,iy) = 0
               end do
            end do
         end do

         !Read crowding table
         call oldfile( 1, crowd1 )
 10      read(1,*,end=20) xcoo,ycoo,(w(imag),dw(imag),imag=1,nmag)

         ! If the star was detected in band j, convert w(j)
         ! to the intrinsic mag by subtracting dw(j)
         do j=1,nmag
            if (dw(j).lt.9.9) w(j) = w(j) - dw(j)
         end do
         
         call axes(  w,  x,  y )
         call delta_axes( dw, dx, dy )

         ! Determine the crowd bin to which the star belongs in each CMD.
         do 12 icmd=1,ncmd
            ix = int((x(icmd) - xmin(icmd))/dbinx) + 1
            iy = int((y(icmd) - ymin(icmd))/dbiny) + 1
            if (ix.lt.1 .or. ix.gt.nx(icmd)) goto 12
            if (iy.lt.1 .or. iy.gt.ny(icmd)) goto 12
            
            if ( n(icmd,ix,iy) .ge. MAXPOP ) then
c               write(*,*) "ERROR: box population exceeds MAXPOP!"
c               write(*,*) icmd, ix, iy
            else
               n(icmd,ix,iy) = n(icmd,ix,iy) + 1
               dx0(icmd, ix, iy, n(icmd,ix,iy) ) = dx(icmd)
               dy0(icmd, ix, iy, n(icmd,ix,iy) ) = dy(icmd)
            endif
 12      continue

         goto 10

 20      close(1)
      endif
C *** End First-run block

      ! Determine the crowd bin to which the target model star belongs
      ! Drop stars which lie outside the crowding region
      do icmd=1,ncmd
         dropflag(icmd) = 0

         ix = int((x(icmd) - xmin(icmd))/dbinx) + 1
         iy = int((y(icmd) - ymin(icmd))/dbiny) + 1

         if (ix.lt.1.or.ix.gt.nx(icmd) .or.
     x       iy.lt.1.or.iy.gt.ny(icmd)) then
            dropflag(icmd) = 1
            dx1 = 0.0
            dy1 = 0.0
            goto 200 ! skip phot err stuff
         endif

         ! Draw randomly from among the artificial stars 
         ! in this crowding bin
         call random( p, 0.0D0, dble(n(icmd,ix,iy)) )
         ii = int(p) + 1

         dx1 = dx0(icmd,ix,iy,ii)
         dy1 = dy0(icmd,ix,iy,ii)
         
         if ( dx1.gt.90.0 .or. dy1.gt.90.0 ) then
            dropflag(icmd) = 1
            dx1 = 0.0
            dy1 = 0.0
         endif

 200     x(icmd) = x(icmd) + dx1
         y(icmd) = y(icmd) + dy1
      end do

      return
      end
