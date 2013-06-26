      subroutine bincombine(icmd,nstars,ndrop,dx)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     This file is optional, and inactive by default.  If there are CMD 
c     regions that are sparsely populated with artificial stars (for 
c     examply, you used isochrones to select intrinsic photometry), then
c     many boxes will have too few artificial stars to reliably determine
c     error statistics.  bincombine will combine adjacent boxes to 
c     increase the number of stars.  You lose photometric resolution in
c     doing so; however, since the sparse boxes are typically at bright
c     magnitudes, photometric resolution isn't necessary.
c
c     If you use bincombine, this file will need to be customized to fit
c     your particular CMD geometry.  I include a version used for our
c     MCPS data as a reference.  Once you have done this, uncomment 
c     the 'call bincombine' line in mkcfn.f, and run make.

      include 'synth.h'

      integer i,icmd,ix,iy,ix2,iy2,num,istar
      integer ndrop(BINX,BINY),nstars(BINX,BINY)
      real dx(2,BINX,BINY,MAXPOP),dx2(2,BINX,BINY,MAXPOP)

      integer nx(MCMD),ny(MCMD)
      real dbinx,dbiny,xmin(MCMD),ymin(MCMD)
      common /c_crowd/ dbinx,dbiny,nx,ny,xmin,ymin

c     Copy dx to dx2
      do i=1,2
         do ix=1,nx(icmd)
            do iy=1,ny(icmd)
               do istar=1,nstars(ix,iy)
                  dx2(i,ix,iy,istar) = dx(i,ix,iy,istar)
               end do
            end do
         end do
      end do

c     Combine the top five rows of each CMD into one large bin.  
      num = 0
      nd = 0
      do ix=1,nx(icmd)
         do iy=1,8              !! top 8 rows
            do inum=1,nstars(ix,iy) !! loop through stars in each bin
               do ix2=1,nx(icmd)
                  do iy2=1,8    !! assign the current star to each bin
                     dx(1,ix2,iy2,num+inum) = dx2(1,ix,iy,inum)
                     dx(2,ix2,iy2,num+inum) = dx2(2,ix,iy,inum)
                  end do
               end do
            end do
            num = num + nstars(ix,iy)
            nd = nd + ndrop(ix,iy)
         end do
      end do
      
      do ix=1,nx(icmd)
         do iy=1,8
            nstars(ix,iy) = num
            ndrop(ix,iy) = nd
         end do
      end do
      

      do iy=9,ny(icmd)
         num = 0
         nd = 0
         do ix=1,nx(icmd)
            do inum=1,nstars(ix,iy)
               do ix2=1,nx(icmd)
                  dx(1,ix2,iy,num+inum) = dx2(1,ix,iy,inum)
                  dx(2,ix2,iy,num+inum) = dx2(2,ix,iy,inum)
               end do
            end do
            num = num + nstars(ix,iy)
            nd = nd + ndrop(ix,iy)
         end do

         do ix=1,nx(icmd)
            nstars(ix,iy) = num
            ndrop(ix,iy) = nd
         end do
      end do
            
      return
      
      !!!Ignore the rest of the file (old junk)


      if (icmd.eq.1) then          ! U-B, B:     
         num = 0
         do iy=6,9              ! combine next 4 rows into one bin
            do ix=1,nx(icmd)
               do inum=1,nstars(ix,iy)
                  do iy2=6,9
                     do ix2=1,nx(icmd)
                        dx(1,ix2,iy2,num+inum) = dx2(1,ix,iy,inum)
                        dx(2,ix2,iy2,num+inum) = dx2(2,ix,iy,inum)
                     end do
                  end do
               end do
               num = num + nstars(ix,iy)
            end do
         end do
         do iy=6,9
            do ix=1,nx(icmd)
               nstars(ix,iy) = num
            end do
         end do
         
         do irow=10,12,2        ! combine next 2 pairs of rows into 2 bins each
            num = 0
            do iy=irow,irow+1
               do ix=1,3
                  do inum=1,nstars(ix,iy)
                     do iy2=irow,irow+1
                        do ix2=1,3
                           dx(1,ix2,iy2,num+inum) = dx2(1,ix,iy,inum)
                           dx(2,ix2,iy2,num+inum) = dx2(2,ix,iy,inum)
                        end do
                     end do
                  end do
                  num = num + nstars(ix,iy)
               end do
            end do
            do iy=irow,irow+1
               do ix=1,3
                  nstars(ix,iy) = num
               end do
            end do
            
            num = 0
            do iy=irow,irow+1
               do ix=4,nx(icmd)
                  do inum=1,nstars(ix,iy)
                     do iy2=irow,irow+1
                        do ix2=4,nx(icmd)
                           dx(1,ix2,iy2,num+inum) = dx2(1,ix,iy,inum)
                           dx(2,ix2,iy2,num+inum) = dx2(2,ix,iy,inum)
                        end do
                     end do
                  end do
                  num = num + nstars(ix,iy)
               end do
            end do
            do iy=irow,irow+1
               do ix=4,nx(icmd)
                  nstars(ix,iy) = num
               end do
            end do
         end do
         
         do iy=14,ny(icmd)         ! 2 boxes in each of the remaining rows: 
            num = 0
            do ix=1,3
               do inum=1,nstars(ix,iy)
                  do ix2=1,3
                     dx(1,ix2,iy,num+inum) = dx2(1,ix,iy,inum)
                     dx(2,ix2,iy,num+inum) = dx2(2,ix,iy,inum)
                  end do
               end do
               num = num + nstars(ix,iy)
            end do
            do ix=1,3
               nstars(ix,iy) = num
            end do
            
            num = 0
            do ix=4,nx(icmd)
               do inum=1,nstars(ix,iy)
                  do ix2=4,nx(icmd)
                     dx(1,ix2,iy,num+inum) = dx2(1,ix,iy,inum)
                     dx(2,ix2,iy,num+inum) = dx2(2,ix,iy,inum)
                  end do
               end do
               num = num + nstars(ix,iy)
            end do
            do ix=4,nx(icmd)
               nstars(ix,iy) = num
            end do
         end do
      else                      ! B-V, V; V-I, I:       
         do irow=6,8,2     ! combine next 2 pairs of rows into two bins each
            num = 0
            do iy=irow,irow+1    
               do ix=1,2
                  do inum=1,nstars(ix,iy)
                     do iy2=irow,irow+1
                        do ix2=1,2
                           dx(1,ix2,iy2,num+inum) = dx2(1,ix,iy,inum)
                           dx(2,ix2,iy2,num+inum) = dx2(2,ix,iy,inum)
                        end do
                     end do
                  end do
                  num = num + nstars(ix,iy)
               end do
            end do
            do iy=irow,irow+1
               do ix=1,2
                  nstars(ix,iy) = num
               end do
            end do
            
            num = 0
            do iy=irow,irow+1
               do ix=3,nx(icmd)
                  do inum=1,nstars(ix,iy)
                     do iy2=irow,irow+1
                        do ix2=3,nx(icmd)
                           dx(1,ix2,iy2,num+inum) = dx2(1,ix,iy,inum)
                           dx(2,ix2,iy2,num+inum) = dx2(2,ix,iy,inum)
                        end do
                     end do
                  end do
                  num = num + nstars(ix,iy)
               end do
            end do
            do iy=irow,irow+1
               do ix=3,nx(icmd)
                  nstars(ix,iy) = num
               end do
            end do
         end do
         
         do iy=10,15            ! 2 boxes in each row until V=20: 
            num = 0
            do ix=1,2
               do inum=1,nstars(ix,iy)
                  do ix2=1,2
                     dx(1,ix2,iy,num+inum) = dx2(1,ix,iy,inum)
                     dx(2,ix2,iy,num+inum) = dx2(2,ix,iy,inum)
                  end do
               end do
               num = num + nstars(ix,iy)
            end do
            do ix=1,2
               nstars(ix,iy) = num
            end do
            
            num = 0
            do ix=3,nx(icmd)
               do inum=1,nstars(ix,iy)
                  do ix2=3,nx(icmd)
                     dx(1,ix2,iy,num+inum) = dx2(1,ix,iy,inum)
                     dx(2,ix2,iy,num+inum) = dx2(2,ix,iy,inum)
                  end do
               end do
               num = num + nstars(ix,iy)
            end do
            do ix=3,nx(icmd)
               nstars(ix,iy) = num
            end do
         end do
         
c     Original resolution for the remaining rows.
         do iy=16,ny(icmd)
            do ix=1,nx(icmd)
               do inum=1,nstars(ix,iy)
                  dx(1,ix,iy,inum) = dx2(1,ix,iy,inum)
                  dx(2,ix,iy,inum) = dx2(2,ix,iy,inum)
               end do
            end do
         end do
      endif 

      return
      end
