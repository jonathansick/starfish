      program mkgrid
      integer ncmd, nx(3), ny
      ncmd = 3
      nx(1) = 70
      nx(2) = 50
      nx(3) = 60
      ny = 200
      nbin = 5
      
      do icmd=1,ncmd
         open(unit=11,file="grid"//char(icmd+48)//".index",type="unknown")

         do ix=1,nx(icmd)
            do iy=1,ny
               nbox = 1 + int((ix-1)/nbin) + 
     x              int((iy-1)/nbin)*nx(icmd)/nbin
               
               write(11,5) ix,iy,nbox
            end do
         end do
         close(11)
      end do


 5    format(3i5)

      stop
      end
