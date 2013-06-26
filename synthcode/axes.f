      subroutine axes(starmag,x,y)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Constructs CMD axes according to the equation strings 
c     (xeqn(ncmd), yeqn(ncmd)) found in synth.dat.
c
      include 'synth.h'

      character*8 str
      integer icmd,xy,col(MCMD,2,2),iflag
      real x(MCMD),y(MCMD),sum
      real starmag(MMAG)
      logical number,add(MCMD,2)
c     number = .true.  :  expect a number
c     number = .false. :  expect +/-

      save add, col, iflag

      character*8 xeqn(MCMD),yeqn(MCMD)
      common /c_eqn/ xeqn, yeqn

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

c     On first call, parse the eqn strings and store numbers, operators
c     in ixcol/iycol, xadd/yadd
      if (iflag.ne.1) then
         iflag = 1
         
         do icmd=1,ncmd
            do xy=1,2
               if (xy.eq.1) str = xeqn(icmd)
               if (xy.eq.2) str = yeqn(icmd)
               jcol = 1
               col(icmd,xy,2) = 0  ! default to no second argument
               number = .true.

               do i=1,len(str)
                  if (str(i:i).ne.' ') then
                     value = ichar(str(i:i)) - 48
                     if (number) then
                        if (value.lt.0.or.value.gt.9) then ! not a digit
                           write(*,*) str
                           pause 'parse error: digit expected'
                        else
                           col(icmd,xy,jcol) = value
                           jcol = jcol + 1
                           if (jcol.eq.3) goto 10 ! done
                        endif
                     else
                        add(icmd,xy) = .true.
                        if (value.eq.-3) add(icmd,xy) = .false. ! read '-'
                     endif
                     number = .not. number
                  endif
               end do
 10            continue
            end do
         end do
      endif

      do icmd=1,ncmd
         do xy=1,2
            if (col(icmd,xy,1).eq.0)
     x           pause 'Error parsing magnitudes in axes.f '

            sum = starmag(col(icmd,xy,1))
            if (col(icmd,xy,2).ne.0) then
               if (add(icmd,xy)) then
                  sum = sum + starmag(col(icmd,xy,2))
               else
                  sum = sum - starmag(col(icmd,xy,2))
               endif
            endif

            if (xy.eq.1) x(icmd) = sum
            if (xy.eq.2) y(icmd) = sum
         end do
      end do

      return
      end
