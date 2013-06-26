      subroutine grid(gtype,icmd,npix,nx,i,j,n)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     grid: Given pixel coords (i,j), returns corresponding box number
C     n, based either on a uniform 5x5 grid, or a custom grid, which 
C     must be coded by the user (see ubbgrid.f/bvvgrid.f/viigrid.f 
C     for an example).
C

      include 'sfh.h'

      integer gtype, icmd, i, j, n, npix, nx(MCMD)

c     Standard single regular grid
      if (gtype.eq.0) then
         n = 1 + int((i-1)/npix) + nx(icmd)*int((j-1)/npix)

c     custom grid
      else if (gtype.eq.1) then
         ! To use a custom grid, you need to define your own
         ! gridding code that determines to which CMD box 
         ! a star belongs, given its photometry
         ! (see the files ubbgrid.f, bvvgrid.f, and viigrid.f
         ! for an example).
         !
         ! Once this is done, replace the following code block 
         ! with call(s) to your custom grid subroutine(s).
         !
         ! You will also need to place the files defining these 
         ! subroutines in the Makefile.
         !
!         if (icmd.eq.1) then
!            call ubbgrid(i,j,n)
!         else if (icmd.eq.2) then
!            call bvvgrid(i,j,n)
!         else if (icmd.eq.3) then
!            call viigrid(i,j,n)
!         else
!            pause 'Invalid icmd value in grid.f'
!         endif

         ! Comment or remove the following block if you define 
         ! a custom grid:
         write(*,*) "Error: no custom grid defined."
         write(*,*) "Use normal grid instead."
         write(*,*) "(see grid.f for more information)"
         stop

      else
         write(*,*) 'Invalid grid type value in grid.f...check sfh.dat'
         stop
      endif

      return
      end
