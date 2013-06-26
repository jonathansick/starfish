      subroutine logp( logfile, pmin, ymin, lamda )

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Log the current parameter-space location to the logfile

      include 'sfh.h'

      character*40 logfile
      double precision pmin(MP),ymin
      real lamda

      integer np
      common /cnp/ np

      integer iverb
      common /verb/ iverb

      call unkfile( 11, logfile )

      do i=1,np
         write(11,2) pmin(i)
      end do
      write(11,3) lamda, ymin

      close(11)

 2    format(g14.7)
 3    format(2g12.6)

      return
      end
