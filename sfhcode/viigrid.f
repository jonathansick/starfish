      subroutine viigrid(i,j,n)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     Viigrid: An example custom gridding of the V-I, I CMD plane.
C     This gridding is designed to be fine where the LMC CMDs are dense, 
C     and coarse where LMC stars are sparse.
C
      integer i,j,n

      n = 0

      if (i.le.12.and.j.gt.168) then
         n = 1
      else if (i.le.12.and.j.gt.120) then
         n = 2 + int((168-j)/16)
      else if (i.le.12.and.j.gt.100) then
         n = 5 + int((120-j)/10) + 2*int((i-1)/6)
      else if (i.le.21.and.j.gt.168) then
         n = 9
      else if (i.le.21.and.j.gt.136) then
         n = 10
      else if (i.le.21.and.j.gt.100) then
         n = 11
      else if (i.gt.27.and.i.le.43.and.j.gt.176) then
         n = 12
      else if (i.gt.27.and.i.le.35.and.j.gt.128) then
         n = 13 + int((176-j)/8)
      else if (i.gt.35.and.i.le.43.and.j.gt.136) then
         n = 19 + int((176-j)/8)
      else if (i.gt.43.and.j.gt.168) then
         n = 24
      else if (i.gt.43.and.j.gt.152) then
         n = 25
      else if (i.gt.43.and.i.le.50.and.j.gt.136) then
         n = 26 + int((152-j)/8)
      else if (i.gt.50.and.j.gt.136) then
         n = 28
      else if (i.gt.35.and.i.le.47.and.j.gt.128) then
         n = 29 + int((i-36)/4)
      else if (i.gt.27.and.i.le.43.and.j.gt.112) then
         n = 32 + int((i-28)/4) + 4*int((128-j)/8)
      else if (i.gt.27.and.i.le.31.and.j.gt.100) then
         n = 40 + int((112-j)/8)
      else if (i.gt.31.and.i.le.35.and.j.gt.100) then
         n = 42 + int((112-j)/4)
      else if (i.gt.35.and.i.le.39.and.j.gt.100) then
         n = 45 + int((112-j)/6)
      else if (i.gt.4.and.i.le.8.and.j.gt.6) then
         if (j.gt.70) then
            n = 47 + int((100-j)/10)
         else if (j.gt.22) then
            n = 50 + int((70-j)/8)
         else 
            n = 56 + int((22-j)/4)
         endif
      else if (i.gt.8.and.i.le.12.and.j.gt.6) then
         if (j.gt.70) then
            n = 60 + int((100-j)/10)
         else if (j.gt.54) then
            n = 63 + int((70-j)/8)
         else 
            n = 65 + int((54-j)/4)
         endif
      else if (i.gt.12.and.i.le.21.and.j.gt.84.and.j.le.100) then
         n = 77
      else if (i.gt.21.and.i.le.37.and.j.gt.92.and.j.le.100) then
         if (i.le.25) then
            n = 78
         else if (i.gt.25.and.i.le.33) then
            n = 79 + int((i-26)/4) + 2*int((100-j)/4)
         else 
            n = 83
         endif
      else if (i.gt.21.and.i.le.37.and.j.gt.84.and.j.le.92) then
         if (i.le.25) then
            n = 84
         else if (i.gt.25.and.i.le.33) then
            n = 85 + int((i-26)/4) + 2*int((92-j)/4)
         else 
            n = 89
         endif
      else if (i.gt.12.and.i.le.20.and.j.gt.70.and.j.le.84) then
         n = 90
      else if (i.gt.12.and.i.le.20.and.j.gt.62.and.j.le.70) then
         n = 91
      else if (i.gt.12.and.i.le.20.and.j.gt.54.and.j.le.62) then
         n = 92 + int((i-13)/4)
      else if (i.gt.20.and.i.le.24.and.j.gt.68.and.j.le.84) then
         n = 94 + int((84-j)/8)
      else if (i.gt.20.and.i.le.24.and.j.gt.54.and.j.le.68) then
         n = 96
      else if (i.gt.24.and.i.le.26.and.j.gt.68.and.j.le.84) then
         n = 97 + int((84-j)/4)
      else if (i.gt.26.and.i.le.30.and.j.gt.80.and.j.le.84) then
         n = 101 + int((84-j)/2) + 2*int((i-27)/2)
      else if (i.gt.26.and.i.le.30.and.j.gt.72.and.j.le.80) then
         n = 105 + 80 - j + 8*(i-27)
      else if (i.gt.24.and.i.le.30.and.j.gt.68.and.j.le.72) then
         n = 137 + int((72-j)/2) + 2*int((i-27)/2)
      else if (i.gt.30.and.i.le.32.and.j.gt.68.and.j.le.84) then
         n = 141 + int((84-j)/4)
      else if (i.gt.24.and.i.le.32.and.j.gt.62.and.j.le.68) then
         n = 145 + int((68-j)/3) + 2*int((i-25)/4)
      else if (i.gt.24.and.i.le.32.and.j.gt.54.and.j.le.62) then
         n = 149 + int((i-25)/4)
      else if (i.gt.32.and.i.le.40.and.j.gt.54.and.j.le.68) then
         n = 151
      else if (i.gt.12.and.i.le.16.and.j.gt.42.and.j.le.54) then
         n = 152 + int((54-j)/4)
      else if (i.gt.16.and.i.le.32.and.j.gt.46.and.j.le.54) then
         n = 155 + int((i-17)/4)
      else if (i.gt.32.and.i.le.40.and.j.gt.38.and.j.le.54) then
         n = 159
      else if (i.gt.16.and.i.le.24.and.j.gt.42.and.j.le.46) then
         n = 160 + int((i-17)/4)
      else if (i.gt.12.and.i.le.20.and.j.gt.38.and.j.le.42) then
         n = 162 + int((i-13)/2) + 4*int((42-j)/2)
      else if (i.gt.20.and.i.le.24.and.j.gt.38.and.j.le.42) then
         n = 170
      else if (i.gt.24.and.i.le.32.and.j.gt.38.and.j.le.46) then
         n = 171 + int((i-25)/4)
      else if (i.gt.12.and.i.le.14.and.j.gt.36.and.j.le.38) then
         n = 173
      else if (i.gt.12.and.i.le.14.and.j.gt.20.and.j.le.36) then
         n = 174 + 36 - j + 16*(i-13)
      else if (i.gt.14.and.i.le.20.and.j.gt.20.and.j.le.38) then
         n = 206 + 38 - j + 18*(i-15)
      else if (i.gt.20.and.i.le.22.and.j.gt.36.and.j.le.38) then
         n = 314
      else if (i.gt.20.and.i.le.22.and.j.gt.20.and.j.le.36) then
         n = 315 + 36 - j + 16*(i-21)
      else if (i.gt.22.and.i.le.24.and.j.gt.20.and.j.le.38) then
         n = 347 + int((38-j)/2)
      else if (i.gt.24.and.i.le.28.and.j.gt.18.and.j.le.38) then
         n = 356 + int((38-j)/4)
      else if (i.gt.28.and.i.le.32.and.j.gt.18.and.j.le.38) then
         n = 361 + int((38-j)/8)
      else if (i.gt.32.and.i.le.40.and.j.gt.14.and.j.le.38) then
         n = 364
      else if (i.gt.12.and.i.le.32.and.j.gt.6.and.j.le.18) then
         n = 365 + int((i-13)/4) + 5*int((18-j)/4)
         if (n.gt.377) n = 374
      else if (i.gt.8.and.i.le.24.and.j.le.6) then
         n = 378 + int((i-16)/8)
      else if (i.gt.32.and.i.le.40.and.j.gt.68.and.j.le.84) then
         n = 380 + int((84-j)/8)
      else if (i.gt.12.and.i.le.24.and.j.gt.18.and.j.le.20) then
         n = 382 + int((i-13)/2)

c     "NON" boxes
      else if (i.le.27.and.j.gt.100) then
         if (j.gt.168) then
            n = 388
         else if (j.gt.136) then
            n = 389
         else 
            n = 390
         endif
      else if (i.le.47.and.j.gt.112) then
         n = 391
      else if (j.gt.112) then
         n = 392
      else if (j.gt.100) then
         n = 393
      else if (i.le.8) then
         if (j.gt.80) then
            n = 394
         else if (j.gt.54) then
            n = 395
         else if (j.gt.30) then
            n = 396
         else if (j.gt.6) then
            n = 397
         else 
            n = 398
         endif
      else if (j.gt.84) then
         n = 399
      else if (j.gt.54) then
         n = 400
      else if (j.gt.38) then
         n = 401
      else if (i.le.32) then
         n = 404
      else if (i.le.40) then
         n = 403
      else 
         n = 402
      endif

c     ****Get rid of boxes****
c      if (n.eq.381) n = 0

      return
      end
