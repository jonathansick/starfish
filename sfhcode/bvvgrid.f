      subroutine bvvgrid(i,j,n)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     bvvgrid: An example custom gridding of the B-V, V CMD plane.
C     This gridding is designed to be fine where our LMC CMDs are dense, 
C     and coarse where they are sparse.
C
      integer i,j,n

      n = 0 
      
      if (i.gt.6.and.i.le.13.and.j.gt.137) then
         if (j.gt.180) then
            n = 1
         else if (j.gt.170) then
            n = 2
         else if (j.gt.160) then
            n = 3
         else if (j.gt.152) then
            n = 4
         else if (j.gt.144) then
            n = 5
         else 
            n = 6
         endif
      else if (i.gt.6.and.i.le.9.and.j.gt.100) then
         if (j.gt.128) then
            n = 7
         else if (j.gt.120) then
            n = 8
         else if (j.gt.112) then
            n = 9
         else if (j.gt.106) then
            n = 10
         else
            n = 11
         endif
      else if (i.gt.9.and.i.le.13.and.j.gt.100) then
         if (j.gt.128) then
            n = 12
         else if (j.gt.120) then
            n = 13
         else if (j.gt.112) then
            n = 14
         else if (j.gt.106) then
            n = 15
         else
            n = 16
         endif
      else if (i.gt.13.and.i.le.21.and.j.gt.100) then
         if (j.gt.180) then
            n = 17
         else if (j.gt.160) then
            n = 18
         else if (j.gt.144) then
            n = 19
         else if (j.gt.128) then
            n = 20
         else if (j.gt.112) then
            n = 21
         else 
            n = 22
         endif
      else if (i.gt.31.and.j.gt.100) then
         if (j.gt.180) then
            n = 23
         else if (j.gt.160) then
            if (i.le.40) then
               n = 24
            else 
               n = 25
            endif
         else if (j.gt.144) then
            if (i.le.37) then
               n = 26
            else if (i.le.43) then
               n = 27
            else 
               n = 28
            endif
         else if (j.gt.134) then
            if (i.le.37) then
               n = 29
            else if (i.le.43) then
               n = 30
            else 
               n = 31
            endif
         else if (j.gt.124) then
            if (i.le.37) then
               n = 32
            else if (i.le.43) then
               n = 33
            else 
               n = 34
            endif
         else if (j.gt.114) then
            if (i.le.40) then
               n = 35
            else 
               n = 36
            endif
         else
            if (i.gt.46) then
               n = 42
            else if (j.gt.106) then
               if (i.le.40) then
                  n = 37
               else if (j.gt.110) then
                  n = 38
               else 
                  n = 39
               endif
            else
               if (i.le.36) then
                  n = 37
               else if (i.le.41) then
                  n = 40
               else 
                  n = 41
               endif
            endif
         endif
      else if (i.gt.5.and.i.lt.21.and.j.gt.64.and.j.le.100) then
         if (j.gt.96) then
            if (i.le.9) then
               n = 43
            else if (i.le.13) then
               n = 44
            else 
               n = 45
            endif
         else if (j.gt.92) then
            if (i.le.9) then
               n = 46
            else if (i.le.13) then
               n = 47
            else 
               n = 45
            endif
         else if (j.gt.88) then
            if (i.le.9) then
               n = 48
            else if (i.le.13) then
               n = 49
            else 
               n = 50
            endif
         else if (j.gt.84) then
            if (i.le.9) then
               n = 51
            else if (i.le.13) then
               n = 52
            else 
               n = 50
            endif
         else if (j.gt.80) then
            if (i.le.10) then
               n = 53
            else if (i.le.14) then
               n = 54
            else 
               n = 55
            endif
         else if (j.gt.76) then
            if (i.le.10) then
               n = 56
            else if (i.le.14) then
               n = 57
            else 
               n = 55
            endif
         else if (j.gt.68) then
            if (i.le.10) then
               n = 58
            else if (i.le.13) then
               n = 59
            else if (i.le.16) then
               n = 60
            else 
               n = 61
            endif
         else 
            if (i.le.10) then
               n = 62
            else if (i.le.13) then
               n = 63
            else if (i.le.17) then
               n = 64
            else 
               n = 65
            endif
         endif
      else if (i.gt.34.and.i.le.46.and.j.gt.96.and.j.le.100) then
         n = 66 + int((i-35)/4)
      else if (i.gt.32.and.i.le.44.and.j.gt.92.and.j.le.96) then
         n = 69 + int((i-33)/4)
      else if (i.gt.28.and.i.le.42.and.j.gt.88.and.j.le.92) then
         if (i.le.34) then
            n = 72
         else if (i.le.38) then
            n = 73
         else 
            n = 74
         endif
      else if (i.gt.28.and.i.le.42.and.j.gt.84.and.j.le.88) then
         if (i.le.32) then
            n = 75
         else if (i.le.35) then
            n = 76
         else if (i.le.38) then
            n = 77
         else 
            n = 78
         endif
      else if (i.gt.28.and.i.le.40.and.j.gt.80.and.j.le.84) then
         if (i.le.31) then
            n = 79
         else if (i.le.37) then
            n = 80 + int((i-32)/2)
         else 
            n = 83
         endif
      else if (i.gt.28.and.i.le.40.and.j.gt.76.and.j.le.80) then
         if (i.le.31) then
            n = 84
         else if (i.le.37) then
            n = 85 + int((i-32)/2)
         else 
            n = 88
         endif
      else if (i.gt.27.and.i.le.40.and.j.gt.72.and.j.le.76) then
         if (i.le.37) then
            n = 89 + int((i-28)/2)
         else 
            n = 94
         endif
      else if (i.gt.26.and.i.le.40.and.j.gt.68.and.j.le.72) then
         if (i.le.36) then
            n = 95 + int((i-27)/2)
         else 
            n = 100
         endif
      else if (i.gt.26.and.i.le.40.and.j.gt.64.and.j.le.68) then
         if (i.le.34) then
            n = 101 + int((i-27)/2)
         else 
            n = 105
         endif
      else if (i.gt.5.and.i.le.21.and.j.gt.60.and.j.le.64) then
         if (i.le.10) then
            n = 106
         else if (i.le.13) then
            n = 107
         else if (i.le.16) then
            n = 108
         else if (i.le.18) then
            n = 109
         else 
            n = 110
         endif
      else if (i.gt.21.and.i.le.37.and.j.gt.62.and.j.le.64) then
         if (i.le.24) then
            n = 111
         else if (i.le.34) then
            n = 112 + int((i-25)/2)
         else 
            n = 117
         endif
      else if (i.gt.37.and.i.le.40.and.j.gt.52.and.j.le.64) then
         if (j.gt.60) then
            n = 118
         else if (j.gt.56) then
            n = 179
         else 
            n = 180
         endif
      else if (i.gt.21.and.i.le.24.and.j.gt.60.and.j.le.62) then
         n = 119
      else if (i.gt.24.and.i.le.26.and.j.gt.60.and.j.le.62) then
         n = 120
      else if (i.gt.26.and.i.le.32.and.j.gt.54.and.j.le.62) then
         n = 121 + 62 - j + 8*(i-27)
      else if (i.gt.32.and.i.le.34.and.j.gt.52.and.j.le.62) then
         n = 169 + int((62-j)/2)
      else if (i.gt.34.and.i.le.37.and.j.gt.52.and.j.le.62) then
         n = 174 + int((62-j)/2)
      else if (i.gt.5.and.i.le.18.and.j.gt.56.and.j.le.60) then
         if (i.le.10) then
            n = 181
         else if (i.le.13) then
            n = 182
         else if (i.le.16) then
            n = 183
         else 
            n = 184
         endif
      else if (i.gt.18.and.i.le.26.and.j.gt.58.and.j.le.60) then
         n = 185 + int((i-19)/2)
      else if (i.gt.18.and.i.le.26.and.j.gt.56.and.j.le.58) then
         n = 189 + int((i-19)/2)
      else if (i.gt.5.and.i.le.16.and.j.gt.54.and.j.le.56) then
         if (i.le.10) then
            n = 193
         else if (i.le.13) then
            n = 194
         else 
            n = 195
         endif
      else if (i.gt.16.and.i.le.26.and.j.gt.54.and.j.le.56) then
         n = 196 + int((i-17)/2)
      else if (i.gt.5.and.i.le.16.and.j.gt.52.and.j.le.54) then
         if (i.le.10) then
            n = 201
         else if (i.le.13) then
            n = 202
         else 
            n = 203
         endif
      else if (i.gt.16.and.i.le.32.and.j.gt.52.and.j.le.54) then
         n = 204 + int((i-17)/2)
      else if (i.gt.6.and.i.le.10.and.j.gt.49.and.j.le.52) then
         n = 212
      else if (i.gt.10.and.i.le.22.and.j.gt.49.and.j.le.52) then
         n = 213 + int((i-11)/2)
      else if (i.gt.6.and.i.le.10.and.j.gt.46.and.j.le.49) then
         n = 219
      else if (i.gt.10.and.i.le.22.and.j.gt.46.and.j.le.49) then
         n = 220 + int((i-11)/2)
      else if (i.gt.6.and.i.le.10.and.j.gt.43.and.j.le.46) then
         n = 226
      else if (i.gt.10.and.i.le.22.and.j.gt.43.and.j.le.46) then
         n = 227 + int((i-11)/2)
      else if (i.gt.6.and.i.le.10.and.j.gt.40.and.j.le.43) then
         n = 233
      else if (i.gt.10.and.i.le.22.and.j.gt.40.and.j.le.43) then
         n = 234 + int((i-11)/2)
      else if (i.gt.7.and.i.le.11.and.j.gt.37.and.j.le.40) then
         n = 240
      else if (i.gt.11.and.i.le.13.and.j.gt.37.and.j.le.40) then
         n = 241
      else if (i.gt.13.and.i.le.18.and.j.gt.37.and.j.le.40) then
         n = 242 + i - 14
      else if (i.gt.18.and.i.le.22.and.j.gt.37.and.j.le.40) then
         n = 247 + int((i-19)/2)
      else if (i.gt.7.and.i.le.11.and.j.gt.35.and.j.le.37) then
         n = 249
      else if (i.gt.11.and.i.le.13.and.j.gt.35.and.j.le.37) then
         n = 250
      else if (i.gt.13.and.i.le.18.and.j.gt.35.and.j.le.37) then
         n = 251 + i - 14
      else if (i.gt.18.and.i.le.22.and.j.gt.35.and.j.le.37) then
         n = 256 + int((i-19)/2)
      else if (i.gt.7.and.i.le.11.and.j.gt.33.and.j.le.35) then
         n = 258
      else if (i.gt.11.and.i.le.13.and.j.gt.33.and.j.le.35) then
         n = 259
      else if (i.gt.13.and.i.le.18.and.j.gt.33.and.j.le.35) then
         n = 260 + i - 14
      else if (i.gt.18.and.i.le.22.and.j.gt.33.and.j.le.35) then
         n = 265 + int((i-19)/2)
      else if (i.gt.22.and.i.le.26.and.j.gt.33.and.j.le.52) then
         if (j.gt.45) then
            n = 267
         else 
            n = 268 + int((45-j)/6)
         endif
      else if (i.gt.26.and.i.le.28.and.j.gt.33.and.j.le.52) then
         if (j.gt.50) then
            n = 270
         else if (j.gt.48) then
            n = 271
         else
            n = 272 + int((48-j)/3)
         endif
      else if (i.gt.28.and.i.le.30.and.j.gt.33.and.j.le.52) then
         if (j.gt.50) then
            n = 277
         else if (j.gt.48) then
            n = 278
         else
            n = 279 + int((48-j)/3)
         endif
      else if (i.gt.30.and.i.le.32.and.j.gt.33.and.j.le.52) then
         if (j.gt.50) then
            n = 284
         else if (j.gt.48) then
            n = 285
         else
            n = 286 + int((48-j)/3)
         endif
      else if (i.gt.32.and.i.le.35.and.j.gt.33.and.j.le.52) then
         if (j.gt.50) then
            n = 291
         else if (j.gt.48) then
            n = 292
         else
            n = 293 + int((48-j)/3)
         endif
      else if (i.gt.7.and.i.le.12.and.j.gt.11.and.j.le.33) then
         n = 298 + int((33-j)/2)
      else if (i.gt.7.and.i.le.14.and.j.gt.9.and.j.le.11) then
         n = 309
      else if (i.gt.9.and.i.le.14.and.j.le.9) then
         n = 310 + int((9-j)/2)
         if (j.eq.1) n = n - 1
      else if (i.gt.12.and.i.le.14.and.j.gt.11.and.j.le.33) then
         n = 314 + int((33-j)/2)
      else if (i.eq.15.and.j.gt.19.and.j.le.33) then
         n = 325 + 33 - j
      else if (i.eq.15.and.j.gt.11.and.j.le.19) then
         n = 339 + int((19-j)/2)
      else if (i.eq.16.and.j.gt.11.and.j.le.33) then
         n = 343 + 33 - j
      else if (i.gt.16.and.i.le.21.and.j.gt.7.and.j.le.33) then
         n = 365 + 33 - j + 26*(i-17)
      else if (i.eq.22.and.j.gt.31.and.j.le.33) then
         n = 495
      else if (i.eq.22.and.j.gt.7.and.j.le.31) then
         n = 496 + 31 - j
      else if (i.gt.22.and.i.le.24.and.j.gt.27.and.j.le.33) then
         n = 520 + int((33-j)/2)
      else if (i.gt.22.and.i.le.24.and.j.gt.7.and.j.le.27) then
         n = 523 + 27 - j + 20*(i-23)
      else if (i.gt.24.and.i.le.26.and.j.gt.23.and.j.le.33) then
         n = 563 + int((33-j)/2)
      else if (i.gt.24.and.i.le.26.and.j.gt.7.and.j.le.23) then
         n = 568 + 23 - j + 16*(i-25)
      else if (i.gt.26.and.i.le.35.and.j.le.33) then
         n = 600 + int((33-j)/2) + 16*int((i-27)/2)
         if (i.eq.35) n = n - 16
         if (j.eq.1) n = n - 1
      else if (i.gt.14.and.i.le.16.and.j.le.11) then
         n = 664 + int((11-j)/2)
         if (j.eq.1) n = n - 1
      else if (i.gt.16.and.i.le.18.and.j.le.7) then
         n = 669 + int((7-j)/2)
         if (j.eq.1) n = n - 1
      else if (i.gt.18.and.i.le.20.and.j.le.7) then
         n = 672 + int((7-j)/2)
         if (j.eq.1) n = n - 1
      else if (i.gt.20.and.i.le.26.and.j.gt.5.and.j.le.7) then
         n = 675 + i - 21
         if (i.eq.26) n = n - 1
      else if (i.gt.20.and.i.le.26.and.j.gt.3.and.j.le.5) then
         n = 680 + i - 21
         if (j.eq.1) n = n - 1
      else if (i.gt.20.and.i.le.26.and.j.le.3) then
         n = 685 + int((i-21)/2)

c     "NON" boxes:
      else if (j.gt.100.and.i.le.6) then
         n = 703
      else if (j.gt.100) then
         n = 704
      else if (i.le.9) then
         if (j.gt.52) then
            n = 688
         else if (j.gt.40) then
            n = 689
         else if (j.gt.9) then
            n = 690
         else 
            n = 691
         endif
      else if (i.le.28) then
         if (j.gt.76) then
            n = 692
         else if (j.gt.72) then
            n = 693
         else if (j.gt.68) then
            n = 694
         else if (j.gt.64) then
            n = 695
         endif
      else if (i.le.34.and.j.gt.96) then
         n = 696
      else if (i.le.32.and.j.gt.92) then
         n = 697
      else if (i.gt.40) then
         if (j.gt.96) then
            n = 698
         else if (j.gt.92) then
            n = 699
         else if (j.gt.89) then
            n = 700
         else
            n = 701
         endif
      else if (i.gt.35.and.j.le.52) then
         n = 702
      endif

      return
      end

