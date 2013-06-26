      subroutine ubbgrid(i,j,n)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     ubbgrid: An example custom gridding of the U-B, B CMD plane.
C     This gridding is designed to be fine where the LMC CMDs are dense, 
C     and coarse where LMC stars are sparse.
C
      integer i,j,n

      n = 0

      if (i.gt.52.and.j.gt.146) then
         n = 1
      else if (i.gt.52.and.j.gt.128) then
         n = 2
      else if (i.gt.52.and.j.gt.100) then
         n = 3
      else if (i.gt.35.and.j.gt.146) then
         n = 4
      else if (i.gt.35.and.j.gt.128) then
         n = 5
      else if (i.gt.35.and.j.gt.100) then
         n = 6
      else if (i.gt.13.and.i.le.24.and.j.gt.176) then
         n = 7
      else if (i.gt.13.and.i.le.24.and.j.gt.146) then
         n = 8
      else if (i.gt.5.and.i.le.13.and.j.gt.176) then
         n = 9
      else if (i.gt.5.and.i.le.13.and.j.gt.156) then
         n = 10
      else if (i.gt.5.and.i.le.9.and.j.gt.146) then
         n = 11
      else if (i.gt.9.and.i.le.13.and.j.gt.146) then
         n = 12
      else if (i.gt.7.and.i.le.10.and.j.gt.140) then
         n = 13
      else if (i.gt.10.and.i.le.14.and.j.gt.140) then
         n = 14
      else if (i.gt.7.and.i.le.11.and.j.gt.134) then
         n = 15
      else if (i.gt.11.and.i.le.15.and.j.gt.134.and.j.le.140) then
         n = 16
      else if (i.gt.7.and.i.le.11.and.j.gt.128) then
         n = 18
      else if (i.gt.11.and.i.le.17.and.j.gt.128.and.j.le.134) then
         n = 19
      else if (i.gt.14.and.i.le.24.and.j.gt.128) then
         n = 17
      else if (i.gt.7.and.i.le.12.and.j.gt.122) then
         n = 20
      else if (i.gt.12.and.i.le.18.and.j.gt.122) then
         n = 21
      else if (i.gt.7.and.i.le.19.and.j.gt.116.and.j.le.122) then
         if (i.le.10) then
            n = 22
         else if (i.le.13) then
            n = 23
         else if (i.le.16) then
            n = 24
         else 
            n = 25
         endif
      else if (i.gt.7.and.i.le.19.and.j.gt.100.and.j.le.116) then
         if (i.le.10) then
            n = 27
         else if (i.le.13) then
            n = 28
         else if (i.le.16) then
            n = 29
         else 
            n = 30
         endif
      else if (i.gt.18.and.i.le.24.and.j.gt.100) then
         n = 26
      else if (i.gt.9.and.i.le.24.and.j.gt.94) then
         if (i.le.13) then
            n = 31
         else if (i.le.16) then
            n = 32
         else if (i.le.19) then
            n = 33
         else 
            n = 34
         endif
      else if (i.gt.9.and.i.le.24.and.j.gt.88) then
         if (i.le.13) then
            n = 35
         else if (i.le.16) then
            n = 36
         else if (i.le.19) then
            n = 37
         else 
            n = 38
         endif
      else if (i.gt.9.and.i.le.24.and.j.gt.82) then
         if (i.le.13) then
            n = 39
         else if (i.le.16) then
            n = 40
         else if (i.le.19) then
            n = 41
         else 
            n = 42
         endif
      else if (i.gt.54.and.j.gt.78.and.j.le.100) then
         n = 44
      else if (i.gt.35.and.j.gt.78.and.j.le.100) then
         n = 43
      else if (i.gt.12.and.i.le.24.and.j.gt.76) then
         if (i.le.15) then
            n = 45
         else if (i.le.18) then
            n = 46
         else if (i.le.21) then
            n = 47
         else 
            n = 48
         endif
      else if (i.gt.12.and.i.le.24.and.j.gt.70) then
         if (i.le.15) then
            n = 49
         else if (i.le.18) then
            n = 50
         else if (i.le.21) then
            n = 51
         else 
            n = 52
         endif
      else if (i.gt.12.and.i.le.24.and.j.gt.64) then
         if (i.le.15) then
            n = 61
         else if (i.le.18) then
            n = 62
         else if (i.le.21) then
            n = 63
         else 
            n = 64
         endif
      else if (i.gt.12.and.i.le.24.and.j.gt.60) then
         if (i.le.15) then
            n = 77
         else if (i.le.18) then
            n = 78
         else if (i.le.21) then
            n = 79
         else 
            n = 80
         endif
      else if (i.gt.12.and.i.le.24.and.j.gt.56) then
         if (i.le.15) then
            n = 87
         else if (i.le.18) then
            n = 88
         else if (i.le.21) then
            n = 89
         else 
            n = 90
         endif
      else if (i.gt.35.and.i.le.66.and.j.gt.74) then 
         if (i.le.54) then
            n = 53
         else if (i.le.58) then
            n = 54
         else if (i.le.62) then
            n = 55
         else 
            n = 56
         endif
      else if (i.gt.35.and.i.le.66.and.j.gt.70) then 
         if (i.le.50) then 
            n = 53
         else if (i.le.54) then
            n = 57
         else if (i.le.58) then
            n = 58
         else if (i.le.62) then
            n = 59
         else 
            n = 60
         endif
      else if (i.gt.35.and.i.le.42.and.j.gt.64) then
         n = 65
      else if (i.gt.42.and.i.le.66.and.j.gt.66) then
         if (i.le.46) then
            n = 66
         else if (i.le.50) then
            n = 67
         else if (i.le.54) then
            n = 68
         else if (i.le.58) then
            n = 69
         else if (i.le.62) then
            n = 70
         else 
            n = 71
         endif
      else if (i.gt.42.and.i.le.62.and.j.gt.62) then
         if (i.le.46) then
            n = 72
         else if (i.le.50) then
            n = 73
         else if (i.le.54) then
            n = 74
         else if (i.le.58) then
            n = 75
         else 
            n = 76
         endif
      else if (i.gt.34.and.i.le.38.and.j.gt.60.and.j.le.64) then
         n = 81
      else if (i.gt.38.and.i.le.42.and.j.gt.60) then
         n = 82
      else if (i.gt.34.and.i.le.38.and.j.gt.56.and.j.le.60) then
         n = 91
      else if (i.gt.38.and.i.le.42.and.j.gt.56) then
         n = 92
      else if (i.gt.42.and.i.le.58.and.j.gt.58) then
         if (i.le.46) then
            n = 83
         else if (i.le.50) then
            n = 84
         else if (i.le.54) then
            n = 85
         else 
            n = 86
         endif
      else if (i.gt.42.and.i.le.55.and.j.gt.54) then
         if (i.le.46) then
            n = 93
         else if (i.le.50) then
            n = 94
         else 
            n = 95
         endif
      else if (i.gt.18.and.i.le.42.and.j.gt.52.and.j.le.56) then
         if (i.le.21) then
            n = 96
         else if (i.le.24) then
            n = 97
         else if (i.le.27) then
            n = 98
         else if (i.le.30) then
            n = 99
         else if (i.le.33) then
            n = 100
         else if (i.le.36) then
            n = 101
         else if (i.le.39) then
            n = 102
         else 
            n = 103
         endif
      else if (i.gt.42.and.i.le.55.and.j.gt.50) then
         if (i.le.46) then
            n = 104
         else if (i.le.50) then
            n = 105
         else 
            n = 106
         endif
      else if (i.gt.18.and.i.le.42.and.j.gt.49.and.j.le.52) then
         if (i.le.21) then
            n = 107
         else if (i.le.24) then
            n = 108
         else if (i.le.27) then
            n = 109
         else if (i.le.30) then
            n = 110
         else if (i.le.33) then
            n = 111
         else if (i.le.36) then
            n = 112
         else if (i.le.39) then
            n = 113
         else 
            n = 114
         endif
      else if (i.gt.18.and.i.le.50.and.j.gt.46.and.j.le.50) then
         if (i.le.21) then
            n = 115
         else if (i.le.24) then
            n = 116
         else if (i.le.27) then
            n = 117
         else if (i.le.30) then
            n = 118
         else if (i.le.33) then
            n = 119
         else if (i.le.36) then
            n = 120
         else if (i.le.39) then
            n = 121
         else if (i.le.42) then
            n = 122
         else if (i.le.46) then
            n = 123
         else 
            n = 124
         endif
      else if (i.gt.18.and.i.le.21.and.j.gt.43) then
         n = 125
      else if (i.gt.21.and.i.le.24.and.j.gt.43) then
         n = 126
      else if (i.gt.24.and.i.le.36.and.j.gt.44.and.j.le.46) then
         if (i.le.26) then
            n = 127
         else if (i.le.28) then
            n = 128
         else if (i.le.30) then
            n = 129
         else if (i.le.32) then
            n = 130
         else if (i.le.34) then
            n = 131
         else 
            n = 132
         endif
      else if (i.gt.36.and.i.le.45.and.j.gt.43) then
         if (i.le.39) then
            n = 133
         else if (i.le.42) then
            n = 134
         else 
            n = 135
         endif
      else if (i.gt.36.and.i.le.45.and.j.gt.40) then
         if (i.le.39) then
            n = 138
         else if (i.le.42) then
            n = 139
         else 
            n = 140
         endif
      else if (i.gt.45.and.i.le.50.and.j.gt.34) then
         if (j.gt.42) then
            n = 136
         else if (j.gt.38) then
            n = 141
         else 
            n = 148
         endif
      else if (i.gt.34.and.i.le.36.and.j.gt.42.and.j.le.44) then
         n = 137
      else if (i.gt.39.and.i.le.45.and.j.gt.37) then
         if (i.le.42) then
            n = 142
         else 
            n = 143
         endif
      else if (i.gt.18.and.i.le.21.and.j.gt.19) then
         if (j.gt.40) then
            n = 149
         else if (j.gt.37) then
            n = 150
         else if (j.gt.34) then
            n = 151
         else if (j.gt.31) then
            n = 152
         else if (j.gt.28) then
            n = 153
         else if (j.gt.25) then
            n = 154
         else if (j.gt.22) then
            n = 155
         else 
            n = 156
         endif
      else if (i.gt.21.and.i.le.23.and.j.gt.25) then
         if (j.gt.41) then
            n = 157
         else if (j.gt.39) then
            n = 158
         else if (j.gt.37) then
            n = 159
         else if (j.gt.35) then
            n = 160
         else if (j.gt.33) then
            n = 161
         else if (j.gt.31) then
            n = 162
         else if (j.gt.29) then
            n = 163
         else if (j.gt.27) then
            n = 164
         else 
            n = 165
         endif
      else if (i.gt.21.and.i.le.23.and.j.gt.13) then
         n = 166 + (25-j) + 12*(i-22)
      else if (i.eq.24.and.j.gt.13) then
         n = 190 + 43 - j
      else if (i.gt.24.and.i.le.29.and.j.gt.13.and.j.le.44) then
         n = 220 + 44 - j + 31*(i-25)
      else if (i.gt.29.and.i.le.31.and.j.gt.15.and.j.le.44) then
         n = 375 + 44 - j + 29*(i-30)
      else if (i.eq.32.and.j.gt.22.and.j.le.44) then
         n = 433 + 44 - j
      else if (i.gt.32.and.i.le.34.and.j.gt.34.and.j.le.44) then
         n = 455 + 44 - j + 10*(i-33)
      else if (i.gt.34.and.i.le.36.and.j.gt.36.and.j.le.42) then
         n = 475 + 42 - j + 6*(i-35)
      else if (i.gt.36.and.i.le.39.and.j.gt.37) then
         n = 487 + 40 - j + 3*(i-37)
      else if (i.gt.34.and.i.le.45.and.j.gt.34.and.j.le.37) then
         if (i.le.36) then
            n = 144
         else if (i.le.39) then
            n = 145
         else if (i.le.42) then
            n = 146
         else 
            n = 147
         endif
      else if (i.gt.32.and.i.le.34.and.j.gt.22.and.j.le.34) then
         if (j.gt.32) then
            n = 496
         else if (j.gt.30) then
            n = 497
         else if (j.gt.28) then
            n = 498
         else if (j.gt.26) then
            n = 499
         else if (j.gt.24) then
            n = 500
         else 
            n = 501
         endif
      else if (i.gt.34.and.i.le.43.and.j.gt.22.and.j.le.34) then
         n = 502 + int((34-j)/3) + 4*int((i-35)/3)
      else if (i.gt.43.and.i.le.47.and.j.gt.22) then
         n = 514 + int((34-j)/4)
      else if (i.gt.31.and.i.le.33.and.j.gt.16.and.j.le.22) then
         n = 517 + int((22-j)/2)
      else if (i.gt.33.and.i.le.42.and.j.gt.16.and.j.le.22) then
         n = 520 + int((i-34)/3) + 3*int((22-j)/3)
      else if (i.gt.14.and.i.le.18.and.j.gt.19.and.j.le.25) then
         n = 526 + int((25-j)/3)
      else if (i.gt.7.and.i.le.10.and.j.le.10) then
         n = 528 + int((10-j)/3)
         if (j.eq.1) n = n - 1
      else if (i.gt.10.and.i.le.19.and.j.le.19) then
         n = 531 + int((19-j)/3) + 6*int((i-11)/3)
         if (j.eq.1) n = n - 1
      else if (i.gt.19.and.i.le.21.and.j.gt.13.and.j.le.19) then
         n = 549 + int((19-j)/2)
      else if (i.gt.19.and.i.le.25.and.j.le.13) then
         n = 552 + int((13-j)/3) + 4*int((i-20)/3)
         if (j.eq.1) n = n - 1
      else if (i.gt.25.and.i.le.28.and.j.gt.4.and.j.le.13) then
         n = 560 + int((13-j)/3)
      else if (i.gt.28.and.i.le.31.and.j.gt.7.and.j.le.16) then
         n = 563 + int((16-j)/3)
      else if (i.gt.31.and.i.le.34.and.j.gt.13.and.j.le.16) then
         n = 566
      else if (i.gt.31.and.i.le.34.and.j.gt.10.and.j.le.13) then
         n = 567
      else if (i.gt.34.and.i.le.37.and.j.gt.13.and.j.le.16) then
         n = 568

c     "NON" boxes:
      else if (i.le.7.and.j.gt.100) then
         if (j.gt.176) then
            n = 569
         else if (j.gt.156) then
            n = 570
         else if (j.gt.146) then
            n = 571
         else 
            n = 572
         endif
      else if (j.gt.100) then
         n = 573
      else if (i.le.14) then
         if (j.gt.82) then
            n = 574
         else if (j.gt.56) then
            n = 575
         else if (j.gt.36) then
            n = 576
         else if (j.gt.19) then
            n = 577
         else if (j.gt.10) then
            n = 578
         else 
            n = 579
         endif
      else if (i.le.18) then
         n = 580
      else if (i.le.35.and.j.gt.56) then
         if (j.gt.88) then
            n = 581
         else if (j.gt.76) then
            n = 582
         else if (j.gt.64) then
            n = 583
         else 
            n = 584
         endif
      else if (j.gt.66) then
         n = 585
      else if (j.gt.54) then
         if (i.le.58) then
            n = 586
         else if (i.le.62) then
            n = 587
         else 
            n = 588
         endif
      else if (j.gt.34) then
         if (i.le.55) then
            n = 589
         else 
            n = 590
         endif
      else if (i.le.37.and.j.gt.10) then
         n = 591
      else if (i.le.37.and.j.gt.7) then
         n = 592
      else if (i.le.28) then
         n = 593
      else if (i.le.37) then
         n = 594
      else if (i.le.42) then
         n = 595
      else if (i.le.47) then
         n = 596
      else 
         n = 597
      endif

      return
      end




