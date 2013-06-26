C     Random: generates a (pseudo)random number.  
C
C     This file is part of StarFISH
C     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
C     
C -=> StarFISH is free software; you can redistribute it 
C -=> and/or modify it under the terms of the GNU General Public 
C -=> License (GPL) as published by the Free Software Foundation;
C -=> either version 2 of the License, or (at your option) any 
C -=> later version.
C
C     Revised Apr 2004 to use Maclaren's DPRAND routine
C     (no longer using a 100-element search/replace array)
C
      subroutine random(a,xl,xu)

      double precision x,a,xl,xu

      x = dprand()
      a = x*(xu - xl) + xl

      return
      end


C   Copyright (C) 1992  N.M. Maclaren
C   Copyright (C) 1992  The University of Cambridge

C   This software may be reproduced and used freely, provided that all
C   users of it agree that the copyright holders are not liable for any
C   damage or injury caused by use of this software and that this
C   condition is passed onto all subsequent recipients of the software,
C   whether modified or not.



      SUBROUTINE SDPRND (ISEED)
      DOUBLE PRECISION XMOD, YMOD, POLY(101), OTHER, OFFSET, X
      PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0)
      INTEGER ISEED, INDEX, IX, IY, IZ, I
      LOGICAL INITAL
      SAVE INITAL
      COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
      DATA INITAL/.TRUE./
C     
C   ISEED should be set to an integer between 0 and 9999 inclusive;
C   a value of 0 will initialise the generator only if it has not
C   already been done.
C
      IF (INITAL .OR. ISEED .NE. 0) THEN
         INITAL = .FALSE.
      ELSE
         RETURN
      END IF
C     
C   INDEX must be initialised to an integer between 1 and 101
C   inclusive, POLY(1...N) to integers between 0 and 1000009710
C   inclusive (not all 0), and OTHER to a non-negative proper fraction
C   with denominator 33554432.  It uses the Wichmann-Hill generator to
C   do this.
C
      IX = MOD(ABS(ISEED),10000)+1
      IY = 2*IX+1
      IZ = 3*IX+1
      DO 10 I = -10,101
         IF (I .GE. 1) POLY(I) = AINT(XMOD*X)
         IX = MOD(171*IX,30269)
         IY = MOD(172*IY,30307)
         IZ = MOD(170*IZ,30323)
         X = MOD(DBLE(IX)/30269.0D0+DBLE(IY)/30307.0D0+
     1        DBLE(IZ)/30323.0D0,1.0D0)
 10   CONTINUE
      OTHER = AINT(YMOD*X)/YMOD
      OFFSET = 1.0D0/YMOD
      INDEX = 1
      END

      DOUBLE PRECISION FUNCTION DPRAND()
      DOUBLE PRECISION XMOD, YMOD, XMOD2, XMOD4, TINY, POLY(101),
     1     OTHER, OFFSET, X, Y
      PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0,
     1     XMOD2 = 2000019422.0D0, XMOD4 = 4000038844.0D0,
     2     TINY = 1.0D-17)
      INTEGER INDEX, N
      LOGICAL INITAL
      SAVE INITAL
      COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
      DATA INITAL/.TRUE./
C     
C   This returns a uniform (0,1) random number, with extremely good
C   uniformity properties.  It assumes that double precision provides
C   at least 33 bits of accuracy, and uses a power of two base.
C
      IF (INITAL) THEN
         CALL SDPRND (0)
         INITAL = .FALSE.
      END IF
C     
C   See [Knuth] for why this implements the algorithm described in
C   the paper.  Note that this code is tuned for machines with fast
C   double precision, but slow multiply and divide; many, many other
C   options are possible.
C
      N = INDEX-64
      IF (N .LE. 0) N = N+101
      X = POLY(INDEX)+POLY(INDEX)
      X = XMOD4-POLY(N)-POLY(N)-X-X-POLY(INDEX)
      IF (X .LT. 0.0D0) THEN
         IF (X .LT. -XMOD) X = X+XMOD2
         IF (X .LT. 0.0D0) X = X+XMOD
      ELSE
         IF (X .GE. XMOD2) THEN
            X = X-XMOD2
            IF (X .GE. XMOD) X = X-XMOD
         END IF
         IF (X .GE. XMOD) X = X-XMOD
      END IF
      POLY(INDEX) = X
      INDEX = INDEX+1
      IF (INDEX .GT. 101) INDEX = INDEX-101
C
C   Add in the second generator modulo 1, and force to be non-zero.
C   The restricted ranges largely cancel themselves out.
C
 10   Y = 37.0D0*OTHER+OFFSET
      OTHER = Y-AINT(Y)
      IF (OTHER .EQ. 0.0D0) GO TO 10
      X = X/XMOD+OTHER
      IF (X .GE. 1.0D0) X = X-1.0D0
      DPRAND = X+TINY
      END
