c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c *** Note: all subroutines herein except dchi and gamma_q are in 
c *** the public domain.
c *** They are from the SLATEC mathematical library, available at:
c *** http://www.netlib.org/slatec/
c

      function dchi(ndof,p)
C     dchi:  Finds the delta-chi**2 value corresponding to the 
C     p confidence interval for any number of degrees of freedom.

      parameter(JMAX=40)
      integer JMAX
      integer ndof,i,j,k
      real p,q,x1,x2,xacc,hdof
      real dchi,f,fmid,xmid

      x1 = 0.1
      x2 = 1000.0
      xacc = 0.001

      q = 1 - p
      hdof = real(ndof)/2.0

      fmid = gamma_q( hdof, x2/2.0 ) - q
      f = gamma_q( hdof, x1/2.0 ) - q

      if (f*fmid.ge.0.) pause 'dchi must be initially bracketed! '

      if (f.lt.0.) then
         dchi = x1
         dx = x2 - x1
      else
         dchi = x2
         dx = x1 - x2
      endif

      do 11 j=1,JMAX
         dx = dx*0.5
         xmid = dchi + dx
         fmid = gamma_q( hdof, xmid/2.0 )
         if (fmid.le.0.) dchi = xmid
         if (abs(dx).lt.xacc .or. fmid.eq.0.) return
 11   continue

      write(*,*) 'too many bisections...'
      write(*,*) dchi
      pause

      return
      end

      function gamma_q( a, x )
c     Return the complement of the incomplete gamma function, 
c     normalized such that as x -> inf, gamma_q -> 1
c     (normalization accomplished by dividing by Gamma(a))
      real a, x, gamma_q

      gamma_q = 1.0 - gami( a, x )/gamma( a ) 

      return
      end


*DECK GAMI
      FUNCTION GAMI (A, X)
C***BEGIN PROLOGUE  GAMI
C***PURPOSE  Evaluate the incomplete Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      SINGLE PRECISION (GAMI-S, DGAMI-D)
C***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate the incomplete gamma function defined by
C
C GAMI = integral from T = 0 to X of EXP(-T) * T**(A-1.0) .
C
C GAMI is evaluated for positive values of A and non-negative values
C of X.  A slight deterioration of 2 or 3 digits accuracy will occur
C when GAMI is very large or very small, because logarithmic variables
C are used.  GAMI, A, and X are single precision.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNGAM, GAMIT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  GAMI
C***FIRST EXECUTABLE STATEMENT  GAMI
!      IF (A .LE. 0.0) CALL XERMSG ('SLATEC', 'GAMI',
!     +   'A MUST BE GT ZERO', 1, 2)
!      IF (X .LT. 0.0) CALL XERMSG ('SLATEC', 'GAMI',
!     +   'X MUST BE GE ZERO', 2, 2)
      if ( a.le.0.0 ) then
         write(*,*) "gami: Fatal Error:  A must be greater than 0.0"
         stop
      endif
      if ( x.lt.0.0 ) then
         write(*,*) "gami: Fatal Error:  X must be greater than "//
     x        "or equal to 0.0"
         stop
      endif
C
      GAMI = 0.0
      IF (X.EQ.0.0) RETURN
C
C THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
      FACTOR = EXP (ALNGAM(A) + A*LOG(X) )
C
      GAMI = FACTOR * GAMIT(A, X)
C
      RETURN
      END


*DECK ALGAMS
      SUBROUTINE ALGAMS (X, ALGAM, SGNGAM)
C***BEGIN PROLOGUE  ALGAMS
C***PURPOSE  Compute the logarithm of the absolute value of the Gamma
C            function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      SINGLE PRECISION (ALGAMS-S, DLGAMS-D)
C***KEYWORDS  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION,
C             FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluates the logarithm of the absolute value of the gamma
C function.
C     X           - input argument
C     ALGAM       - result
C     SGNGAM      - is set to the sign of GAMMA(X) and will
C                   be returned at +1.0 or -1.0.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNGAM
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ALGAMS
C***FIRST EXECUTABLE STATEMENT  ALGAMS
      ALGAM = ALNGAM(X)
      SGNGAM = 1.0
      IF (X.GT.0.0) RETURN
C
      INT = MOD (-AINT(X), 2.0) + 0.1
      IF (INT.EQ.0) SGNGAM = -1.0
C
      RETURN
      END


*DECK ALNGAM
      FUNCTION ALNGAM (X)
C***BEGIN PROLOGUE  ALNGAM
C***PURPOSE  Compute the logarithm of the absolute value of the Gamma
C            function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      SINGLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C)
C***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C ALNGAM(X) computes the logarithm of the absolute value of the
C gamma function at X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  GAMMA, R1MACH, R9LGMC, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  ALNGAM
      LOGICAL FIRST
      EXTERNAL GAMMA
      SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL, FIRST
      DATA SQ2PIL / 0.9189385332 0467274E0/
      DATA SQPI2L / 0.2257913526 4472743E0/
      DATA PI     / 3.1415926535 8979324E0/
      DATA FIRST  /.TRUE./
C***FIRST EXECUTABLE STATEMENT  ALNGAM
      IF (FIRST) THEN
         XMAX = R1MACH(2)/LOG(R1MACH(2))
         DXREL = SQRT (R1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.10.0) GO TO 20
C
C LOG (ABS (GAMMA(X))) FOR  ABS(X) .LE. 10.0
C
      ALNGAM = LOG (ABS (GAMMA(X)))
      RETURN
C
C LOG (ABS (GAMMA(X))) FOR ABS(X) .GT. 10.0
C
! 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'ALNGAM',
!     +   'ABS(X) SO BIG ALNGAM OVERFLOWS', 2, 2)
 20   if ( y.gt.xmax ) then
         write(*,*) "alngam: Overflow Error:  "//
     x     "abs(x) too large."
         stop
      endif
C
      IF (X.GT.0.) ALNGAM = SQ2PIL + (X-0.5)*LOG(X) - X + R9LGMC(Y)
      IF (X.GT.0.) RETURN
C
      SINPIY = ABS (SIN(PI*Y))
!      IF (SINPIY .EQ. 0.) CALL XERMSG ('SLATEC', 'ALNGAM',
!     +   'X IS A NEGATIVE INTEGER', 3, 2)
      if ( sinpiy.eq.0.) then
         write(*,*) "alngam: Fatal Error:  X cannot be negative."
         stop
      endif
C
!      IF (ABS((X-AINT(X-0.5))/X) .LT. DXREL) CALL XERMSG ('SLATEC',
!     +   'ALNGAM', 'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR ' //
!     +   'NEGATIVE INTEGER', 1, 1)
      if (abs((x-aint(x-0.5))/x) .lt. dxrel ) 
     x     write(*,*) "alngam: Warning:  X is very near to a "//
     x     "negative-integer value."
C
      ALNGAM = SQPI2L + (X-0.5)*LOG(Y) - X - LOG(SINPIY) - R9LGMC(Y)
      RETURN
C
      END


*DECK GAMIT
      REAL FUNCTION GAMIT (A, X)
C***BEGIN PROLOGUE  GAMIT
C***PURPOSE  Calculate Tricomi's form of the incomplete Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      SINGLE PRECISION (GAMIT-S, DGAMIT-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
C             SPECIAL FUNCTIONS, TRICOMI
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C   Evaluate Tricomi's incomplete gamma function defined by
C
C   GAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) *
C             T**(A-1.)
C
C   for A .GT. 0.0 and by analytic continuation for A .LE. 0.0.
C   GAMMA(X) is the complete gamma function of X.
C
C   GAMIT is evaluated for arbitrary real values of A and for non-
C   negative values of X (even though GAMIT is defined for X .LT.
C   0.0), except that for X = 0 and A .LE. 0.0, GAMIT is infinite,
C   which is a fatal error.
C
C   The function and both arguments are REAL.
C
C   A slight deterioration of 2 or 3 digits accuracy will occur when
C   GAMIT is very large or very small in absolute value, because log-
C   arithmic variables are used.  Also, if the parameter  A  is very
C   close to a negative integer (but not a negative integer), there is
C   a loss of accuracy, which is reported if the result is less than
C   half machine precision.
C
C***REFERENCES  W. Gautschi, A computational procedure for incomplete
C                 gamma functions, ACM Transactions on Mathematical
C                 Software 5, 4 (December 1979), pp. 466-481.
C               W. Gautschi, Incomplete gamma functions, Algorithm 542,
C                 ACM Transactions on Mathematical Software 5, 4
C                 (December 1979), pp. 482-489.
C***ROUTINES CALLED  ALGAMS, ALNGAM, GAMR, R1MACH, R9GMIT, R9LGIC,
C                    R9LGIT, XERCLR, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
C***END PROLOGUE  GAMIT
      LOGICAL FIRST
      SAVE ALNEPS, SQEPS, BOT, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  GAMIT
      IF (FIRST) THEN
         ALNEPS = -LOG(R1MACH(3))
         SQEPS = SQRT(R1MACH(4))
         BOT = LOG(R1MACH(1))
      ENDIF
      FIRST = .FALSE.
C
!      IF (X .LT. 0.0) CALL XERMSG ('SLATEC', 'GAMIT', 'X IS NEGATIVE',
!     +   2, 2)
      if ( x.lt.0.0) then
         write(*,*) "gamit: Fatal Error:  X is negative."
         stop
      endif
C
      IF (X.NE.0.0) ALX = LOG(X)
      SGA = 1.0
      IF (A.NE.0.0) SGA = SIGN (1.0, A)
      AINTA = AINT (A+0.5*SGA)
      AEPS = A - AINTA
C
      IF (X.GT.0.0) GO TO 20
      GAMIT = 0.0
      IF (AINTA.GT.0.0 .OR. AEPS.NE.0.0) GAMIT = GAMR(A+1.0)
      RETURN
C
 20   IF (X.GT.1.0) GO TO 40
      IF (A.GE.(-0.5) .OR. AEPS.NE.0.0) CALL ALGAMS (A+1.0, ALGAP1,
     1  SGNGAM)
      GAMIT = R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      RETURN
C
 40   IF (A.LT.X) GO TO 50
      T = R9LGIT (A, X, ALNGAM(A+1.0))
!      IF (T.LT.BOT) CALL XERCLR
      GAMIT = EXP(T)
      RETURN
C
 50   ALNG = R9LGIC (A, X, ALX)
C
C EVALUATE GAMIT IN TERMS OF LOG(GAMIC(A,X))
C
      H = 1.0
      IF (AEPS.EQ.0.0 .AND. AINTA.LE.0.0) GO TO 60
      CALL ALGAMS (A+1.0, ALGAP1, SGNGAM)
      T = LOG(ABS(A)) + ALNG - ALGAP1
      IF (T.GT.ALNEPS) GO TO 70
      IF (T.GT.(-ALNEPS)) H = 1.0 - SGA*SGNGAM*EXP(T)
      IF (ABS(H).GT.SQEPS) GO TO 60
!      CALL XERCLR
!      CALL XERMSG ('SLATEC', 'GAMIT', 'RESULT LT HALF PRECISION', 1, 1)
      write(*,*) "gamit: Warning:  result is less than half precision."
C
 60   T = -A*ALX + LOG(ABS(H))
!      IF (T.LT.BOT) CALL XERCLR
      GAMIT = SIGN (EXP(T), H)
      RETURN
C
 70   T = T - A*ALX
!      IF (T.LT.BOT) CALL XERCLR
      GAMIT = -SGA*SGNGAM*EXP(T)
      RETURN
C
      END


*DECK GAMMA
      FUNCTION GAMMA (X)
C***BEGIN PROLOGUE  GAMMA
C***PURPOSE  Compute the complete Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      SINGLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C GAMMA computes the gamma function at X, where X is not 0, -1, -2, ....
C GAMMA and X are single precision.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, GAMLIM, INITS, R1MACH, R9LGMC, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  GAMMA
      DIMENSION GCS(23)
      LOGICAL FIRST
      SAVE GCS, PI, SQ2PIL, NGCS, XMIN, XMAX, DXREL, FIRST
      DATA GCS   ( 1) / .0085711955 90989331E0/
      DATA GCS   ( 2) / .0044153813 24841007E0/
      DATA GCS   ( 3) / .0568504368 1599363E0/
      DATA GCS   ( 4) /-.0042198353 96418561E0/
      DATA GCS   ( 5) / .0013268081 81212460E0/
      DATA GCS   ( 6) /-.0001893024 529798880E0/
      DATA GCS   ( 7) / .0000360692 532744124E0/
      DATA GCS   ( 8) /-.0000060567 619044608E0/
      DATA GCS   ( 9) / .0000010558 295463022E0/
      DATA GCS   (10) /-.0000001811 967365542E0/
      DATA GCS   (11) / .0000000311 772496471E0/
      DATA GCS   (12) /-.0000000053 542196390E0/
      DATA GCS   (13) / .0000000009 193275519E0/
      DATA GCS   (14) /-.0000000001 577941280E0/
      DATA GCS   (15) / .0000000000 270798062E0/
      DATA GCS   (16) /-.0000000000 046468186E0/
      DATA GCS   (17) / .0000000000 007973350E0/
      DATA GCS   (18) /-.0000000000 001368078E0/
      DATA GCS   (19) / .0000000000 000234731E0/
      DATA GCS   (20) /-.0000000000 000040274E0/
      DATA GCS   (21) / .0000000000 000006910E0/
      DATA GCS   (22) /-.0000000000 000001185E0/
      DATA GCS   (23) / .0000000000 000000203E0/
      DATA PI /3.14159 26535 89793 24E0/
C SQ2PIL IS LOG (SQRT (2.*PI) )
      DATA SQ2PIL /0.91893 85332 04672 74E0/
      DATA FIRST /.TRUE./
C
C LANL DEPENDENT CODE REMOVED 81.02.04
C
C***FIRST EXECUTABLE STATEMENT  GAMMA
      IF (FIRST) THEN
C
C ---------------------------------------------------------------------
C INITIALIZE.  FIND LEGAL BOUNDS FOR X, AND DETERMINE THE NUMBER OF
C TERMS IN THE SERIES REQUIRED TO ATTAIN AN ACCURACY TEN TIMES BETTER
C THAN MACHINE PRECISION.
C
         NGCS = INITS (GCS, 23, 0.1*R1MACH(3))
C
         CALL GAMLIM (XMIN, XMAX)
         DXREL = SQRT (R1MACH(4))
C
C ---------------------------------------------------------------------
C FINISH INITIALIZATION.  START EVALUATING GAMMA(X).
C
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.10.0) GO TO 50
C
C COMPUTE GAMMA(X) FOR ABS(X) .LE. 10.0.  REDUCE INTERVAL AND
C FIND GAMMA(1+Y) FOR 0. .LE. Y .LT. 1. FIRST OF ALL.
C
      N = X
      IF (X.LT.0.) N = N - 1
      Y = X - N
      N = N - 1
      GAMMA = 0.9375 + CSEVL(2.*Y-1., GCS, NGCS)
      IF (N.EQ.0) RETURN
C
      IF (N.GT.0) GO TO 30
C
C COMPUTE GAMMA(X) FOR X .LT. 1.
C
      N = -N
!      IF (X .EQ. 0.) CALL XERMSG ('SLATEC', 'GAMMA', 'X IS 0', 4, 2)
!      IF (X .LT. 0. .AND. X+N-2 .EQ. 0.) CALL XERMSG ('SLATEC', 'GAMMA'
!     1, 'X IS A NEGATIVE INTEGER', 4, 2)
!      IF (X .LT. (-0.5) .AND. ABS((X-AINT(X-0.5))/X) .LT. DXREL) CALL
!     1XERMSG ( 'SLATEC', 'GAMMA',
!     2'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER'
!     3, 1, 1)
      if ( x.eq.0.) then
         write(*,*) "gamma: Fatal Error:  X is 0.0"
         stop
      endif
      if ( x.lt.0. .and. x+n-2 .eq.0.) then
         write(*,*) "gamma: Fatal Error:  X is a negative integer"
         stop
      endif
      if ( x.lt.-0.5 .and. abs((x-aint(x-0.5))/x) .lt. dxrel ) 
     x     write(*,*) "gamma: Warning:  X is near negative integer"

C
      DO 20 I=1,N
        GAMMA = GAMMA / (X+I-1)
 20   CONTINUE
      RETURN
C
C GAMMA(X) FOR X .GE. 2.
C
 30   DO 40 I=1,N
        GAMMA = (Y+I)*GAMMA
 40   CONTINUE
      RETURN
C
C COMPUTE GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
C
! 50   IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'GAMMA',
!     +   'X SO BIG GAMMA OVERFLOWS', 3, 2)
 50   if ( x.gt.xmax ) then
         write(*,*) "gamma: Overflow Error:  X too large"
         stop
      endif
C     
      GAMMA = 0.
!      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'GAMMA',
!     +   'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
!      IF (X.LT.XMIN) RETURN
      if ( x.lt.xmin ) then
         write(*,*) "gamma: Underflow Error:  X too small"
         return
      endif
C
      GAMMA = EXP((Y-0.5)*LOG(Y) - Y + SQ2PIL + R9LGMC(Y) )
      IF (X.GT.0.) RETURN
C
!      IF (ABS((X-AINT(X-0.5))/X) .LT. DXREL) CALL XERMSG ('SLATEC',
!     +   'GAMMA',
!     +   'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
      if ( abs((x-aint(x-0.5))/x) .lt. dxrel ) 
     x     write(*,*) "gamma: Warning:  X near negative integer value"
C
      SINPIY = SIN (PI*Y)
!      IF (SINPIY .EQ. 0.) CALL XERMSG ('SLATEC', 'GAMMA',
!     +   'X IS A NEGATIVE INTEGER', 4, 2)
      if ( sinpiy .eq. 0. ) then
         write(*,*) "gamma: Fatal Error:  X is a negative integer"
         stop
      endif
C
      GAMMA = -PI / (Y*SINPIY*GAMMA)
C
      RETURN
      END


*DECK GAMR
      FUNCTION GAMR (X)
C***BEGIN PROLOGUE  GAMR
C***PURPOSE  Compute the reciprocal of the Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      SINGLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C)
C***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C GAMR is a single precision function that evaluates the reciprocal
C of the gamma function for single precision argument X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALGAMS, GAMMA, XERCLR, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  GAMR
      EXTERNAL GAMMA
C***FIRST EXECUTABLE STATEMENT  GAMR
      GAMR = 0.0
      IF (X.LE.0.0 .AND. AINT(X).EQ.X) RETURN
C
!      CALL XGETF (IROLD)
!      CALL XSETF (1)
      IF (ABS(X).GT.10.0) GO TO 10
      GAMR = 1.0/GAMMA(X)
!      CALL XERCLR
!      CALL XSETF (IROLD)
      RETURN
C
 10   CALL ALGAMS (X, ALNGX, SGNGX)
!      CALL XERCLR
!      CALL XSETF (IROLD)
      GAMR = SGNGX * EXP(-ALNGX)
      RETURN
C
      END


*DECK GAMLIM
      SUBROUTINE GAMLIM (XMIN, XMAX)
C***BEGIN PROLOGUE  GAMLIM
C***PURPOSE  Compute the minimum and maximum bounds for the argument in
C            the Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A, R2
C***TYPE      SINGLE PRECISION (GAMLIM-S, DGAMLM-D)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Calculate the minimum and maximum legal bounds for X in GAMMA(X).
C XMIN and XMAX are not the only bounds, but they are the only non-
C trivial ones to calculate.
C
C             Output Arguments --
C XMIN   minimum legal value of X in GAMMA(X).  Any smaller value of
C        X might result in underflow.
C XMAX   maximum legal value of X in GAMMA(X).  Any larger value will
C        cause overflow.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  GAMLIM
C***FIRST EXECUTABLE STATEMENT  GAMLIM
      ALNSML = LOG(R1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5)*XLN - XMIN - 0.2258 + ALNSML)
     1    / (XMIN*XLN + 0.5)
        IF (ABS(XMIN-XOLD).LT.0.005) GO TO 20
 10   CONTINUE
!      CALL XERMSG ('SLATEC', 'GAMLIM', 'UNABLE TO FIND XMIN', 1, 2)
      write(*,*) "gamlim: Fatal Error:  Unable to find xmin"
      stop
C
 20   XMIN = -XMIN + 0.01
C
      ALNBIG = LOG(R1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5)*XLN - XMAX + 0.9189 - ALNBIG)
     1    / (XMAX*XLN - 0.5)
        IF (ABS(XMAX-XOLD).LT.0.005) GO TO 40
 30   CONTINUE
!      CALL XERMSG ('SLATEC', 'GAMLIM', 'UNABLE TO FIND XMAX', 2, 2)
      write(*,*) "gamlim: Fatal Error:  Unable to find xmax"
      stop
C
 40   XMAX = XMAX - 0.01
      XMIN = MAX (XMIN, -XMAX+1.)
C
      RETURN
      END


*DECK CSEVL
      FUNCTION CSEVL (X, CS, N)
C***BEGIN PROLOGUE  CSEVL
C***PURPOSE  Evaluate a Chebyshev series.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      SINGLE PRECISION (CSEVL-S, DCSEVL-D)
C***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Evaluate the N-term Chebyshev series CS at X.  Adapted from
C  a method presented in the paper by Broucke referenced below.
C
C       Input Arguments --
C  X    value at which the series is to be evaluated.
C  CS   array of N terms of a Chebyshev series.  In evaluating
C       CS, only half the first coefficient is summed.
C  N    number of terms in array CS.
C
C***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
C                 Chebyshev series, Algorithm 446, Communications of
C                 the A.C.M. 16, (1973) pp. 254-256.
C               L. Fox and I. B. Parker, Chebyshev Polynomials in
C                 Numerical Analysis, Oxford University Press, 1968,
C                 page 56.
C***ROUTINES CALLED  R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900329  Prologued revised extensively and code rewritten to allow
C           X to be slightly outside interval (-1,+1).  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CSEVL
      REAL B0, B1, B2, CS(*), ONEPL, TWOX, X
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  CSEVL
      IF (FIRST) ONEPL = 1.0E0 + R1MACH(4)
      FIRST = .FALSE.
!      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'CSEVL',
!     +   'NUMBER OF TERMS .LE. 0', 2, 2)
!      IF (N .GT. 1000) CALL XERMSG ('SLATEC', 'CSEVL',
!     +   'NUMBER OF TERMS .GT. 1000', 3, 2)
!      IF (ABS(X) .GT. ONEPL) CALL XERMSG ('SLATEC', 'CSEVL',
!     +   'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
      if ( n.lt.1 ) then
         write(*,*) "csevl: Fatal Error:  Num. of terms less than 1"
         stop
      endif
      if ( n.gt.1000) then 
         write(*,*) "csevl: Fatal Error:  Num. of terms greater "//
     x        "than 1000"
         stop
      endif
      if ( abs(x).gt.onepl ) 
     x     write(*,*) "csevl: Error:  X outside interval (-1,+1)"
C
      B1 = 0.0E0
      B0 = 0.0E0
      TWOX = 2.0*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
C
      CSEVL = 0.5E0*(B0-B2)
C
      RETURN
      END


*DECK INITS
      FUNCTION INITS (OS, NOS, ETA)
C***BEGIN PROLOGUE  INITS
C***PURPOSE  Determine the number of terms needed in an orthogonal
C            polynomial series so that it meets a specified accuracy.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      SINGLE PRECISION (INITS-S, INITDS-D)
C***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
C             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Initialize the orthogonal series, represented by the array OS, so
C  that INITS is the number of terms needed to insure the error is no
C  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
C  machine precision.
C
C             Input Arguments --
C   OS     single precision array of NOS coefficients in an orthogonal
C          series.
C   NOS    number of coefficients in OS.
C   ETA    single precision scalar containing requested accuracy of
C          series.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891115  Modified error message.  (WRB)
C   891115  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  INITS
      REAL OS(*)
C***FIRST EXECUTABLE STATEMENT  INITS
!      IF (NOS .LT. 1) CALL XERMSG ('SLATEC', 'INITS',
!     +   'Number of coefficients is less than 1', 2, 1)
      if ( nos.lt.1 ) write(*,*) "inits: Error:  Number of "//
     x     "coefficients is less than 1"
C
      ERR = 0.
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(OS(I))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
C
!   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'INITS',
!     +   'Chebyshev series too short for specified accuracy', 1, 1)
 20   if ( i.eq.nos ) 
     x     write(*,*) "inits: Error:  Chebyshev series too short"
      INITS = I
C
      RETURN
      END


*DECK R1MACH
      REAL FUNCTION R1MACH (I)
C***BEGIN PROLOGUE  R1MACH
C***PURPOSE  Return floating point machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   R1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument, and can be referenced as follows:
C
C        A = R1MACH(I)
C
C   where I=1,...,5.  The (output) value of A above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
C   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   R1MACH(3) = B**(-T), the smallest relative spacing.
C   R1MACH(4) = B**(1-T), the largest relative spacing.
C   R1MACH(5) = LOG10(B)
C
C   Assume single precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(11) = T, the number of base-B digits.
C   I1MACH(12) = EMIN, the smallest exponent E.
C   I1MACH(13) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890213  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added CONVEX -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C***END PROLOGUE  R1MACH
C
!      INTEGER SMALL(2)
!      INTEGER LARGE(2)
!      INTEGER RIGHT(2)
!      INTEGER DIVER(2)
!      INTEGER LOG10(2)
C
      REAL RMACH(5)
      SAVE RMACH
C
!      EQUIVALENCE (RMACH(1),SMALL(1))
!      EQUIVALENCE (RMACH(2),LARGE(1))
!      EQUIVALENCE (RMACH(3),RIGHT(1))
!      EQUIVALENCE (RMACH(4),DIVER(1))
!      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7EFFFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1) / 16#00800000 /
C     DATA LARGE(1) / 16#7FFFFFFF /
C     DATA RIGHT(1) / 16#33800000 /
C     DATA DIVER(1) / 16#34000000 /
C     DATA LOG10(1) / 16#3E9A209B /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA RMACH(1) / Z400800000 /
C     DATA RMACH(2) / Z5FFFFFFFF /
C     DATA RMACH(3) / Z4E9800000 /
C     DATA RMACH(4) / Z4EA800000 /
C     DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
C
C     DATA RMACH(1) / O1771000000000000 /
C     DATA RMACH(2) / O0777777777777777 /
C     DATA RMACH(3) / O1311000000000000 /
C     DATA RMACH(4) / O1301000000000000 /
C     DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA RMACH(1) / Z"3001800000000000" /
C     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
C     DATA RMACH(3) / Z"3FD2800000000000" /
C     DATA RMACH(4) / Z"3FD3800000000000" /
C     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA RMACH(1) / 00564000000000000000B /
C     DATA RMACH(2) / 37767777777777777776B /
C     DATA RMACH(3) / 16414000000000000000B /
C     DATA RMACH(4) / 16424000000000000000B /
C     DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn COMPILER OPTION
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7FFFFFFF' /
C     DATA RMACH(3) / Z'34800000' /
C     DATA RMACH(4) / Z'35000000' /
C     DATA RMACH(5) / Z'3F9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 OR -pd8 COMPILER OPTION
C
C     DATA RMACH(1) / Z'0010000000000000' /
C     DATA RMACH(2) / Z'7FFFFFFFFFFFFFFF' /
C     DATA RMACH(3) / Z'3CC0000000000000' /
C     DATA RMACH(4) / Z'3CD0000000000000' /
C     DATA RMACH(5) / Z'3FF34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CRAY
C
C     DATA RMACH(1) / 200034000000000000000B /
C     DATA RMACH(2) / 577767777777777777776B /
C     DATA RMACH(3) / 377224000000000000000B /
C     DATA RMACH(4) / 377234000000000000000B /
C     DATA RMACH(5) / 377774642023241175720B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC RMACH(5)
C
C     DATA SMALL /    20K,       0 /
C     DATA LARGE / 77777K, 177777K /
C     DATA RIGHT / 35420K,       0 /
C     DATA DIVER / 36020K,       0 /
C     DATA LOG10 / 40423K,  42023K /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA RMACH(1) / '00000080'X /
C     DATA RMACH(2) / 'FFFF7FFF'X /
C     DATA RMACH(3) / '00003480'X /
C     DATA RMACH(4) / '00003500'X /
C     DATA RMACH(5) / '209B3F9A'X /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FLOAT
C
C     DATA RMACH(1) / '00800000'X /
C     DATA RMACH(2) / '7F7FFFFF'X /
C     DATA RMACH(3) / '33800000'X /
C     DATA RMACH(4) / '34000000'X /
C     DATA RMACH(5) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1) /       128 /
C     DATA LARGE(1) /    -32769 /
C     DATA RIGHT(1) /     13440 /
C     DATA DIVER(1) /     13568 /
C     DATA LOG10(1) / 547045274 /
C
C     DATA SMALL(1) / Z00000080 /
C     DATA LARGE(1) / ZFFFF7FFF /
C     DATA RIGHT(1) / Z00003480 /
C     DATA DIVER(1) / Z00003500 /
C     DATA LOG10(1) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*4 IS THE DEFAULT REAL)
C
C     DATA SMALL(1) / '00800000'X /
C     DATA LARGE(1) / '7F7FFFFF'X /
C     DATA RIGHT(1) / '33800000'X /
C     DATA DIVER(1) / '34000000'X /
C     DATA LOG10(1) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
C     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA RMACH(1) / O402400000000 /
C     DATA RMACH(2) / O376777777777 /
C     DATA RMACH(3) / O714400000000 /
C     DATA RMACH(4) / O716400000000 /
C     DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1) / 00004000000B /
C     DATA LARGE(1) / 17677777777B /
C     DATA RIGHT(1) / 06340000000B /
C     DATA DIVER(1) / 06400000000B /
C     DATA LOG10(1) / 07646420233B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA RMACH(1) / Z00100000 /
C     DATA RMACH(2) / Z7FFFFFFF /
C     DATA RMACH(3) / Z3B100000 /
C     DATA RMACH(4) / Z3C100000 /
C     DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
!     DATA SMALL(1) / 1.18E-38      /
!     DATA LARGE(1) / 3.40E+38      /
!     DATA RIGHT(1) / 0.595E-07     /
!     DATA DIVER(1) / 1.19E-07      /
!     DATA LOG10(1) / 0.30102999566 /
      rmach(1) = 1.18E-38
      rmach(2) = 3.40E+38
      rmach(3) = 0.595E-07
      rmach(4) = 1.19E-07
      rmach(5) = 0.30102999566
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
C
C     DATA RMACH(1) / "000400000000 /
C     DATA RMACH(2) / "377777777777 /
C     DATA RMACH(3) / "146400000000 /
C     DATA RMACH(4) / "147400000000 /
C     DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1) /    8388608 /
C     DATA LARGE(1) / 2147483647 /
C     DATA RIGHT(1) /  880803840 /
C     DATA DIVER(1) /  889192448 /
C     DATA LOG10(1) / 1067065499 /
C
C     DATA RMACH(1) / O00040000000 /
C     DATA RMACH(2) / O17777777777 /
C     DATA RMACH(3) / O06440000000 /
C     DATA RMACH(4) / O06500000000 /
C     DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /   128,     0 /
C     DATA LARGE(1), LARGE(2) / 32767,    -1 /
C     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
C     DATA DIVER(1), DIVER(2) / 13568,     0 /
C     DATA LOG10(1), LOG10(2) / 16282,  8347 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
C     DATA DIVER(1), DIVER(2) / O032400, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA RMACH(1) / Z'0010000000000000' /
C     DATA RMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA RMACH(3) / Z'3CA0000000000000' /
C     DATA RMACH(4) / Z'3CB0000000000000' /
C     DATA RMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
C
C     DATA RMACH(1) / O000400000000 /
C     DATA RMACH(2) / O377777777777 /
C     DATA RMACH(3) / O146400000000 /
C     DATA RMACH(4) / O147400000000 /
C     DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA SMALL(1), SMALL(2) /     0,    256/
C     DATA LARGE(1), LARGE(2) /    -1,   -129/
C     DATA RIGHT(1), RIGHT(2) /     0,  26880/
C     DATA DIVER(1), DIVER(2) /     0,  27136/
C     DATA LOG10(1), LOG10(2) /  8347,  32538/
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
!      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'R1MACH',
!     +   'I OUT OF BOUNDS', 1, 2)
      if ( i.lt.1.or.i.gt.5 ) then
         write(*,*) "r1mach: Fatal Error:  i out of bounds"
         stop
      endif
C
      R1MACH = RMACH(I)
      RETURN
C
      END


*DECK R9LGMC
      FUNCTION R9LGMC (X)
C***BEGIN PROLOGUE  R9LGMC
C***SUBSIDIARY
C***PURPOSE  Compute the log Gamma correction factor so that
C            LOG(GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X
C            + R9LGMC(X).
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      SINGLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
C             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the log gamma correction factor for X .GE. 10.0 so that
C  LOG (GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + R9LGMC(X)
C
C Series for ALGM       on the interval  0.          to  1.00000D-02
C                                        with weighted error   3.40E-16
C                                         log weighted error  15.47
C                               significant figures required  14.39
C                                    decimal places required  15.86
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  R9LGMC
      DIMENSION ALGMCS(6)
      LOGICAL FIRST
      SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
      DATA ALGMCS( 1) /    .1666389480 45186E0 /
      DATA ALGMCS( 2) /   -.0000138494 817606E0 /
      DATA ALGMCS( 3) /    .0000000098 108256E0 /
      DATA ALGMCS( 4) /   -.0000000000 180912E0 /
      DATA ALGMCS( 5) /    .0000000000 000622E0 /
      DATA ALGMCS( 6) /   -.0000000000 000003E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  R9LGMC
      IF (FIRST) THEN
         NALGM = INITS (ALGMCS, 6, R1MACH(3))
         XBIG = 1.0/SQRT(R1MACH(3))
         XMAX = EXP (MIN(LOG(R1MACH(2)/12.0), -LOG(12.0*R1MACH(1))) )
      ENDIF
      FIRST = .FALSE.
C
!      IF (X .LT. 10.0) CALL XERMSG ('SLATEC', 'R9LGMC',
!     +   'X MUST BE GE 10', 1, 2)
      if ( x.lt.10.0 ) then
         write(*,*) "r9lgmc: Fatal Error:  X is less than 10."
         stop
      endif
      IF (X.GE.XMAX) GO TO 20
C
      R9LGMC = 1.0/(12.0*X)
      IF (X.LT.XBIG) R9LGMC = CSEVL (2.0*(10./X)**2-1., ALGMCS, NALGM)/X
      RETURN
C
 20   R9LGMC = 0.0
!      CALL XERMSG ('SLATEC', 'R9LGMC', 'X SO BIG R9LGMC UNDERFLOWS', 2,
!     +   1)
      write(*,*) "r9lgmc: Underflow Error:  X too large"
      RETURN
C
      END


*DECK R9GMIT
      FUNCTION R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
C***BEGIN PROLOGUE  R9GMIT
C***SUBSIDIARY
C***PURPOSE  Compute Tricomi's incomplete Gamma function for small
C            arguments.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      SINGLE PRECISION (R9GMIT-S, D9GMIT-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
C             SPECIAL FUNCTIONS, TRICOMI
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute Tricomi's incomplete gamma function for small X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNGAM, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  R9GMIT
      SAVE EPS, BOT
      DATA EPS, BOT / 2*0.0 /
C***FIRST EXECUTABLE STATEMENT  R9GMIT
      IF (EPS.EQ.0.0) EPS = 0.5*R1MACH(3)
      IF (BOT.EQ.0.0) BOT = LOG(R1MACH(1))
C
!      IF (X .LE. 0.0) CALL XERMSG ('SLATEC', 'R9GMIT',
!     +   'X SHOULD BE GT 0', 1, 2)
      if ( x.le.0.0) then
         write(*,*) "r9gmit: Error:  X <= 0"
         stop
      endif
C
      MA = A + 0.5
      IF (A.LT.0.0) MA = A - 0.5
      AEPS = A - MA
C
      AE = A
      IF (A.LT.(-0.5)) AE = AEPS
C
      T = 1.0
      TE = AE
      S = T
      DO 20 K=1,200
        FK = K
        TE = -X*TE/FK
        T = TE/(AE+FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 30
 20   CONTINUE
!      CALL XERMSG ('SLATEC', 'R9GMIT',
!     +   'NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES', 2, 2)
      write(*,*) "r9gmit: Fatal Error:  no convergence after 200 "//
     x     "terms of the Taylor-S series."
      stop
C
 30   IF (A.GE.(-0.5)) ALGS = -ALGAP1 + LOG(S)
      IF (A.GE.(-0.5)) GO TO 60
C
      ALGS = -ALNGAM(1.0+AEPS) + LOG(S)
      S = 1.0
      M = -MA - 1
      IF (M.EQ.0) GO TO 50
      T = 1.0
      DO 40 K=1,M
        T = X*T/(AEPS-M-1+K)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40   CONTINUE
C
 50   R9GMIT = 0.0
      ALGS = -MA*LOG(X) + ALGS
      IF (S.EQ.0.0 .OR. AEPS.EQ.0.0) GO TO 60
C
      SGNG2 = SGNGAM*SIGN(1.0,S)
      ALG2 = -X - ALGAP1 + LOG(ABS(S))
C
      IF (ALG2.GT.BOT) R9GMIT = SGNG2*EXP(ALG2)
      IF (ALGS.GT.BOT) R9GMIT = R9GMIT + EXP(ALGS)
      RETURN
C
 60   R9GMIT = EXP(ALGS)
      RETURN
C
      END


*DECK R9LGIC
      FUNCTION R9LGIC (A, X, ALX)
C***BEGIN PROLOGUE  R9LGIC
C***SUBSIDIARY
C***PURPOSE  Compute the log complementary incomplete Gamma function
C            for large X and for A .LE. X.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      SINGLE PRECISION (R9LGIC-S, D9LGIC-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
C             LOGARITHM, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the log complementary incomplete gamma function for large X
C and for A .LE. X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  R9LGIC
      SAVE EPS
      DATA EPS / 0.0 /
C***FIRST EXECUTABLE STATEMENT  R9LGIC
      IF (EPS.EQ.0.0) EPS = 0.5*R1MACH(3)
C
      XPA = X + 1.0 - A
      XMA = X - 1.0 - A
C
      R = 0.0
      P = 1.0
      S = P
      DO 10 K=1,200
        FK = K
        T = FK*(A-FK)*(1.0+R)
        R = -T/((XMA+2.0*FK)*(XPA+2.0*FK)+T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
!      CALL XERMSG ('SLATEC', 'R9LGIC',
!     +   'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 1, 2)
      write(*,*) "r9lgic: Fatal Error:  No convergence after 200 "//
     x     "terms of the continued fraction."
      stop
C
 20   R9LGIC = A*ALX - X + LOG(S/XPA)
C
      RETURN
      END


*DECK R9LGIT
      FUNCTION R9LGIT (A, X, ALGAP1)
C***BEGIN PROLOGUE  R9LGIT
C***SUBSIDIARY
C***PURPOSE  Compute the logarithm of Tricomi's incomplete Gamma
C            function with Perron's continued fraction for large X and
C            A .GE. X.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      SINGLE PRECISION (R9LGIT-S, D9LGIT-D)
C***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
C             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the log of Tricomi's incomplete gamma function with Perron's
C continued fraction for large X and for A .GE. X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  R9LGIT
      SAVE EPS, SQEPS
      DATA EPS, SQEPS / 2*0.0 /
C***FIRST EXECUTABLE STATEMENT  R9LGIT
      IF (EPS.EQ.0.0) EPS = 0.5*R1MACH(3)
      IF (SQEPS.EQ.0.0) SQEPS = SQRT(R1MACH(4))
C
!      IF (X .LE. 0.0 .OR. A .LT. X) CALL XERMSG ('SLATEC', 'R9LGIT',
!     +   'X SHOULD BE GT 0.0 AND LE A', 2, 2)
      if ( x.le.0.0 .or. a.lt.x ) then
         write(*,*) "r9lgit:  Fatal Error:  "//
     x     "X should be GT 0.0 and LE A"
         stop
      endif
C
      AX = A + X
      A1X = AX + 1.0
      R = 0.0
      P = 1.0
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.0+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
!      CALL XERMSG ('SLATEC', 'R9LGIT',
!     +   'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 3, 2)
      write(*,*) "r9lgit: Fatal Error:  No convergence after 200"//
     x     "terms of continued fraction."
      stop
C
 30   HSTAR = 1.0 - X*S/A1X
!      IF (HSTAR .LT. SQEPS) CALL XERMSG ('SLATEC', 'R9LGIT',
!     +   'RESULT LESS THAN HALF PRECISION', 1, 1)
      if ( hstar.lt.sqeps ) write(*,*) "r9lgit:  Warning:  Result is"//
     x     "less than half precision"
C
      R9LGIT = -X - ALGAP1 - LOG(HSTAR)
C
      RETURN
      END
