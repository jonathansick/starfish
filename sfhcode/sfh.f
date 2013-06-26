      program sfh

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
C     sfh: Given observed photometry, and a library of synthetic CMDs, 
C     find the set of SFR amplitudes which produce model CMDs that 
C     best match the data.
C
C     The program is flexible; you have the choice of three fitting 
C     statistics: chi-squared, Lorentzian, or Poisson.  
C     
C     Once the best-fit parameters are determined, SFH then determines the 
C     uncertainty on each parameter.  It can do this either by searching 
C     parameter space around the best-fit location, or by holding each 
C     parameter fixed in turn, at several different values, until the 
C     range of values which provide an "acceptable" fit are found. 
C
C     Additional features: You may hold certain amplitudes fixed at 
C     specified values.  You may define the CMDs that are used in the fits, 
C     and exactly how each is gridded.  You may also mask out certain 
C     regions of the CMDs.
C
      include 'sfh.h'

      character*40 datpre,holdfile,maskfile,datfile(MCMD)
      character*40 outfile,logfile,logfile2
      double precision p(MP+1,MP),ppos(MP),pneg(MP)
      double precision pmin(MP),pold(MP),dp
      double precision y(MP+1),ymin,ytry,yerr,yold
      double precision ftol0,ftol
      integer i,j,itry
      integer ntry,uselog,nvec
      real cthr,ptol,conf
      real lamda, lamda0

      real dat(MCMD,2,NDATA),chem(MP),age(MP)
      integer fstat,gtype
      integer nstars(MCMD),pnum(MP),mask(MCMD,MBOX)
      character*40 cmdfile,chifile
      common /cstat/ dat,cmdfile,chifile,nstars,mask,
     x               chem,age,pnum,fstat,gtype

      double precision psynth(MCMD,MP,MBOX)
      double precision pdat(MCMD,MBOX),pmod(MCMD,MBOX)
      common /cfit/ psynth, pmod, pdat

      integer niso,nheld,ihold(MP)
      double precision pheld(MP)
      common /hold/ pheld,ihold,niso,nheld

      integer np
      common /cnp/ np

      integer iverb
      common /verb/ iverb

      character*8 suffix(MCMD)
      integer npix, ncmd, nbox(MCMD)
      real xlo(MCMD), xhi(MCMD), ylo(MCMD), yhi(MCMD), dpix
      common /cmds/ suffix, xlo, xhi, ylo, yhi, nbox, dpix, npix, ncmd
      logical done
c
c *** Read input parameters from STDIN
c     lines beginning with "#" or "*" are skipped
c     only the text up to the first space (" ") char on each line is parsed

c     *** Filenames ***
      call fetchchar40( 5, datpre )
      call fetchchar40( 5, cmdfile )
      call fetchchar40( 5, maskfile )
      call fetchchar40( 5, holdfile )
      call fetchchar40( 5, outfile )
      call fetchchar40( 5, logfile )
      call fetchchar40( 5, logfile2 )
      call fetchchar40( 5, chifile )

c     *** Synthetic CMD parameters ***
      call fetchint( 5, niso )
      call fetchint( 5, ncmd )
      call fetchint( 5, npix )
      call fetchreal( 5, dpix )
      do icmd=1,ncmd
         call fetchchar8( 5, suffix(icmd) )
         call fetchreal( 5, xlo(icmd) )
         call fetchreal( 5, xhi(icmd) )
         call fetchreal( 5, ylo(icmd) )
         call fetchreal( 5, yhi(icmd) )
         call fetchint( 5, nbox(icmd) )
      end do

c     *** Runtime parameters ***
      call fetchint( 5, iseed )
      call fetchint( 5, fstat )
      call fetchint( 5, uselog )
      call fetchint( 5, iplg )
      call fetchint( 5, gtype )
      call fetchint( 5, iverb )
      call fetchreal( 5, lamda0 )
      call fetchreal( 5, conf )
      call fetchreal( 5, cthr )
      call fetchreal( 5, ptol )
      call fetchdouble( 5, ftol0 )
      call fetchint( 5, nvec )
      call fetchint( 5, ntry )

c
c *** Initialize variables
c
      lamda = lamda0
      y(1)  = 1.0D8
      ftol  = 1.0

      do i=1,niso
         ppos(i) = -1.0D30
         pneg(i) =  1.0D30
      end do

      do icmd=1,ncmd
         nstars(icmd) = 0

         do ibox=1,nbox(icmd)
            pmod(icmd,ibox) = 0.0D0
            pdat(icmd,ibox) = 0.0D0
         end do
      end do

c     Initialize random number generator.
      call sdprnd(iseed)

c     Set CMD grid mask
      call oldfile(44,maskfile)
 108  read(44,*,end=109) i,j,mask(i,j)
      goto 108
 109  close(44)

c     Read in fixed amplitudes, if any
c     Don't use oldfile() because we just want to fail silently,
c     if no holdfile
      open(45, file=holdfile, iostat=iflag)
      j = 1
      if ( iflag .eq. 0 ) then
         do i=1,niso
            read(45,*) value
            if ( value.gt.0.0 ) then
               ihold(i) = j
               pheld(j) = value
               j = j + 1
            else 
               ihold(i) = 0
            endif
         end do
      else
         do i=1,niso
            ihold(i) = 0
         end do
      endif

c     np is the number of independent, variable amplitudes
      nheld = j - 1
      np = niso - nheld

      if ( iverb.ge.3 ) 
     x     write(*,*) "Number of free isochrone groups: ", np

c     
c *** Read in input photometry:
c *** dat(icmd,1,i) = color; dat(icmd,2,i) = magnitude
c
      do icmd=1,ncmd
         call strcat(datpre, suffix(icmd), datfile(icmd))
         if (iverb.ge.4) write(*,*) 'Reading ', datfile(icmd)
         call oldfile(2,datfile(icmd))
         do i=1,NDATA
            read(2,*,end=10) dat(icmd,1,i),dat(icmd,2,i)
         end do
 10      close(2)
         nstars(icmd) = i - 1
      end do

c     Open the position archive file
      if (iplg.eq.1) call unkfile(3,logfile2)

c *** 
c *** Begin Convergence.
c ***
      if (iverb.ge.1) write(*,*) 'Begin convergence. '

      ! Find the global minimum with a modified downhill simplex
      ! "amoeba" algorithm.  Once the amoeba finds a minimum, 
      ! its size is reset, and it is allowed to reconverge.  
      ! If it reconverges to the same location, we then search 
      ! the local parameter space along random directions for 
      ! lower chi**2 values.  If lower values are found, 
      ! the amoeba is restarted at the new lowest location.  
      ! Otherwise, the minimum is accepted.

      if (iverb.ge.3) write(*,*) ' Initializing simplex '
         
c     If logfile exists and uselog=1, start at logged position...
      open(unit=11,file=logfile,type="old",iostat=iflag)
      if (iflag.eq.0.and.uselog.eq.1) then
         do jj=1,np
            read(11,2) pmin(jj)
         end do
         read(11,3) lamda, ymin

c     ...otherwise, start at random position (and improve on it 
c     with a random-search loop)
      else
         do i=1,np
            call random(pmin(i),0.D0,100000.D0)
         end do

         ! Begin by searching in nvec random directions to find a low 
         ! initial chi**2 value (this step is a big time saver; the amoeba 
         ! takes much longer to reach similarly low chi**2 values).
         if ( iverb.ge.2 ) write(*,*) ' Improving initial position... '
         call search(pmin,y(1),nvec,lamda)
         call logp( logfile, pmin, y(1), lamda )
      endif
      close(11)

      ! Run the amoeba to find the location in parameter space where  
      ! y (the value of the model fitting statistic) is minimized.  
 13   if (iverb.ge.1) write(*,*) ' Starting the downhill simplex... '
      call simplex(pmin,p,y, lamda)
      call amoeba(p,y,iter,ftol)

      ! record the minimum chi**2 value (ymin) and the 
      ! location of the minimum (pmin(i)).  
      ymin = y(1)
      do j=1,np
         pmin(j) = p(1,j)
      end do
      call logp( logfile, pmin, ymin, lamda )

      ! In order to make sure we are not at a local minimum,
      ! we "nudge" the amoeba by re-expanding it to its original
      ! size and letting it shrink and reconverge.  If it converges
      ! again to the same pmin location (within +/- ptol), 
      ! then we accept the minimum.  The amoeba is re-expanded by
      ! calling simplex.
 15   if (iverb.ge.3) write(*,*) '   Nudging amoeba...'
      call simplex(pmin,p,y,lamda)
      call amoeba(p,y,iter,ftol)
      do j=1,np
         pmin(j) = p(1,j)
      end do
      call logp( logfile, pmin, y(1), lamda )

      if (iverb.ge.2) write(*,*) ' :)  ymin, ftol: ', y(1), ftol

      ! Check for improvement over the last attempt.
      inew = 0
      do j=1,np
         if (abs(p(1,j)-pmin(j)).gt.ptol) then
            inew = inew + 1
         endif
      end do
      
      if (inew.gt.0.and.y(1).lt.ymin) then
         ! we have found a better minimum.  Update pmin/ymin 
         ! and "nudge" again.
         if (iverb.ge.2) write(*,*) ' Current minimum: ', y(1)
         ymin = y(1)
         do j=1,np
            pmin(j) = p(1,j)
         end do
         call logp(logfile,pmin,ymin,lamda)
         goto 15  ! start next nudge
      endif
      
      ! Nudging did not improve the minimum.  
      ! Try reducing ftol and nudge again.

      if (ftol.gt.(ftol0*5.0)) then
         ftol = ftol/10.0
         goto 15  ! start next nudge
      endif
      
      ! ftol has reached its minimum value.
      ! Do one more nudge.
      if (ymin.lt.y(1)) then
         !Record the new global minimum candidate.
         y(1) = ymin
         do j=1,np
            p(1,j) = pmin(j)
         end do 
         call logp(logfile,pmin,ymin,lamda)
      
         goto 15
      end if
      
      ! OK, that's the best we could do with amoeba.
      ! Search local parameter space along random directions 
      ! for lower chi**2 values.
      ! If a smaller chi**2 is found, we loop back and start
      ! the amoeba again at the new location.
      if (iverb.ge.1) write(*,*) ' Starting random-search loop...'
      do itry=1,ntry
         if (iverb.ge.3) write(*,*) 'Round ', itry
         ytry = y(1)
         do j=1,np
            pold(j) = p(1,j)
         end do
         
         call search(p,ytry,nvec,lamda0)
         
         if (ytry.lt.(y(1)-cthr)) then
            y(1) = ytry
            do j=1,np
               pmin(j) = p(1,j)
            end do
            call logp(logfile,pmin,ymin,lamda)
            
            if (iwrt.ge.2) 
     x           write(*,*) ' random search found new minimum: ', y(1)
            goto 13  ! back to the amoeba!
         else
            do j=1,np
               p(1,j) = pold(j)
            end do
         endif
      end do
      
      !That's it; we have iterated through a full cycle of 
      !the amoeba and the random-search without finding a 
      !lower minimum.  Accept the global minimum.
      do j=1,np
         pmin(j) = p(1,j)
      end do
      call logp(logfile,pmin,ymin,lamda)

      if ( iverb.ge.1 ) then
         write(*,*) 'Convergence completed.  Final minimum: ', ymin
      endif
 
c *** Done with minimization, compute errors and write output files

      call fitstat(pmin,ymin,1) ! <--- 1 = write chifile

      if (iverb.ge.1) write(*,*) 'Determining fit uncertainty: ' 
      yerr = ymin + dchi(np,conf)
      call errs(pmin,ymin,yerr,nvec,lamda0,ppos,pneg)

      if (iverb.ge.3) write(*,*) 'Writing output amplitudes: '
      call newfile(4,outfile)

      jj = 0
      jfree = 1
      do j=1,niso
         if ( ihold(j) .eq. 0 ) then  ! free amplitude
            do k=1,pnum(j)
               jj = jj + 1
               write(4,6) chem(jj),age(jj),p(1,jfree),
     x              pneg(jfree),ppos(jfree)
               jfree = jfree + 1
            end do
         else  ! fixed amplitude
            do k=1,pnum(j)
               jj = jj + 1
               write(4,6) chem(jj),age(jj),(pheld(ihold(j)),ii=1,3)
            end do
         endif
      end do
      
      if (iverb.ge.1) write(*,*) ' '
      if (iverb.ge.1) write(*,*) 'Final minimum: ', y(1)

      close(3)
      close(4)

 2    format(g14.7)
 3    format(2g12.6)
 6    format(f7.4,f7.2,3(2x,g15.8))

      stop
      end

