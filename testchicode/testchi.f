      program testchi

c     Testchi provides an interactive interface to the code that calculates
c     the chi**2 of a given SFH model.  The user may change the parameters 
c     of the fit at runtime and see the effect of the changes on the chi**2
c
c     Jason Harris,  Revised May 24, 2000
c
      
      include '../sfhcode/sfh.h'

      character*40 datpre,maskfile,datfile(MCMD)
      character*40 ampfile
      double precision x(MP),y
      integer i,j,ntry

      real dat(MCMD,2,NDATA),chem(MP),age(MP)
      integer fstat,gtype
      integer nstars(MCMD),pnum(MP),mask(MCMD,MBOX)
      character*40 cmdfile,chifile
      common /cstat/ dat,cmdfile,chifile,nstars,mask,
     x               chem,age,pnum,fstat,gtype

      double precision psynth(MCMD,MP,MBOX)
      double precision pdat(MCMD,MBOX),pmod(MCMD,MBOX)
      common /cfit/ psynth, pmod, pdat

      character*8 suffix(MCMD)
      integer npix, ncmd, nbox(MCMD)
      real xlo(MCMD), xhi(MCMD), ylo(MCMD), yhi(MCMD), dpix
      common /cmds/ suffix, xlo, xhi, ylo, yhi, nbox, dpix, npix, ncmd

      integer niso,nheld,ihold(MP)
      double precision pheld(MP)
      common /hold/ pheld,ihold,niso,nheld

      integer np
      common /cnp/ np

      integer iverb
      common /verb/ iverb

c     *** Filenames ***
      call fetchchar40( 5, datpre )
      call fetchchar40( 5, cmdfile )
      call fetchchar40( 5, maskfile )
      call fetchchar40( 5, ampfile )
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
      call fetchint( 5, fstat )
      call fetchint( 5, gtype )
      call fetchint( 5, iverb )

      do icmd=1,ncmd
         nstars(icmd) = 0

         do ibox=1,nbox(icmd)
            pmod(icmd,ibox) = 0.0D0
            pdat(icmd,ibox) = 0.0D0
         end do
      end do
      
c     Set CMD grid mask
      call oldfile(44,maskfile)
 108  read(44,*,end=109) i,j,mask(i,j)
      goto 108
 109  close(44)

c     fitstat expects the hold arrays, but we don't need them 
c     for this program, so just fill null arrays.
      do i=1,niso
         ihold(i) = 0
         pheld(i) = 0.0
      end do
      nheld = 0

c     np is the number of independent, variable amplitudes
      np = niso 
      
c *** Read in input photometry:
c *** dat(icmd,1,i) = color; dat(icmd,2,i) = magnitude
      do icmd=1,ncmd
         call strcat(datpre, suffix(icmd), datfile(icmd))
         if (iverb.gt.0) write(*,*) 'Reading ', datfile(icmd)
         call oldfile(2,datfile(icmd))
         do i=1,NDATA
            read(2,*,end=10) dat(icmd,1,i),dat(icmd,2,i)
         end do
 10      close(2)
         nstars(icmd) = i - 1
      end do

c     Read in SFH amplitudes
      call oldfile(24,ampfile)
      do i=1,np
         read(24,*) rjunk, rjunk, x(i)
         if ( ihold(i).eq.1 ) then
            x(i) = pheld(i)
         endif
      end do
      close(24)

C     Compute the fitting statistic for the model:
      call fitstat( x, y, 1 )  ! 1 = write chifile

      if ( fstat.eq.0 ) then
         write(*,*) 'chi**2 of model: ', y
      else if ( fstat.eq.1 ) then
         write(*,*) 'Lorentzian stat. of model: ', y
      else if ( fstat.eq.2 ) then
         write(*,*) 'Poisson stat. of model: ', y
      endif

      stop
      end

