c     Utility functions for StarFISH
c
c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c

C     oldfile: opens an existing file.  Fails gracefully.
C
C     Jason Harris,  Rev. Feb 1, 2000
C
      subroutine oldfile(iu,filename)
      integer iu,iflag,ans
      character*(*) filename
      
 10   open(unit=iu,file=filename,status="old",iostat=iflag)

      if (iflag.ne.0) then
         write(*,*) 'Error: I could not find the file ', filename
         write(*,*) 'Exiting program.'
         stop
      endif

      return
      end

C     newfile: opens a new file.  Fails gracefully.
C
C     Jason Harris,  Rev. Feb 1, 2000
C
      subroutine newfile(iu,filename)
      integer iu,iflag,ans
      character*(*) filename

 10   open(unit=iu,file=filename,status="new",iostat=iflag)

      if (iflag.ne.0) then
         write(*,*) '*****WARNING*****'
         write(*,*) 'Error: I could not create the file ', filename
         write(*,*) 'Exiting program.'
         stop
      endif

      return
      end

C     unkfile: opens a new file/overwrites old file.  Fails gracefully.
C
C     Jason Harris,  Rev. Feb 1, 2000
C
      subroutine unkfile(iu,filename)
      integer iu,iflag,ans
      character*(*) filename

 10   open(unit=iu,file=filename,status="unknown",iostat=iflag)

      if (iflag.ne.0) then
         write(*,*) 'Error: I could not use the file ', filename
         write(*,*) 'Exiting program.'
         stop
      endif

      return
      end

C     strcat:  concatenate 2 strings, removing trailing whitespace from s1.
C
C     Jason Harris,  Rev. May 21, 2000
C

      subroutine strcat(s1,s2,s)

      integer i
      character*(*) s1, s2, s

      i = len(s1)
      do while (s1(i:i).eq.' ')
         i = i - 1
      end do

      s(1:i) = s1(1:i)
      i = i + 1
      s(i:(i+len(s2))) = s2(1:len(s2))

      do j=i+len(s2)+1,len(s)
         s(j:j) = ' '
      end do

      return 
      end

C     fetchreal: retrieve a real number from the next non-comment line
      subroutine fetchreal( nunit, f )
      integer nunit
      real f
      character*100 line

      call getnextline( nunit, line )
      call getreal( line, f )

      return
      end

C     fetchdouble: retrieve a real number from the next non-comment line
      subroutine fetchdouble( nunit, d )
      integer nunit
      real f
      double precision d
      character*100 line

      call getnextline( nunit, line )
      call getreal( line, f )
      d = dble(f)
      return
      end

C     fetchint: retrieve an integer from the next non-comment line
C     line on STDIN
      subroutine fetchint( nunit, i )
      integer nunit, i
      character*100 line

      call getnextline( nunit, line )
      call getint( line, i )

      return
      end

C     fetchchar40: retrieve a 40-character string from the next 
C     non-comment line
      subroutine fetchchar40( nunit, string )
      integer nunit
      character*40 string
      character*100 line

      call getnextline( nunit, line )
      call getword( line, string, 40 )

      return
      end

C     fetch2char40: retrieve two 40-character strings from the next 
C     non-comment line
      subroutine fetch2char40( nunit, string1, string2 )
      integer nunit
      character*40 string1, string2
      character*100 line
      logical line_ok

      call getnextline( nunit, line )

      line_ok = .false.
      do i=1,len(line)-1
         if ( line(i:i) .ne. ' ' ) then
            line_ok = .true.
            goto 10
         endif
      end do

 10   if ( line_ok .eqv. .false. ) then
         string1 = ''
         string2 = ''
         return
      endif

      call getword( line, string1, 40 )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 20
      end do
 20   line = line(i:len(line))

      call getword( line, string2, 40 )

      return
      end

C     fetch3char40: retrieve three 40-character strings from the next 
C     non-comment line
      subroutine fetch3char40( nunit, string1, string2, string3 )
      integer nunit
      character*40 string1, string2, string3
      character*100 line
      logical line_ok

      call getnextline( nunit, line )

      line_ok = .false.
      do i=1,len(line)-1
         if ( line(i:i) .ne. ' ' ) then
            line_ok = .true.
            goto 10
         endif
      end do

 10   if ( line_ok .eqv. .false. ) then
         string1 = ''
         string2 = ''
         string3 = ''
         return
      endif

      call getword( line, string1, 40 )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 20
      end do
 20   line = line(i:len(line))

      call getword( line, string2, 40 )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 30
      end do
 30   line = line(i:len(line))

      call getword( line, string3, 40 )

      return
      end

C     fetchchar8: retrieve an 8-character string from the next non-comment line
      subroutine fetchchar8( nunit, string )
      integer nunit
      character*8 string
      character*100 line

      call getnextline( nunit, line )
      call getword( line, string, 8 )

      return
      end

C     getnextline: retrieves the next non-comment line from
C     unit nunit
      subroutine getnextline( nunit, line )
      
      integer nunit
      character*100 line

 10   read(nunit,100,end=50) line
      if ( line(1:1) .eq. "#" .or. line(1:1) .eq. "*" ) goto 10
      return

      !End-of-file reached; return empty line
 50   line = ''

 100  format(a100)

      return
      end

C     readisodatline: parse a line from the file iso.dat into 4 values:
C     isochrone age, original iso file, lib file, and V_msto.
      subroutine readisodatline( iunit, age, infile, libfile, vmsto )

      integer i,iunit, iword
      double precision d
      real age, vmsto
      character*40 infile, libfile
      character*100 line
      logical line_ok

      call getnextline( iunit, line )

      line_ok = .false.
      do i=1,len(line)-1
         if ( line(i:i) .ne. ' ' ) then
            line_ok = .true.
            goto 10
         endif
      end do

 10   if ( line_ok .eqv. .false. ) then
         age = -1.0
         infile = ''
         libfile = ''
         vmsto = 99.9
         return
      endif

      call getreal( line, age )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 20
      end do
 20   line = line(i:len(line))

      call getword( line, infile, 40 )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 30
      end do
 30   line = line(i:len(line))

      call getword( line, libfile, 40 )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 40
      end do
 40   line = line(i:len(line))

      call getreal( line, vmsto )

      return
      end

C     readisolockline: parse a line from the iso.lock file into 3 values:
C     the index number, the lib filename stem,
C     and the output pixel file name stem.
      subroutine readisolockline( iunit, index, libfile, synthfile )

      integer i,iunit, index
      character*40 libfile, synthfile
      character*100 line
      logical line_ok

      call getnextline( iunit, line )

      line_ok = .false.
      do i=1,len(line)-1
         if ( line(i:i) .ne. ' ' ) then
            line_ok = .true.
            goto 10
         endif
      end do

 10   if ( line_ok .eqv. .false. ) then
         index = -1
         libfile   = ''
         synthfile = ''
         return
      endif

      call getint( line, index )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 20
      end do
 20   line = line(i:len(line))

      call getword( line, libfile, 40 )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 40
      end do
 40   line = line(i:len(line))

      call getword( line, synthfile, 40 )

      return
      end

C     readcmdline: parse a line from the cmdfile file into 3 
C     values:  the metallicity, the log(age), and the output
C     pixel file name stem.
      subroutine readcmdline( iunit, chem, age, pixbase )

      integer i,iunit
      double precision d
      real chem,age
      character*40 pixbase
      character*100 line
      logical line_ok

      call getnextline( iunit, line )

      line_ok = .false.
      do i=1,len(line)-1
         if ( line(i:i) .ne. ' ' ) then
            line_ok = .true.
            goto 10
         endif
      end do

 10   if ( line_ok .eqv. .false. ) then
         index = -1
         chem = -1.0
         age = -1.0
         pixbase = ''
         return
      endif

      call getreal( line, chem )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 20
      end do
 20   line = line(i:len(line))

      call getreal( line, age )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 30
      end do
 30   line = line(i:len(line))

      call getword( line, pixbase, 40 )

      return
      end

C     readtestpopline: parse a line from the cmdfile file into 4 
C     values:  the SFH amplitude, the metallicity, the log(age), 
C     and the number of isochrones in the isochrone group.
      subroutine readtestpopline( iunit, amp, chem, age, niso )

      integer i,iunit,niso
      double precision d
      real amp,chem,age
      character*100 line
      logical line_ok

      call getnextline( iunit, line )

      line_ok = .false.
      do i=1,len(line)-1
         if ( line(i:i) .ne. ' ' ) then
            line_ok = .true.
            goto 10
         endif
      end do

 10   if ( line_ok .eqv. .false. ) then
         amp = -1.0
         chem = -1.0
         age = -1.0
         niso = -1
         return
      endif

      call getreal( line, amp )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 20
      end do
 20   line = line(i:len(line))

      call getreal( line, chem )
      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 30
      end do
 30   line = line(i:len(line))

      call getreal( line, age )

      ! remove the word we just read from the line string
      ! (including any leading whitespace)
      iword = 0
      do i=1,len(line)
         if ( line(i:i).ne.' ' ) iword = 1
         if ( iword.eq.1 .and. line(i:i).eq.' ' ) goto 40
      end do
 40   line = line(i:len(line))

      call getint( line, niso )

      return
      end
