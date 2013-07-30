      program interp
c     Use the empirical EEPs in the isochrone files as anchor points
c     from which to interpolate between bracketing isochrone pairs.
c
c     First must make sure the set of EEPs match between the two isochrones.
c     In geeneral, there are slight differences in the EEP sequences.  Must
c     come up with a clever way to find corresponding EEP pairs.
c     
c     Some details: the EEP is an integer flag, that is non-zero if
c     the point represents a "node" in the B-V,V CMD.  A node is simply
c     a local extremum in either magnitude or color.  The value of the
c     EEP flag indicates the type of node: 
c       0 = not a node
c       1 = color minimum
c       2 = color maximum
c       3 = mag minimum
c       4 = mag maximum

c     We can attempt to match-up corresponding EEPs in the pair of 
c     isochrones by matching the flag numbers.  However, this isn't 
c     foolproof; sometimes the same node manifests as different
c     kinds of extrema.  As a second check, we will use the offset
c     in mass between candidate EEP matches; the mass offset remains 
c     somewhat constant for all pairs.

      parameter(MMAG=8,NN=400,SS=1)
      character*40 afile,bfile,outfile
      character*16 fstr(MMAG)
      real mass(3,NN),mag(3,MMAG,NN)
      real frac,fa,fb
      real amag(MMAG),bmag(MMAG),amass,bmass
      integer i,i0,nn1,nn2,niso1,niso2,nmag
      integer n1,n2,n3,o1,o2,p1,p2
      integer node1(NN),node2(NN),inode1(NN),inode2(NN),node3

      fstr(1) = '(f13.8,f8.3,i4)'
      fstr(2) = '(f13.8,2f8.3,i4)'
      fstr(3) = '(f13.8,3f8.3,i4)'
      fstr(4) = '(f13.8,4f8.3,i4)'
      fstr(5) = '(f13.8,5f8.3,i4)'
      fstr(6) = '(f13.8,6f8.3,i4)'
      fstr(7) = '(f13.8,7f8.3,i4)'
      fstr(8) = '(f13.8,8f8.3,i4)'

      call fetchint( 5, nmag )
      call fetchint( 5, imag )
      call fetchreal( 5, frac )
 10   call fetch3char40( 5, afile, bfile, outfile )
      if ( afile.eq.'' ) goto 100  ! end-of-file

c     Read in the bracketing isochrones
      i = 1
      i0 = 1
      open(unit=11,file=afile,status="old")
 11   read(11,*,end=12) mass(1,i),(mag(1,j,i),j=1,nmag),node1(i)
      if ( node1(i).gt.0 ) then
         inode1(i0) = i
         i0 = i0 + 1
      endif
      i = i + 1
      goto 11
 12   close(11)
      niso1 = i - 1
      nn1 = i0 - 1

      i = 1
      i0 = 1
      open(unit=12,file=bfile,status="old")
 13   read(12,*,end=14) mass(2,i),(mag(2,j,i),j=1,nmag),node2(i)
      if ( node2(i).gt.0 ) then
         inode2(i0) = i
         i0 = i0 + 1
      endif
      i = i + 1
      goto 13
 14   close(12)
      niso2 = i - 1
      nn2 = i0 - 1

      open(unit=13,file=outfile,status="unknown")

      o1 = 0
      o2 = 0
      p1 = 0
      p2 = 0
      j1 = 1
      j2 = 1

c     20 is the beginning of the main loop through the EEPs ("nodes")
c     in the bracketing isochrones (iso1 and iso2).
c     We attempt to match up nodes in the bracketing isochrones.
c     To be considered a match, a pair of nodes should have the same
c     flag value, and should also be near each other in the CMD.
c     inode1(j1) is the isochrone index for node j1 in iso1
c     inode2(j2) is the isochrone index for node j2 in iso2
c     node1 and node2 are the node values (1-4) for nodes in iso1/iso2
c
c     There are several criteria for finding a match, described below.
c
c     If we are at the end of an isochrone, do one last interpolation
c     (goto 26)
 20   if ( j1.gt.nn1 .or. j2.gt.nn2 ) then
         o1 = p1
         o2 = p2
         p1 = niso1
         p2 = niso2
         goto 26
      endif

c     The current nodes in iso1 and iso2 are identified by j1 and j2
c     Do these nodes have the same value?  If so, see if their position
c     is within 0.5 mag of a perfect match.
c     (The match radius calculation takes into account the possibility 
c     that the previous matched pair were not at the same spot, so 
c     we should expect a similar offset in the next pair)
      if ( node1(inode1(j1)).eq.node2(inode2(j2)) ) then
         v1 = mag(1,imag,inode1(j1))
         v2 = mag(2,imag,inode2(j2))
         c1 = mag(1,imag-1,inode1(j1)) - mag(1,imag,inode1(j1))
         c2 = mag(2,imag-1,inode2(j2)) - mag(2,imag,inode2(j2))
         if ( o1.ne.0 .and. o2.ne.0 ) then
            vo1 = mag(1,imag,o1)
            vo2 = mag(2,imag,o2)
            co1 = mag(1,imag-1,o1) - mag(1,imag,o1)
            co2 = mag(2,imag-1,o2) - mag(2,imag,o2)
            r = sqrt( SS*(c2-c1 + co1-co2)**2 + (v2-v1+vo1-vo2)**2)
         else
            r = sqrt( SS*(c2-c1)**2 + (v2-v1)**2)
         endif

c        If the match is good, accept it by going to 25
c        otherwise, we'll keep looking...
         if ( r .lt. 0.5 ) then
            goto 25
         endif
      endif 

c     No match yet.  We will next try incrementing to the next node 
c     in either iso1 or iso2.  However, it is also possible that 
c     both increments will produce a value-match.  If so, adopt the 
c     best match between the two.
c
c     Note that we always return to the top of the loop (20) if we find a
c     possible match; the only way to get to 25 is in the above block.
      if ( node1(inode1(j1)).eq.node2(inode2(j2+1)) .and. 
     x     node1(inode1(j1+1)).eq.node2(inode2(j2)) ) then
         
         v1 = mag(1,imag,inode1(j1))
         v2 = mag(2,imag,inode2(j2+1))
         c1 = mag(1,imag-1,inode1(j1)) - mag(1,imag,inode1(j1))
         c2 = mag(2,imag-1,inode2(j2+1)) - mag(2,imag,inode2(j2+1))
         if ( o1.ne.0 .and. o2.ne.0 ) then
            vo1 = mag(1,imag,o1)
            vo2 = mag(2,imag,o2)
            co1 = mag(1,imag-1,o1) - mag(1,imag,o1)
            co2 = mag(2,imag-1,o2) - mag(2,imag,o2)
            r1 = sqrt( SS*(c2-c1 + co1-co2)**2 + (v2-v1 + vo1-vo2)**2)
         else
            r1 = sqrt( SS*(c2-c1)**2 + (v2-v1)**2)
         endif

         v1 = mag(1,imag,inode1(j1+1))
         v2 = mag(2,imag,inode2(j2))
         c1 = mag(1,imag-1,inode1(j1+1)) - mag(1,imag,inode1(j1+1))
         c2 = mag(2,imag-1,inode2(j2)) - mag(2,imag,inode2(j2))
         if ( o1.ne.0 .and. o2.ne.0 ) then
            vo1 = mag(1,imag,o1)
            vo2 = mag(2,imag,o2)
            co1 = mag(1,imag-1,o1) - mag(1,imag,o1)
            co2 = mag(2,imag-1,o2) - mag(2,imag,o2)
            r2 = sqrt( SS*(c2-c1 + co1-co2)**2 + (v2-v1 + vo1-vo2)**2)
         else
            r2 = sqrt( SS*(c2-c1)**2 + SS*(v2-v1)**2)
         endif

c        Increment whichever index (j1 or j2) produces the best match
c        and return to 20 to evaluate this new pair.
         if ( r1 < r2 ) then
            inode2(j2) = -1
            j2 = j2 + 1
         else 
            inode1(j1) = -1
            j1 = j1 + 1
         endif

         goto 20
      endif

c     Try skipping ahead in list 1.  If that makes a value match, 
c     increment j1 and return to 20
      if ( node1(inode1(j1+1)).eq.node2(inode2(j2)) ) then
c         write(*,*) "skipping inode1(", j1, ")"
         inode1(j1) = -1
         j1 = j1 + 1
         goto 20
      endif

c     Try skipping ahead in list 2.  If that makes a value match, 
c     increment j2 and return to 20
      if ( node1(inode1(j1)).eq.node2(inode2(j2+1)) ) then
c         write(*,*) "skipping inode2(", j2, ")"
         inode2(j2) = -1
         j2 = j2 + 1
         goto 20
      endif

c     Ok, we could not find a value-match between j1 and j2,
c     nor between j1 and (j2+1) nor (j1+1) and j2.
c     It's possible that j1 and j2 are really a match, but
c     have different values for some reason.  Because we are 
c     giving up the requirement for a value-match, we require 
c     the photometric match radius to be smaller here.
      v1 = mag(1,imag,inode1(j1))
      v2 = mag(2,imag,inode2(j2))
      c1 = mag(1,imag-1,inode1(j1)) - mag(1,imag,inode1(j1))
      c2 = mag(2,imag-1,inode2(j2)) - mag(2,imag,inode2(j2))
      if ( o1.ne.0 .and. o2.ne.0 ) then
         vo1 = mag(1,imag,o1)
         vo2 = mag(2,imag,o2)
         co1 = mag(1,imag-1,o1) - mag(1,imag,o1)
         co2 = mag(2,imag-1,o2) - mag(2,imag,o2)
         r = sqrt( SS*(c2-c1 + co1-co2)**2 + (v2-v1 + vo1-vo2)**2)
      else
         r = sqrt( SS*(c2-c1)**2 + (v2-v1)**2)
      endif
      if ( r .lt. 0.25 ) then
         goto 25
      endif

c     I give up!  Skip the current pair of nodes.
      inode1(j1) = -1
      inode2(j2) = -1
      j1 = j1 + 1
      j2 = j2 + 1
      goto 20


c     25: we have now identified the next matched pair of nodes.
c     We will now interpolate between the bracketing isochrones in 
c     the interval between the last node and the current node.

c     Two-stage interpolation.  We identify n3 points along the interval
c     for which we will interpolate photometry.  n3 = 0.5(n1+n2), where
c     n1 is the number of points along iso1's interval and n2 is the same 
c     for iso2.  

c     We first interpolate along each bracketing isochrone to obtain n3
c     points which are uniformly distributed along the interval.
c     We then interpolate between these pairs of points to get the 
c     interpolated photometry
 25   o1 = p1
      o2 = p2
      p1 = inode1(j1) 
      p2 = inode2(j2) 
 26   n1 = p1 - o1
      n2 = p2 - o2
      n3 = int(0.5*(n1 + n2))

c     Uncomment to see data on matched nodes
c      write(*,*) node1(p1), node2(p2), r, mass(1,p1)-mass(2,p2)

      do i=1,n3
c        fa is distance along isoc1 interval, runs from 0 to n1
c        later, we will subtract ia so that fa is the fractional 
c        distance between isoc1 points (ia) and (ia+/-1).
         fa = real(o1) + real(i)*real(n1)/real(n3) 
         ia = int(fa)           

         if (ia.ge.p1) then
            ia = p1 
            fa = fa - ia 
            do k=1,nmag
               amag(k)=mag(1,k,ia)+fa*(mag(1,k,ia)-mag(1,k,ia-1))
            end do
            amass = mass(1,ia) + fa*(mass(1,ia)-mass(1,ia-1))
         else if (ia.lt.1) then
            ia = 1
            fa = fa - ia
            do k=1,nmag
               amag(k)=mag(1,k,ia)+fa*(mag(1,k,ia+1)-mag(1,k,ia))
            end do
            amass = mass(1,ia) + fa*(mass(1,ia+1)-mass(1,ia))
         else
            fa = fa - ia 
            do k=1,nmag
               amag(k) = (1-fa)*mag(1,k,ia) + fa*mag(1,k,ia+1)
            end do
            amass = (1-fa)*mass(1,ia) + fa*mass(1,ia+1)
         endif
            
         fb = real(o2) + real(i)*real(n2)/real(n3) 
         ib = int(fb)

         if (ib.ge.p2) then
            ib = p2 
            fb = fb - ib
            do k=1,nmag
               bmag(k)=mag(2,k,ib)+fb*(mag(2,k,ib)-mag(2,k,ib-1))
            end do 

            bmass = mass(2,ib) + fb*(mass(2,ib)-mass(2,ib-1))
         else if (ib.lt.1) then
            ib = 1
            fb = fb - ib
            do k=1,nmag
               bmag(k)=mag(2,k,ib)+fb*(mag(2,k,ib+1)-mag(2,k,ib))
            end do
            bmass = mass(2,ib) + fb*(mass(2,ib+1)-mass(2,ib))
         else
            fb = fb - ib
            do k=1,nmag
               bmag(k) = (1-fb)*mag(2,k,ib) + fb*mag(2,k,ib+1)
            end do

            bmass = (1-fb)*mass(2,ib) + fb*mass(2,ib+1)
         endif                     

         do k=1,nmag
            mag(3,k,i) = (1.0-frac)*amag(k) + frac*bmag(k)
         end do
         mass(3,i) = (1.0-frac)*amass + frac*bmass

c         write(*,*) amag(imag),bmag(imag),mag(3,imag,i),frac

         node3 = 0
         if ( i.eq.n3 ) node3 = node1(p1)
         write(13,fstr(nmag)) 
     x        mass(3,i),(mag(3,k,i),k=1,nmag),node3

      end do

c     Are we all finished?
      if ( p1.eq.niso1 .or. p2.eq.niso2 ) goto 99

c     Done with this interval; try to find the next matched node pair
      j1 = j1 + 1
      j2 = j2 + 1
      goto 20




 99   close(13)
      goto 10

 100  stop
      end
