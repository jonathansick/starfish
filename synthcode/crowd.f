      subroutine crowd(x,y,dropflag)

c     This file is part of StarFISH
c     (C) 2004 by Jason Harris  <jharris@as.arizona.edu>
c     
c -=> StarFISH is free software; you can redistribute it 
c -=> and/or modify it under the terms of the GNU General Public 
c -=> License (GPL) as published by the Free Software Foundation;
c -=> either version 2 of the License, or (at your option) any 
c -=> later version.
c
c     Apply photometric error offsets to current model star.  Also, 
c     determine whether the star is detected or not (to simulate data's 
c     incompleteness).  Non-detected stars are called "dropped" stars.
c
c     The photometric errors are drawn from an empirical  photometric 
c     error model, based on the results of artificial star tests.
c
c     The empirical modeling of the photometric errors is a bit 
c     complicated, and deserves further explanation:  

c     For simplicity, consider a 1-dimensional distribution of photometric 
c     errors.  In an ideal world, this distribution would be centered at 
c     zero (i.e., no systematic errors) and would have a Gaussian profile 
c     with some small width.  Unfortunately, real photometric errors are 
c     not usually so nice; they have asymmetric non-Gaussian tails, and 
c     often peak at some nonzero value.  Rather than try to impose some 
c     kind of complex analytic model on such an error distribution, I have 
c     opted to empirically represent the actual distribution as closely as 
c     I can.  I do this by binning the distribution along the delta-mag 
c     axis, and storing the delta-mag histogram.
c
c     However, this needs to be done separately for different regions of 
c     the CMD, because the photometric error characteristics depend on 
c     both brightness and color.  Therefore, we divide the CMD(s) into 
c     "crowding bins", such that the artificial stars in each bin can be
c     considered to have the same photometric error properties.  In 
c     addition, to avoid "quantum jumps" in the error characteristics
c     at crowding-bin boundaries, we will need to interpolate between
c     the delta-mag distributions in the bracketing crowding bins
c     (more on this below).
c
c     Furthermore, there are photometric errors in both the X (color)
c     and Y (brightness) directions of the CMDs.  In previous versions 
c     of StarFISH, we treated these independently, using two 
c     one-dimensional error histograms for each crowding bin.  However, 
c     this makes it impossible to reproduce correlated (dX,dY) errors
c     in the model.  We are now using a two-dimensional binning of the
c     errors; in other words, a "pixelization" of the (dX,dY) parameter 
c     space.
c
c     Out interpolation strategy requires us to represent the (dX,dY)
c     pixelization as a 1-dimensional cumulative function.  To 
c     understand this, it helps to consider the 1-dimensional case 
c     first.  Consider two histograms, between which you need to 
c     interpolate to find a third histogram.  One way to do so is to 
c     construct the *cumulative* distribution function for each 
c     histogram.  This makes it trivial to interpolate, because the 
c     interpolated cumulative function always lies between the bracketing
c     cumulative functions.  For example, if you wanted a histogram 
c     that was X% of the way between histograms A and B, you first convert
c     A and B to cumulative functions Ca and Cb.  Each cumulative 
c     function covers some range in delta-mag, but the vertical range 
c     is always between 0 and 1 (because we normalize them to the total 
c     number of stars in each histogram).  The interpolated cumulative 
c     function is then just the set of points which are X% of the way 
c     from Ca to Cb, measured horizontally (along the delta-mag axis).
c
c     That's the picture for one-dimensional histograms, but it can be 
c     generalized to our two-dimensional "pixelizations".  The 
c     cumulative function for a pixel-array (i.e., two-dimensional 
c     histogram) is a one-dimensional vector; the cumulative pixel value 
c     sum in a nested loop spanning both rows and columns in the 
c     pixel-array.  For example, consider the pixel array:
c
c        0   1   1   0
c        1   2   2   1
c        1   2   2   1
c        0   1   1   0
c
c     The cumulative vector starts out by just accumulating across 
c     the top row:  ( 0, 1, 2, 2, ... ); and then it continues with the 
c     next row:  ( ..., 3, 5, 7, 8, ... ); and finally the last two rows: 
c     ( ..., 9, 11, 13, 14, 14, 15, 16, 16 ).  
c     
c     However, there's an additional wrinkle introduced by the two-D 
c     generalization.  We need to be able to use the cumulative function 
c     to select a position in the (dX,dY) space stochastically.  We do this
c     by choosing a random number between 0 and 1, which identifies a
c     vertical position on the (normalized) cumulative function.  We 
c     then determine the delta-mag position where the cumulative function 
c     intersects that vertical position.  This point will generally lie 
c     between two of the actual values of the cumulative function vector.
c     In the 1-D case, this is no problem; you simply interpolate between
c     the bracketing cumfcn values to find the right delta-mag value.  
c     However, in 2-D, we can only interpolate along rows, because adjacent 
c     entries in the cumfcn vector are (generally) adjacent along rows.  
c     To avoid striping artifacts in the final (dX,dY) distribution, we 
c     will scatter each interpolated point vertically by +/- half of the 
c     pixel size.
c
c     So, givan all that, the procedure for crowd is the following:
c     We are given a star whose position in each cmd is (x(ncmd),y(ncmd)).
c     For each CMD, we identify the crowding bins which bracket its 
c     position.  Next, we interpolate linearly between the dropout 
c     fractions in each bracketing crowding bin, to the position (x,y).
c     This interpolated dropout fraction gives the probability that the 
c     input star is undetected in the current CMD.  We draw a random number 
c     between 0 and 1; if it is less than the dropout probability, then 
c     the star is flagged as a dropout, and the subroutine exits.  If the 
c     star is not a dropout, we continue with the photometric error modeling.
c     We select another random number between 0 and 1, and use this number 
c     to select (dX,dY) pairs from the cumulative function of each of the 
c     four bracketing crowding bins (as described above).  We then 
c     interpolate linearly between these (dX,dY) pairs, based on the 
c     position (x,y) relative to the bracketing bin centers.  Finally, 
c     the interpolated (dx0,dy0) offsets are applied to the input 
c     (x,y) position.
c
      include 'synth.h'

      integer iflag,dropflag(MCMD)
      integer ix,iy
      real x(MCMD),y(MCMD),dx,dy
      real dm(2,2,2),xbinmag,ybinmag, rr
      double precision f, dr

      integer nmag,ncmd
      common /c_nmag/ nmag,ncmd

      integer nx(MCMD),ny(MCMD)
      real dbinx,dbiny,xmin(MCMD),ymin(MCMD)
      common /c_crowd/ dbinx,dbiny,nx,ny,xmin,ymin

      real pdrop(MCMD,BINX,BINY)
      real cumerr(MCMD,BINX,BINY,NBINS*NBINS)
      common /c_cum/ pdrop,cumerr

      save iflag
      
C     If we haven't yet read the crowding lookup table, read it now.
C     If it doesn't exist yet, readcrowd will create it.
      if (iflag.eq.0) then
         iflag = 1
         call readcrowd()
      endif

C     Loop over each CMD
      do icmd=1,ncmd
         ! Initialize the dropflag variable (will be set to 1 if 
         ! the star drops out)
         dropflag(icmd) = 0
         
         ! ix,iy identify the crowding bin containing 
         ! the input point (x,y)
         ix = int((x(icmd) - xmin(icmd))/dbinx) + 1
         iy = int((y(icmd) - ymin(icmd))/dbiny) + 1

         ! Points that lie outside the bounds of the 
         ! crowding bins are dropped
         if ( ix.lt.1.or.ix.gt.nx(icmd) .or.
     x        iy.lt.1.or.iy.gt.ny(icmd) ) then
            dropflag(icmd) = 1
            dx = 0.0
            dy = 0.0
            goto 200 ! skip phot errs stuff
         endif

         if ( interp_errs.eq.1 ) then
            ! pick bracketing crowding bins for bin-interpolation
            jx = ix + 1
            jy = iy + 1
         
            ! If jx|jy is beyond nx|ny, 
            ! we need to use (ix-1)|(iy-1) instead:
            if ( jx.gt.nx(icmd) ) then
               jx = nx(icmd)
               ix = jx - 1
            endif
            if ( jy.gt.ny(icmd) ) then
               jy = ny(icmd)
               iy = jy - 1
            endif

            ! corner position of crowding bin (ix,iy)
            xbinmag = (ix-1)*dbinx + xmin(icmd) 
            ybinmag = (iy-1)*dbiny + ymin(icmd) 
            ! input position relative to crowding bin (ix,iy) 
            dxbin = (x(icmd)-xbinmag)/dbinx 
            dybin = (y(icmd)-ybinmag)/dbiny 

            ! Select delta-mags from each of the four bracketing 
            ! crowding bins using a single random number to index 
            ! the cumulative error distribution in each bin.
            ! The four pairs of dm values (dx,dy) are returned in dm.
            call bin_dmags( icmd, ix, iy, dm )

            ! Interpolate the dm values and dropout rates from the 
            ! bin centers to the input position (x,y).
            ! The interpolated values are returned through (dx,dy) and p0.
            call bin_interp( icmd, ix, iy, dxbin, dybin, 
     x                       dm, pdrop, dx, dy, p0 )
         else
            ! Choose a random number for indexing the cumulative 
            ! error functions.  Because we are recasting double precision
            ! to real, it is possible to have a value that is equal
            ! to 1.0 (the original double was less than 1.0).  When
            ! this occurs, draw a different random number.
 199        call random( dr, 0.0D0, 1.0D0 )
            rr = real(dr)
            if ( rr.eq.1.0 ) goto 199

            p0 = pdrop(icmd,ix,iy) ! just use the bin's dropout fraction

            call dmags( icmd, ix, iy, rr, dx, dy )
         endif

         ! Now that we have a dropout fraction, we can test for dropout
         call random(f,0.0D0,1.0D0)
         if (real(f).lt.p0) then
            dropflag(icmd) = 1
            dx = 0.0
            dy = 0.0
         endif

         ! Finally, apply the computed error offsets.
 200     x(icmd) = x(icmd) + dx
         y(icmd) = y(icmd) + dy
      end do

      return
      end

