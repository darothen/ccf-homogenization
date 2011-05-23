      subroutine hofnout(itarg, ntstn, nchg, ichgmo)

c  Originally the MLRTest2.f part of the HOFN Reference Series Changepoint
c     algorithm by Matt Menne.

c  Major rework for Graphics output of the HOFN Pairwise Changepoint algorithm
c     03 Apr 2008

c  This version differs from our original version in that the q-array 
c  is scaled by an annual cyle of standard deviations of the monthly 
c  means.  All stations, candidate and neighbors, are scaled by the 
c  same monthly standard deviation.  In this sense we can think of 
c  this annual cycle as representative of each station but not unique 
c  to it.  The purpose for doing this is to minimize the effect that 
c  generally higher variance of monthly means in winter than in summer 
c  could have on producing a non-stationary q-time series.
c
c  The first step is to calculate the correlation coefficients between
c  the candidate and each of its neighbors.  These calculations are 
c  performed in subroutines associated with the formerly separate program
c  RHODRVR.  These subroutines, now a part of this program, include:
c  hfrstdif, depave, & hcorrel.
c
c  A first difference filter is applied to the time series of monthly
c  mean temperatures.  Then the departures from the annual cycle of
c  monthly means of the differenced data for each station are computed.
c  Lastly, the candidate-neighbor correlation coefficients are calculated.
c  Missing data are taken into account in the calculations.
c
c  Correlation coefficients between the candidate and neighbors. 
c  effectively serve as station weights.
c  The annual cycle of monthly means is computed for each station
c  using an integer number of years.
c
c  The annual cycle of representive standard deviations of
c  monthly means is provided in subroutine MONSD.
c 
c  The q and z time series are computed in subroutine QZOFT2.
c  The q series is the monthly differences between the candidate
c  and the mean of its neighbors weighted by their respective
c  correlation coefficients and divided by the representative standard
c  deviations calculate in MONSD.  The z series is the standardized q
c  series.

      INCLUDE 'inhomog.parm.mthly.incl'
      INCLUDE 'inhomog.comm.mthly.incl'
      INCLUDE 'inhomog.restart.mthly.incl'

c     candidate-network array
      dimension rTarray(nmo,nstns)
      
c     Correlation between neighbors and candidate (see correl)
      dimension rho(nstns)
      
c     Mean annual cycle (monthly) for each month/station
      dimension rMac(12, nstns)
      
c     standard deviation of all stations for each month
      dimension rCNsd(12)
      
c     raw anomalies for plotting with XMGR
      dimension rAnom(nmo,nstns)
      
c     difference between cand & composite neighbors (q) and 
c        standardized (z) series
      dimension q(nmo),z(nmo)
      
c     Non standardized candidate and reference series
      dimension refSer(nmo),canSer(nmo)

c     temp arrays for data transfer
      dimension rDatOut(nmo,nstns),rDatTmp(nmo,nstns)

c     input parameters
c       montemp = cand & net temperature data, 12 monthly + annual
c       cCNids = cand & net coop numbers
      character*6 CCNids(nstns), ntstn(maxstns)
      
      character*4 cel3220(2)/'TMAX','TMIN'/

      integer ichgmo(maxstns, ninh)      
      dimension rMnDiff(ninh)
      
      integer beg2, end2

      mlag = 15
      
      do istn = 1, nstns
        ccnids(istn) = ''
      enddo  
      
      do ichg = 1, ninh
        rMnDiff(ichg) = 0.0
      enddo
      
      do imo = 1, nmo
        do  istn = 1, nstns
          rTarray(imo, istn) = amiss
        enddo
      enddo	  
      
c     decode the various faces of ihyear
c       if == 0 set to begyr else make ABS(ihyear)
c       Note: has already been checked to be >= begyr
      ifrstyr = ihyear
      if(ihyear .lt. 0) then
        ifrstyr = ihyear * -1
      else if(ihyear .eq. 0) then
        ifrstyr = begyr
      endif
      

C     populate the arrays needed for HOFN output processing
      irecent = 0
      do iy = ifrstyr+1, endyr
        do im = 1, 12
          call iym2imo(iy,im,imo)
          rtarray(imo,1) = temp(itarg, iy, im)
          if(iy .ge. 2006 .and. rtarray(imo,1) .ne. amiss) then
            irecent = 1
          endif
        enddo
      enddo
      
      if(irecent .eq. 0) then
        print *, ntstn(itarg), ' No data more recent than 2006'
        return
      endif  

      cCNIDs(1) = ntstn(itarg)
c     use 7 closest neighbors with 80% data overlap w/targ
      nNeigh = 2
      iStn = 2

      do while (iStn .le. 8 .and. nNeigh .le. nstns) 
        ipair = nindx(itarg,nNeigh)
c       skip out if there are no more neighbors
        if(ipair .eq. 0) goto 50
        
        iCnt = 0
        iDat = 0
     	  do iY = ifrstyr+1, EndYr
          do iM = 1,12
            if(temp(itarg,iY,iM) .ne. amiss) then
              iDat = iDat + 1
              if(temp(ipair,iy,im).ne. amiss) iCnt = iCnt + 1
            endif  
          end do
        end do
     
        if ((real(iCnt)/real(iDat)) .ge. 0.8) then
          cCNIDs(iStn) = ntstn(ipair)
          do iy = ifrstyr+1, endyr
            do im = 1, 12
              call iym2imo(iy,im,imo)
              rtarray(imo,iStn) = temp(ipair, iy, im)
            enddo
          enddo
          print *, 'OVR ',iStn,nNeigh,ipair,' ',ntstn(ipair),iCnt,iDat
       	  iStn = iStn + 1
          nNeigh = nNeigh + 1
        else
          print *,'XOVRLP ',nNeigh,ipair,' ',ntstn(ipair),iCnt,iDat
          nNeigh = nNeigh + 1
        endif
     
      enddo
      
   50 istn = istn - 1
      if(istn .lt. 4) then
        print *, ntstn(itarg), ' Not enough neighbors: ', istn
        return
      endif
      
c     Output raw temp for cand & neigh
      do iy = ifrstyr+1, endyr
        do im = 1, 12
          call iym2imo(iy,im,imo)
          if(rtarray(imo,1) .ne. amiss) then
            do iCNStn = 1, iStn
              if(rtarray(imo,iCNStn) .ne. amiss)
     *          write(47,'(a6,1x,a6,1x,a4,1x,i6,1x,f7.2)') cCNIds(1),
     *            cCNIds(iCNStn), cel3220(inel), iY*100+iM,
     *            rtarray(imo,iCNStn)
            enddo
          endif  
        enddo
      enddo
      
c     These are indices are used for rTarray
      iBegFul = 2
      iEndFul = nmo
      
      do iCNStn = 1, iStn  !begin station loop
c       Call to First Difference subroutine is done on a station by station basis
        call hfrstdif(nmo,rTarray(1,iCNStn),rDatTmp(1,iCNStn),amiss)
      enddo !end station loop

      do iCNStn = 1, iStn  !begin station loop
c       The annual cycle of monthly mean temps of the differenced
c       data for a station is computed, then the departures from the
c       annual cycle 
        call depave(1,nmo,rDatTmp(1,iCNStn),rDatOut(1,iCNStn),amiss)
      enddo !end station loop
  
      nyr = numyr

c     The temperature departure data (after differencing) are now
c     ready for calculating correlation coefficients between the
c     candidate and neighbors.
      call hcorrel(iStn,nstns,nmo,rDatOut,rho,1,nmo,amiss)
c      print *,'RHO: ',rho

c     Compute the mean annual cycle of the monthly means for each station
c     and the standard deviation of the monthly anomalies for the network
c     (the candidate and its neighbors)

c     use this subroutine to estimate variance of each monthly time series
      call monsd(istn,nstns,iBegFul,iEndFul,rTarray,rMac,rCNsd,amiss)

c     All the data needed by subroutine qzoft2 are now available
c     qzoft2 calculates the difference series q & z used by the algorithms
      call qzoft2(istn,nstns,iBegFul,iEndFul,rTarray,rMac,rCNsd,rho,
     &  q,z,rAnom,refSer,canSer,nmo,amiss)


C         PRINT EVERYTHING HERE !!!!!!!!!!!
      print *,'RMAC: ',(rMac(im,1),im=1,12)
      print *,' z:', z


c     Output the candidate, reference, and normalized series 
      do iy = ifrstyr+1, endyr
        do im = 1, 12
          call iym2imo(iy,im,imo)
          if(z(imo) .ne. amiss)
     *      write(36,'(a6,1x,a4,1x,i6,1x,3f7.2)') cCNIds(1),
     *      cel3220(inel), iY*100+iM, canSer(imo), refSer(imo), z(imo)
        enddo
      enddo
      
c     set the end of the first segment at the end POR
      do ichg = 1, nchg-1
        if(ichg .eq. 1) then
          beg2 = ichgmo(itarg,ichg)
        else  
          beg2 = ichgmo(itarg,ichg) + 1
        endif
        end2 = ichgmo(itarg, ichg+1)
        
        zSum = 0.0
        qSum = 0.0
        n = 0
        print *, 'imo: ', itarg, ichg, beg2, end2
        do imo = beg2, end2
          if(z(imo) .ne. amiss) then
            zSum = zSum + z(imo)
            qSum = qSum + (canSer(imo) - refSer(imo))
            n = n + 1
          endif  
        enddo
        
c       write out the mean of the normalized segment series
        if(n .gt. 0) then
          zMean = zSum / n
          qMean = qSum / n
          rMnDiff(ichg) = qMean
          do imo = beg2, end2
            call imo2iym(iy,im,imo)
            if(z(imo) .ne. amiss) then
              write(46,'(a6,1x,a4,1x,i6,1x,f7.2)') cCNids(1),
     &          cel3220(inel), (iy*100)+im, zMean
            endif  
          enddo
        endif  
      enddo
      
c     write out the changepoints offset from one segment to the next
      if(nchg .gt. 2) then
        do ichg = 1, nchg-1
          call imo2iym(iy,im,ichgmo(itarg,ichg+1))
          write(39,'(a6,1x,a4,1x,i6,f7.2)') cCNids(1),
     *      cel3220(inel), (iy*100)+im, rMnDiff(iChg+1)-rMnDiff(iChg)
        enddo
      endif  

      return
      end

c
c     =======================================================================

      subroutine qzoft2 (istn,ns,nx1,nx2,x,rMac,yxsd,rho,q,z,
     *  rAnom,refSeries,canSeries,nmo,rmiss)

c***********************************************************************
c
c     Subroutine qzoft computes q and z time series as defined by Eqs. (2)
c     and (3) in Alexandersson and Moberg (Int. J. Climatol., 17, 25-34, 
c     1997) in which the input variable is monthly mean temperature.
c     Variable q is called the "differenced" series and z the "standarized"
c     series.  The z series is used in the test statistics ts (shift or
c     step change) and tt (trend).
c
c     What is different about sbqzoft2 with respect to sbqzoft is that
c     the q array is scaled by the annual cycle of a standard deviations of
c     monthly means.  This was done because, particularly in mid-latitude 
c     stations, there is greater variance of monthly means in winter than
c     in summer.  Without accounting for this variation with month of year
c     the q values are nonstationary.  Whether the impact on he SNIT T-
c     statistic of not scaling q is minor or major, it makes good sense,
c     mathematically, to create as stationary a data set as possible.
c     This modification was done in July, 2001.
c 
c     September 21, 1999; 30 November, 1999; 10 July, 2001
c
c                               -on input-
c
c        ns    The number of stations (candidate & neighbors)
c
c       nx1    Beginning of the monthly series
c
c       nx2    Ending of the monthly series
c
c         x    The total months by ns array of raw data 
c
c      rMac    The 12xm array of average monthly mean temperatures
c               for each month of the year for the m neighbor stations
c
c       rho     The m by 1 array of Pearson correlation coefficients
c               between the candidate station and the neighbor 
c               stations
c
c      yxsd     The 12 by 1 array of monthly standard deviations from
c               subroutine monsd
c
c                               -on output-
c
c       q       The n by 1 array of differences between the candidate
c               and the average of the neighbors as defined by Eq.(2)
c               referenced above
c
c       z       The n by 1 array of the standardized q series as
c               defined by Eq. (3) referenced above
c
c              The following variables are passed for print output
c
c      rAnom    Monthly anomalies for all stations
c
c     refSer    Weighted Network monthly anomalies
c
c     canSer    Monthly anomalies for candidate station
c
c                               -comments-
c
c               Arrays y, x, q, and z all begin in January of a
c               given year and end in the current month of
c               analysis in a later year
c
c               The annual cycle of ensemble average mean monthly
c               temperatures and the Pearson correlation coefficients
c               are calculated or recalculated annually only when
c               a new batch of 12 monthly temperatures,
c               January thru December, are available.
c
c               Only the y and x arrays on input have missing data
c
c************************************************************************
      dimension x(nmo,ns)
      dimension rMac(12,ns),rho(ns)
      dimension yxsd(12)

c     the concatenated arrays
      dimension q(nmo),z(nmo)
      dimension rAnom(nmo,ns)
      dimension refSeries(nmo),canSeries(nmo)
c
      sumq = 0.0
      sumq2 = 0.0
      numq = 0
      
      do im = 1, nmo
        q(im) = rMiss
        z(im) = rMiss
        refSeries(im) = rMiss
        canSeries(im) = rMiss
        do is = 1, ns
          rAnom(im,is) = rMiss
        enddo  
      enddo  
c
c     Compute sum of squares of correlation coefficients
c     A value is computed for each nth value to account for the 
c     reapportionment of the weights where there are missing
c     observations. The rhosqrd series is "compressed" with no
c     observations with the candidate missing or less than
c     2 stations in the network

      do ix =  1, nmo
        call imo2iym(iy,im,ix)
        
        if(x(ix,1) .ne. rMiss) then
          rAnom(ix,1) = x(ix,1) - rMac(im,1)
          rhosqrd = 0.0
          nrho = 0
          sum = 0.0
          do is = 2, istn
            if(x(ix,is) .ne. rMiss) then
c             Compute raw anomalies to be used for viewing in XMGR
              rAnom(ix,is) = x(ix,is) - rMac(im,is)
              sum = sum + rho(is)**2*(rAnom(ix,is))
              rhosqrd = rhosqrd + rho(is)**2
              nrho = nrho + 1
            else
              rAnom(ix,is) = rMiss  
            endif
          enddo
c         see if there are at least 1 stations in the network
          if(nrho .ge. 1) then
c           Compute the q series 
            q(ix)=(rAnom(ix,1)-sum/rhosqrd)/yxsd(im)
            sumq = sumq + q(ix)
            sumq2 = sumq2 + q(ix) * q(ix)
            numq = numq + 1
            refSeries(ix) = sum/rhosqrd
            canSeries(ix) = rAnom(ix,1)
          endif
        endif    
      end do

c     Compute the mean and standard deviation of the q series
      if(numq .gt. 1) then
        qmn = sumq/numq
        qsd = sqrt((sumq2 - (sumq * sumq / numq)) / (numq - 1))
      endif  
c
c     Compute z series
c
      do i = 1,nmo
        if(q(i) .ne. rMiss) z(i) = (q(i)-qmn)/qsd
      enddo

      return
      end

c     =======================================================================

      subroutine monsd (istn,ns,nx1,nx2,tdata,rmac,yxsd,rmiss)
      
c***************************************************************************
c
c     Subroutine monsd produces a 12 by 1 array of standard deviations of 
c     monthly means that is representative of all stations (candidate and 
c     neighbors). It is used to account for the annual cycle in the variance
c     of monthly means, which is larger in winter than summer, especially at
c     middle latitude stations.
c 
c     The standard deviations are used in subroutine qzoft2 to scale the
c     q-variable.
c
c     10 July, 2001
c
c                               -on input-
c
c  tdata  The month by station array of mean monthly temperatures.
c
c    mac  The (12, ns) array of 12 mean monthly means for the candidate 
c           and neighbors.
c
c    istn actual number of stations (including candidates)
c
c     ns  maximum number of stations (candidate plus neighbors).
c
c    nx1  Beginning month in the time series
c
c    nx2  Ending month in the time series
c
c                               -on output-
c
c   rmac  The 12 x ns array of monthly means for each station
c
c   yxsd  The 1-d arrary comprising the annual cycle of standard deviations
c         common to all stations to be used in qzoft.
c
c***************************************************************************
c
      dimension rmac(12,ns),yxsd(12)
      dimension tdata(nx2-nx1+2,ns)
      dimension rsum(12), rnum(12)
      
      do im = 1, 12
        rSum(im) = 0.0
        rNum(im) = 0.0
      enddo  

c     Compute mean annual cycle of monthly means using an
c     integer number of years = nyears
      do is = 1, istn
        do imo = nx1,nx2
          call imo2iym(iy, im, imo)
          if(tdata(imo,is).ne.rMiss) then
            rSum(im) = rSum(im)+tdata(imo,is)
            rNum(im) = rNum(im)+1
          endif
        enddo
        do im = 1,12
          if(rNum(im) .ne. 0) then
            rMac(im,is) = rSum(im)/rNum(im)
          else
            rMac(im,is) = rMiss
          endif
          rSum(im) = 0.0
          rNum(im) = 0
          enddo
      enddo  

c     Compute the standard deviation for each of the 12 months that is repre-
c     sentative of all stations.  The departures in the formula are calculated
c     with respect to the mean monthly means for each station as given in array
c     yxmn above. 
      do im = 1,12
        rSum(im) = 0.0
        rNum(im) = 0
      enddo
      
      do imo = nx1, nx2  
        call imo2iym(iy, im, imo)
        do is = 1,istn
          if(tdata(imo,is) .ne. rMiss) then
            rSum(im) = rSum(im) + (tdata(imo,is) - rmac(im,is))**2
            rNum(im) = rNum(im) + 1
          endif
        enddo
      enddo
      
      do im = 1, 12  
        yxsd(im) = sqrt(rSum(im)/rNum(im))      !the 12 stnd deviations
      enddo
      return
      end 

c     =======================================================================
c     =======================================================================
c     =======================================================================

      subroutine hfrstdif (nx,x,y,rmiss)
c*****************************************************************
c
c     Subroutine hfrstdif computes a first difference time series of
c     monthly mean temperatures taking into account missing monthly
c     values.  7 September, 1999; 26 November, 1999
c
c                               -on input-
c
c       nx      The number (odd) of monthly mean temperatures
c
c       x       The 2-d (ns, nx) monthly mean temperatures
c
c                               -on output-
c
c       ny      The number (even) of first differences ny = nx -1
c
c       y       The 2-d (ns, ny) first differences
c
c******************************************************************
      dimension x(nx)
      dimension y(nx-1)
      kount = 0

c     initialize
      do im = 1, nx - 1
        y(im) = rMiss
      enddo

c
c     find first good value
c
      do im1 = 1,nx
        if (x(im1).ne.rMiss) goto 5
      enddo
      goto 99    
c
c     do first difference filter accounting for missing data
c
    5 kount = 1
      itot = 0
      do im = im1+1,nx
        if (x(im) .eq. rMiss) then
          kount = 0
        else
          kount = kount + 1
          if (kount.gt.1) then
            y(im-1) = (x(im) - x(im-1))/2
            itot = itot + 1
          endif  
        endif
      enddo
c      print *,' Frstdif: ', itot  
 99   return
      end     

c     =======================================================================

      subroutine depave (nx1,nx2,x,y,rmiss)
c****************************************************************
c
c     Subroutine depave (1) computes the average annual cycle of  
c     monthly mean temperatures from a time series
c     of monthly mean temperatures comprising an integer number of 
c     years and (2) calculates the departures from the average
c     annual cycle.  The input time series can have missing data.
c     Richard Jones method is used (see p. 33 of Chapter 1, Fourier
c     Analysis, Metr. 5323).  8 September, 1999; 27 November, 1999
c
c                               -on input-
c
c       nx       The number of monthly mean temperatures 
c       
c       x        The 1-d (nx) monthly mean temperatures
c
c                               -internal-
c
c       amean    The 1-d array of the 12 average monthly mean 
c                temperatures
c
c                               -on output-
c
c       y        The 1-d (ny) departures from the mean
c                annual cycle for a candidate and neighbors
c 
c       nogyrs   The 1-d array of the number of years used in 
c                computing each of the 12 average monthly mean
c                temperatures
c
c       Note     While reference is made to monthly temperatures, 
c                the subroutine can be applied to any time series of
c                monthly values that have an annual cycle.
c
c*****************************************************************

      dimension amean(12),x(nx2-nx1+1),y(nx2-nx1+1)
      dimension sum(12)
      integer num(12)
      
      do im = 1, nx2-nx1+1
        y(im) = rMiss
      enddo  
      
      do im = 1, 12
        sum(im) = 0.0
        num(im) = 0
      enddo  
c
c     calculate the average monthly mean temperature for each of the
c     12 months of the year
c
      do imo = nx1,nx2
        call imo2iym(iy,im,imo)
        if (x(imo).ne.rMiss) then
          sum(im) = sum(im) + x(imo)
          num(im) = num(im) + 1
        endif
      enddo  

      do im = 1, 12
        if(num(im) .ge. 4) then
          amean(im) = sum(im)/num(im)
        else
          amean(im) = rMiss
        endif
      enddo
c      print *,' Dpave amean: ', amean
c
c     compute the monthly departures from the annual cycle of average
c     monthly mean temperatures calculated above
c
      yndx = 1
      do imo = nx1,nx2
        call imo2iym(iy,im,imo)
        if (x(imo).eq.rMiss.or.amean(im).eq.rMiss) then
          y(imo) = rMiss
        else
          y(imo) = x(imo) - amean(im)
        endif
      enddo
      return
      end

c     =======================================================================

      subroutine hcorrel(istn,ns,nm,x,ccoef,nx1,nx2,rmiss)
c********************************************************************
c
c     Subroutine hcorrel computes the Pearson correlation coefficent
c     between a candidate station and each of its neighbors.  Account
c     is taken of missing data in the candidate and each neighbor so
c     that the mean, stnd dev, and cross-products are computed only
c     for those months when both candidate and neighbor have good,
c     as opposed to missing, data.  Each correlation coefficient applies
c     to all months of the year.
c     9 September, 1999; 29 November, 1999
c
c                               -on input-
c
c      istn     Actual number of stations (including candidate)
c
c       ns      Maximum number of stations (including candidate)
c
c       nm      number of months of good and missing data in each 
c               time series of monthly departures.
c
c       x       the (nm x ns) array of monthly departures for the candidate
c               and successive neighbor stations - each array sum is zero
c
c                               -on output-
c
c       ccoef   the (ns-1 x 1) array of candidate-neighbor Pearson correlation
c               coefficients
c
c       mon     the (ns-1 x 1) array of the number of months used in calculating
c               each of the m correlation coefficients
c
c       prcnt   the (ns-1 x 1) array of the percentages of the maximum possible
c               number of months that could be used to calculate each
c               of the m correlation coefficients, i.e., (mon(j)/n)x100%
c
c********************************************************************
      dimension x(nm,ns),ccoef(ns)
      dimension xn(nm),yn(nm)
      
      do is = 1, ns
        ccoef(is) = rMiss
      enddo  
      
c
c     Determine those months that both candidate and a neighbor have
c     good data.  They will vary with the particular candidate-neighbor
c     combination.  
c
      ccoef(1) = 1.0
      do 10 is = 2,istn
        kk = 0       !the counter for good candidate-neighbor data
        do im = nx1,nx2
          if (x(im,1).ne.rMiss.and.x(im,is).ne.rMiss) then
            kk = kk + 1
            xn(kk) = x(im,1)
            yn(kk) = x(im,is)
          endif  
        enddo  
c
c       calculate means for candidate and neighbor using only those months
c       that both neighbor and candidate have good data
c
        xmn = 0.0
        ymn = 0.0
        do 30 k = 1,kk
          xmn = xmn + xn(k)
          ymn = ymn + yn(k)
 30     continue
        if(kk .gt. 4) then
          xmn = xmn/kk
          ymn = ymn/kk
        else
          xmn = rMiss
          ymn = rMiss
          goto 10
        endif
c  
c       calculate standard deviations for candidate and neigbhor using
c       only those months that both neighbor and candidate have good
c       data
c
        xsd = 0.0
        ysd = 0.0
        do 40 k = 1,kk
          xsd = xsd + (xn(k)-xmn)**2
          ysd = ysd + (yn(k)-ymn)**2
 40     continue
        xsd = sqrt(xsd/(kk-1))
        ysd = sqrt(ysd/(kk-1))
c
c       calculate Pearson correlation coefficient for each candidate-
c       neighbor combination
c
        sum = 0.0
        do 50 k = 1,kk
          sum = sum + (xn(k)-xmn)*(yn(k)-ymn)
 50     continue
        ccoef(is) = sum/((kk-1)*xsd*ysd)
 10   continue
      return
      end
      
