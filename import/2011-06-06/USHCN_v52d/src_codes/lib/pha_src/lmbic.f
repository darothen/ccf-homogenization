      subroutine lmbic(iopt, xseg, yseg, end1, end2, qcrit, qstat,
     *  qmin, qoff, rmulm, slplm, amplm, philm, sseq, iqtype, mkfinq)

c      using the LMDIF routine from NETLIB that impliments
c        the non-linear curve fitting Levenberg-Marquardt algorithm
c        read in an x,y series file and fit to the following:

c        Y = (B1 + B2*X) + (B3 * SIN(B4 * (X - B5))) + E
c        where:
c         B1 is the intercept of the LINEAR function
c         B2 is the slope of the LINEAR function
c         B3 is the amplitude of the SINE function
c         B4 is the frequency of the SINE function
c         B5 is the phase shift of the SINE function

c     yseg(kmo) - serial input data array (with missing)
c     xseg(kmo) - serial input position (date) array
c     end1 - end of the first segment
c     end2 - end of the second segment

c     output is the Levenberg-Marquardt model fit with BIC test statistic

c     qcrit - 95% threshold of best model
c     qstat - test statistic at changepoint
c     qoff - offset using best model
c     iqtype - type of model
c     rmuq(2) - intercept for each segment of best model
c     slpq(2) - slope for each segment of best model
c     sseq - sum square error of the best model
c     mkfinq(2) - number of obs in each segment

c       IQTYPE 
c          8     lvmr - sloped lines with annual periodicity

c       Model comparison technique:
c       bic  - Bayesian Information Criteria test

      INCLUDE 'inhomog.parm.mthly.incl' 
      INCLUDE 'inhomog.comm.mthly.incl'
      INCLUDE 'inhomog.restart.mthly.incl'
      INCLUDE 'inhomog.MDparm.mthly.incl'

      integer end1, end2, beg2

c     given (mu + slp * xi) + amp * sin(lamda * (xi - phi) = yi 
c          for any given segment, then the variables are:
c       muq is the intecept for each segment
c       slpq is the slope
c       amp is the amplitude of the anuual sin component
c       phi is the phase offset of the sin component
c       lambda = 12 (for monthly data)
c       yseg is the yseged data series (that is: monotonically increasing
c         WITH missing data as required)
c       rfit is the fitted series
      real xseg(kmo), yseg(kmo), rfit(kmo), resid(kmo)
      real rmulm(2), slplm(2), amplm(2), philm(2)
      integer mkfinq(2)
      character*7 cmodel
     
      ibdebug = 2

c     the last index is the end of the series
      nobs = end2
c     the chgpt index is the end of the first segment
      ichgpt = end1
      do imo = end1 + 1, end2
        if(yseg(imo) .gt. amiss + 1.0) goto 10
      enddo  
   10 beg2 = imo
c      print *,' minbic end1, beg2, end2: ', end1, beg2, end2
      
c     save as original
      ichgpt0 = ichgpt
      ichgloop = ichgpt0
      if(iopt .gt. 1) 
     *  write(6, '(" QTYP     QVAL    QRSE     QPF     MU1     MU2",
     *    "  ALPHA1  ALPHA2   MSTAT   MCRIT    MOFF KNT1 KNT2")')

      call lmline(ichgpt, aMiss, xseg, yseg, xmed1, ymed1, slpme1, 
     *  Yintme1, Yampme1, Yphime1, SSElmqme1, mkntme1, iopt)
c      print *,ichgpt, xmed1, ymed1, slpme1, 
c     *  Yintme1, SSEflatme1, SSEslpme1, mkntme1

      call lmline(nobs-ichgpt, aMiss, xseg(ichgpt+1), yseg(ichgpt+1), 
     *  xmed2, ymed2, slpme2, Yintme2, Yampme2, Yphime2, SSElmqme2,
     *  mkntme2, iopt)
c      print *,nobs-ichgpt,  
c     *  xmed2, ymed2, slpme2, Yintme2, SSEflatme2, SSEslpme2, mkntme2

      cmodel = 'LMQSIN2'
      if(mkntme1 .ge. mintest .and. mkntme2 .ge. mintest) then
c       DO NOT USE the periodic component in the series offset estimate
        y1 = (yintme1 + slpme1 * xseg(ichgpt))
c     *   + Yampme1 * sin(lamda12 * (xseg(ichgpt) - Yphime1))
        y2 = (yintme2 + slpme2 * xseg(beg2))
c     *    + Yampme2 * sin(lamda12 * (xseg(beg2) - Yphime2))
        qoff = y1 - y2
        mKnt = mkntme1 + mkntme2
        SSEful = SSElmqme1 + SSElmqme2
        F2n4 = amiss
        ftpr = amiss
        call bayes(mknt, SSEful, 9, qtpr, rsq1, rsq2)
        if(ibdebug .gt. 1) 
     *    write(6,'(a,5f8.2,2f8.3,5f8.2,2i5)') cmodel, qtpr, rsq1, rsq2,
     *    yintme1, yintme2, slpme1, slpme2, yampme1, yampme2, yphime1, 
     *    yphime2, qoff, mkntme1, mkntme2
        qmin = qtpr
        iqtype = 8
        rmulm(1) = yintme1
        rmulm(2) = yintme2
        slplm(1) = slpme1
        slplm(2) = slpme2
        amplm(1) = yampme1
        amplm(2) = yampme2
        philm(1) = yphime1
        philm(2) = yphime2
        sseq = sqrt(SSEful/(mkntme1+mkntme2))
        mkfinq(1) = mkntme1
        mkfinq(2) = mkntme2
        qstat = F2n4
        qcrit = ftpr
        idfree = 9
      else
        if(ibdebug .gt. 1) 
     *    print *,cmodel,' - Unable to fit model - skipping'
      endif  

      return
      end
      

c ===================================================================


      subroutine lmline(inval, rMiss, xData, yData, xmean, ymean, slope,
     *  yInt, yamp, yphi, sselm, nval, iopt)
!-----------------------------------------------------------------------
! Subroutine to calculate levenberg-marquart non-linear regression line
! given an array of dependent (yData) and independent (xData) data 
! values.
!                          -on input-
!  inval - number input values including missing
!  rMiss - Missing value code (assumed to be a large negative number) 
!  xData - independent data (predictor)
!  yData - dependent data (predictand)
!
!                          -on output-
!  xmean - arithmetic average of x
!  ymean - arithmetic average of y
!  slope - slope of linear regression line
!  yInt - y - intercept of regression line
!  yamp - amplitude of sine conponent
!  yphi - phase of sine component
!  sselm - sum square error of l-m fit
!  resid - array of least squares residuals
!
!-----------------------------------------------------------------------      
      INCLUDE 'inhomog.parm.mthly.incl'
      parameter (nparm=4, lwa=nmo*nparm + 5*nparm + nmo)

      integer inval, nval, nslp
      real xmean,ymean,slope,yint,rMiss
      real xdata(nmo), ydata(nmo)
      
c     arrays for the fitting routines
c     paramter array
      double precision beta(nparm), tol
      integer iwa(nparm)
      double precision wa(lwa)

c     fitted array
      double precision fvec(nmo)

c     input data arrays
      double precision rX(nmo), rY(nmo), rX3(nmo), rY3(nmo)
c     try 3-month running average
      common rX3, rY3
c      common rX, rY

      integer maxfev,mode,mp5n,nfev,nprint

      double precision ftol,xtol,gtol,epsfcn,factor,zero
      data factor,zero /1.0d2,0.0d0/

c     user generated routine to calculate delta = obs-calc
      external fcn

      twopi = 2. * 3.141592653

      nval = 0
      xmean = rmiss
      ymean = rmiss
      slope = rmiss
      yint = rmiss
      yamp = rmiss
      yphi = rmiss
      sselm = 0.0
      
      sy = 0
      sx = 0
      ymin = 9999.
      ymax = -9999.
      do i=1,inval
        if ((xdata(i).gt.rMiss + 1).and.(ydata(i).gt.rMiss + 1)) then
          nval = nval + 1
          rx(nval) = xdata(i)
          ry(nval) = ydata(i)
          if(ydata(i) .gt. ymax) ymax = ydata(i)
          if(ydata(i) .lt. ymin) ymin = ydata(i)
          sx = sx + xdata(i)
          sy = sy + ydata(i)
        endif
      enddo        
      if(nval .lt. 5) return
      ymean = sy / nval
      xmean = sx / nval
      
c     try a three month running average
      yc = 0.0
      ycn = 0.0
      do i = 2, nval-1
        ry3(i-1) = (ry(i-1)+ry(i)+ry(i+1))/3.0
        rx3(i-1) = rx(i)
        if(ry(i-1) .lt. ymean .and. ry(i) .gt. ymean) then
          xc = rx(i)
          yc = yc + amod(xc, 12.)
          ycn = ycn + 1
        endif  
      enddo
      nval = nval -2

c     Initial estimate of the parameter array values
      beta(1) = ymean
      beta(2) = 0.0
      beta(3) = .3333 * (ymax - ymin)
      beta(4) = yc / ycn
      if(iopt .gt. 0) write(6,'(" BETA init:",i5,7f9.3)'),nval, xmean,
     *   ymean, beta
      tol = 0.001

c     some wrapper defaults
      maxfev = 200*(nparm + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      epsfcn = zero
      mode = 1
      nprint = 0
      mp5n = nval + 5*nparm

c      call lmdif1(fcn,nval,nparm,beta,fvec,tol,info,iwa,wa,lwa)
      call lmdif(fcn,nval,nparm,beta,fvec,ftol,xtol,gtol,maxfev,epsfcn,
     *  wa(1),mode,factor,nprint,info,nfev,wa(mp5n+1),nval,iwa,
     *  wa(nparm+1),wa(2*nparm+1),wa(3*nparm+1),wa(4*nparm+1),
     *  wa(5*nparm+1))
      if(iopt .gt. 0) write(6,'(" BETA:",2i5,7f9.3)'),nval, info,
     *   xmean, ymean, beta

!     calculate residuals
      do i=1,nval
          ylin = (beta(1) + beta(2)*rx(i))
          x = rx(i)
          ysin = beta(3) * sin((twopi/12.0) * (x - beta(4)))
          yest = ylin + ysin
          resid = yest - ry(i)
          sselm = sselm + resid ** 2
          if(iopt .gt. 0) write(6,'(" RESID: ",6f8.2)') rx(i),
     *      ry(i), ylin, ysin, yest, resid
      enddo 
      
      slope = beta(2)
      yint = beta(1)
      yamp = beta(3)
      yphi = beta(4)

      return
      end 

C     ************************************************************

      subroutine fcn(npts, nparm, beta, fvec, iflag)
      INCLUDE 'inhomog.parm.mthly.incl'

      integer npts,nparm,iflag
      double precision beta(nmo),fvec(nmo)

c     input point arrays
      double precision rX(nmo), rY(nmo)
      common rX, rY

c       Err = Y - (B1 + B2*X) + (B3 * SIN(2pi/12 * (X - B5)))
c     ----------
c     calculate the functions at x and
c     return this vector in fvec.

      twopi = 2. * 3.141592653

      do i = 1, npts
        x = rX(i)
        am = amod(x, 12.0)
        ophase = (twopi / 12.) * (am - beta(4))
        fvec(i) = rY(i) - (beta(1) + beta(2) * rX(i) + 
     *    beta(3) * sin(ophase))
      enddo

c     ----------
      return
      end
