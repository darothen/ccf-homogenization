c     =======================================================================
c
c     version   date     description
c
c    --------------------------- Official Version USHCNv2 ----------------------
c
c      v6b     08jan07   last pass of ESTAMT uses ONLY TPR0 for best adj 
c                          consistency (iopt=6)

      subroutine minbic(iopt, xseg, yseg, end1, end2, qcrit, qstat,
     *  qmin, qoff, rmuq, slpq, sseq, inqtype, iqtype, mkfinq, ifail)

c     routine to select the most appropriate model for a given single
c       changepoint along with the adjacent segments using the
c       Bayesian Information Criteria (BIC)
c
c     For the "best" model to be selected the BIC must be less than 
c       that of simpler models.
c     Also, models with slopes:
c       1) a minimum of 5 years (60 months) of raw data is required and
c       2) the t-test for the slope parameter(s) must be significant
c     kmo - maximum size of the data/position arrays
c     amiss - missing value

c     input:
c     iopt - Process option 
c       1 = best model including SLR's & TPR's 
c             (used on interim chgpts)
c       2 = best chgpt amplitude TPR's only (except for TPR0 vs SLR1)
c             (used in splitmerge for model ID and comparison)
c             !!!!! NOTE: as of v6a - returns ONLY TPR0 fit!!!!!!
c       3 = best model plus segment intercept test - looking for break slope
c       4 = same as iopt=2 with slope AND critical value failure test
c             (used in estamt to get best amp est from all available neighbors)
c       5 = same as iopt=1 with slope AND critical value failure test
c             (used in estamt to test unsupport metadata chgpts)
c       6 = last pass of ESTAMT uses ONLY TPR0 for best adjustment consistency
c       NOTE: for iopt=4-6 (ESTAMT options) minslptest == minslpknt
c                otherwise minslptest == mintest (== 2)
c       NOTE: all iopt
c             include critical value failure test (no slope test 11may06)
c        For each iopt set these two on/off (1/0) arrays for testing
c          ictest for significant threshold and istest for realistic slope
c     yseg(kmo) - serial input data array (with missing)
c     xseg(kmo) - serial input position (date) array
c     end1 - end of the first segment
c     end2 - end of the second segment
c
c     output from the best (lowest) model using the BIC test statistic
c     qcrit - 95% threshold of best model
c     qstat - test statistic at changepoint
c     qoff - offset using best model
c     inqtype - type of stat test model used for current changepoints
c     iqtype - type of model
c     rmuq(2) - intercept for each segment of best model
c     slpq(2) - slope for each segment of best model
c     sseq - sum square error of the best model
c     mkfinq(2) - number of obs in each segment
c
C     Chgptmodels is a consolidation of all the statistical tests and their
c       attendant routines for the UCP. 
c     An effort has been made to standardize the tests for simplicity sake
c
c
c     Current models:
c
c       IQTYPE 
c              Homogeneous linear types:
c         1    slr0 - simple flat line
c         2    slr1 - sloped straight line
c
c              Two phase regression types:
c         3    tpr0 - amplitude only, 0 sloped segments 
c         4    tpr1 - amplitude shift with equal sloped segments
c         5    tpr2 - amplitude shift and non-equal sloped segments
c         6    tpr3 - amplitude shift with flat-to-slope segments
c         7    tpr4 - amplitude shift with slope-to-flat segments
c
c       Model comparison technique:
c       bic  - Bayesian Information Criteria test


      INCLUDE 'inhomog.parm.mthly.incl' 
      INCLUDE 'inhomog.comm.mthly.incl'
      INCLUDE 'inhomog.restart.mthly.incl'
      INCLUDE 'inhomog.MDparm.mthly.incl'

      integer end1, end2, beg2

c     given mu + slp * xi = yi for any given segment, then the variables
c       may be generalized such that:
c       muq is the intecept for each segment
c       slpq is the slope
c       yseg is the yseged data series (that is: monotonically increasing
c         WITH missing data as required)
c       rfit is the fitted series
      real xseg(kmo), yseg(kmo), rfit(kmo), resid(kmo)
      real rmuslr(2), rmutpr(2), slptpr2(2)
      real rmuq(2), slpq(2), rmut(2), slpt(2)
      integer mkfinq(2), mkfint(2)
      character*7 cmodel
      
c     For each iopt set these two on/off (1/0) arrays for testing
c       ictest for significant threshold and istest for realistic slope
      integer ictest(6)/1,1,1,1,1,0/,istest(6)/0,0,0,1,1,0/
      
      mintest = 2
c      mintest = minlen
      if(iopt .lt. 4) then 
        minslptest = mintest
      else
        minslptest = minslpknt
      endif  
      
c     isloop keeps track of the loop-back for the segment intercept test
c     iopt == 3 will test the BIC Q value ONLY when:
c       1) One of the Model types "flat-slope", "slope-flat", and "slope-slope" 
c             is chosen as the best model
c       2) The estimated intercept is within the series limits
      isloop = 0
      
      ifail = 0
      
      ibdebug = 0
c      if(iopt .ge. 4) ibdebug = 2
      if(ibdebug .ge. 1) 
     *  print *,' CHGPTMODEL ingtype:',inqtype,' iopt:',iopt 
      
c     lsqopt (=1) prints the calculations using least squares NOT the solution
c       which always comes from the Kendall-Theil paired-ratio approach
      lsqopt = 0

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
      if(ibdebug .gt. 1) 
     *  write(6, '(" QTYP     QVAL    QRSE     QPF     MU1     MU2",
     *    "  ALPHA1  ALPHA2   MSTAT   MCRIT    MOFF KNT1 KNT2")')

c     run full series as one segment (straight line)
      qoff = 0.0
      if(lsqopt .eq. 1) 
     *  call lsqline(nobs, aMiss, xseg, yseg, xmean, ymean, slope, 
     *    Yint, sseflat, sseslope, mknt)
c      print *,nobs, xmean, ymean, slope, 
c     *  Yint, sseflat, sseslope, mknt
      call kthline(nobs, aMiss, xseg, yseg, xmed, ymed, slpmed, 
     *  Yintmed, ssefltmed, sseslpmed, mknt)
c      print *,nobs, xmean, ymean, slope, 
c     *  Yint, sseflat, sseslope, mknt
c     save the "reduced" SSE for the TPR models
      SSEred = sseslope
      SSEredmed = sseslpmed
      qmin = 9999.
      tmin = 9999.
      estloop = 99.

c     if iopt == 1 or 5 then find the best of ANY model
      if(iopt .eq. 1 .or. iopt .eq. 5) goto 20
      
c     for all other iopt values, if TPR0 include the SLR1 test
      if(inqtype .eq. 3) go to 50

c     For all other TPR types use sloped line fit to get SSEred but
c       NOT for a better model
      if(lsqopt .eq. 1 .and. ibdebug .gt. 1) 
     *  print *,' SLR1 yint, slope, SSEred:',yint, slope, SSEred
c     use the kendall-theil method (with slope calc)
      if(ibdebug .gt. 1) 
     *  print *,' KTH1 yint, slope, SSEred:',yintmed,slpmed,SSEredmed
      iqtype = 2
      rmuq(1) = Yintmed
      rmuq(2) = Yintmed
      slpq(1) = slpmed
      slpq(2) = slpmed
      sseq = sqrt(SSEredmed/mknt)
      mkfinq(1) = mknt
      mkfinq(2) = 0
      qstat = 0.0
      qcrit = 99.0
      ittype = 2
      rmut(1) = Yintmed
      rmut(2) = Yintmed
      slpt(1) = slpmed
      slpt(2) = slpmed
      sset = sqrt(SSEredmed/mknt)
      mkfint(1) = mknt
      mkfint(2) = 0
      tstat = 0.0
      tcrit = 99.0
      goto 100

c     The 2 straight line BIC models use the residual of the entire 
c       series as one segment - first the flat straight line
   20 if(lsqopt .eq. 1) then
        cmodel = 'FITSLR0'
        if(mknt .ge. mintest) then
          call bayes(mknt, sseflat, 1, qslr0, rsq1, rsq2)
          if(ibdebug .gt. 1) 
     *      write(6,'(a,4f8.2,5(" -------"),f8.2,i5," ----")')
     *        cmodel, qslr0, rsq1, rsq2, ymean, qoff, mknt
        else
          if(ibdebug .gt. 1) 
     *      print *,cmodel,' - Unable to fit model - skipping'
        endif  
      endif  

c     use the kendall-theil method (without slope calc)
      cmodel = 'KTHSLR0'
      if(mknt .ge. mintest) then
        call bayes(mknt, ssefltmed, 1, qslr0, rsq1, rsq2)
        qmin = qslr0
        iqtype = 1
        rmuq(1) = ymed
        rmuq(2) = ymed
        slpq(1) = 0.0
        slpq(2) = 0.0
        sseq = sqrt(sseflat/mknt)
        mkfinq(1) = mknt
        mkfinq(2) = 0
        qstat = 0.0
        qcrit = 99.0
        tmin = qslr0
        ittype = 1
        rmut(1) = ymed
        rmut(2) = ymed
        slpt(1) = 0.0
        slpt(2) = 0.0
        sset = sqrt(sseflat/mknt)
        mkfint(1) = mknt
        mkfint(2) = 0
        tstat = 0.0
        tcrit = 99.0
        if(ibdebug .gt. 1) 
     *    write(6,'(a,4f8.2,5(" -------"),f8.2,i5," ----")')
     *      cmodel, qslr0, rsq1, rsq2, ymed, qoff, mknt
      else
        if(ibdebug .gt. 1) 
     *    print *,cmodel,' - Unable to fit model - skipping'
      endif  
      
c     the sloped line (also gives back the Sum Square Error for 0 chngpts)
   50 if(lsqopt .eq. 1) then
        cmodel = 'FITSLR1'
        if(mknt .ge. mintest) then
          call bayes(mknt, SSEred, 2, qslr1, rsq1, rsq2)
          if(ibdebug .gt. 1) 
     *      write(6,'(a,4f8.2," -------",f8.3,3(" -------"),f8.2,
     *        i5," ----")') cmodel, qslr1,  rsq1, rsq2, Yint, slope, 
     *        qoff, mknt
        else
          print *,cmodel,' - Unable to fit model - skipping'
        endif  
      endif  

c     use the kendall-theil method (with slope calc)
      cmodel = 'KTHSLR1'
      if(mknt .ge. minslptest) then
        call bayes(nobs, SSEredmed, 2, qslr1, rsq1, rsq2)
        if(ibdebug .gt. 1) 
     *    write(6,'(a,4f8.2," -------",f8.3,3(" -------"),f8.2,
     *      i5," ----")') cmodel, qslr1,  rsq1, rsq2, Yintmed,
     *      slpmed, qoff, mknt
c       following line replaced to reposition the t-test AFTER all BICs
        if(qslr1.lt.qmin) then
          qmin = qslr1
          iqtype = 2
          rmuq(1) = Yintmed
          rmuq(2) = rmuq(1)
          slpq(1) = slpmed
          slpq(2) = slpq(1)
          sseq = sqrt(SSEredmed/mknt)
          mkfinq(1) = mknt
          mkfinq(2) = 0
          qstat = 0.0
          qcrit = 99.0
          if(istest(iopt).eq.0 .or. abs(slpmed).lt.slpmthmax) then
            tmin = qslr1
            ittype = 2
            rmut(1) = Yintmed
            rmut(2) = rmuq(1)
            slpt(1) = slpmed
            slpt(2) = slpq(1)
            sset = sqrt(SSEredmed/mknt)
            mkfint(1) = mknt
            mkfint(2) = 0
            tstat = 0.0
            tcrit = 99.0
          endif
        endif
      else
        if(ibdebug .gt. 1) 
     *    print *,cmodel,' - Unable to fit model - skipping'
      endif

C     for the TPR models - the series is divided into two segments
c       from 1 to chgpt and from chgpt+1 to nobs - for all the models
c       EXCEPT for TPR4 (constant slope with one offset) the individual
c       segment info change be obtained once
c     first segment (starting at 1 for ichpts)
  100 if(lsqopt .eq. 1)
     *  call lsqline(ichgpt, aMiss, xseg, yseg, xmean1, ymean1, slope1,
     *  Yint1, SSEflat1, SSEslope1, mknt1)
c      print *,ichgpt, xmean1, ymean1, slope1, 
c     *  Yint1, SSEflat1, SSEslope1, mknt1
      call kthline(ichgpt, aMiss, xseg, yseg, xmed1, ymed1, slpme1, 
     *  Yintme1, SSEflatme1, SSEslpme1, mkntme1)
c      print *,ichgpt, xmed1, ymed1, slpme1, 
c     *  Yintme1, SSEflatme1, SSEslpme1, mkntme1
c     second segment (starting at ichgpt+1 for nobs-ichgpt)
      if(lsqopt .eq. 1)
     *  call lsqline(nobs-ichgpt, aMiss, xseg(ichgpt+1), yseg(ichgpt+1),
     *    xmean2, ymean2, slope2, Yint2, SSEflat2, SSEslope2, mknt2)
c      print *,nobs-ichgpt,  
c     *  xmean2, ymean2, slope2, Yint2, SSEflat2, SSEslope2, mknt2
      call kthline(nobs-ichgpt, aMiss, xseg(ichgpt+1), yseg(ichgpt+1), 
     *  xmed2, ymed2, slpme2, Yintme2, SSEflatme2, SSEslpme2, mkntme2)
c      print *,nobs-ichgpt,  
c     *  xmed2, ymed2, slpme2, Yintme2, SSEflatme2, SSEslpme2, mkntme2

c     if isloop == 1 then iopt == 3 and the best model was TPR3,4,or 5
c       for each type recalculate offset and if less than previous
c       find new offset and recheck offset
      if(isloop .eq. 1) then
        if(iqtype .eq. 5) then
c         full TPR (slope-slope)
          cmodel = 'KTHBRK2'
          est = yintme1-yintme2 + (slpme1-slpme2)*xseg(ichgpt)
          if(abs(est) .lt. abs(estloop)) then
            qoff = est
            ichgloop = ichgpt
            rmuloop1 = yintme1
            rmuloop2 = yintme2
            slploop1 = slpme1
            slploop2 = slpme2
            mkloop1 = mkntme1
            mkloop2 = mkntme2
            SSEfloop = SSEslpme1 + SSEslpme2
            idloop = 5
            goto 180
          endif
        else if(iqtype .eq. 6) then
c         flat-slope model
          cmodel = 'KTHBRK3'
          est = ymed1 - (yintme2 + xseg(ichgpt) * slpme2)
          if(abs(est) .lt. abs(estloop)) then
            qoff = est
            ichgloop = ichgpt
            rmuloop1 = ymed1
            rmuloop2 = yintme2
            slploop1 = 0.0
            slploop2 = slpme2
            mkloop1 = mkntme1
            mkloop2 = mkntme2
            SSEfloop = SSEflatme1 + SSEslpme2
            idloop = 4
            goto 180
          endif
        else if(iqtype .eq. 7) then
c         slope-flat model
          cmodel = 'KTHBRK4'
          est = (yintme1 + xseg(ichgpt) * slpme1) - ymed2
          if(abs(est) .lt. abs(estloop)) then
            qoff = est
            ichgloop = ichgpt
            rmuloop1 = yintme1
            rmuloop2 = ymed2
            slploop1 = slpme1
            slploop2 = 0.0
            mkloop1 = mkntme1
            mkloop2 = mkntme2
            SSEfloop = SSEslpme1 + SSEflatme2
            idloop = 4
            goto 180
          endif
        else  
          print *,' Incorrect model for SIT loop ', iqtype
          stop
        endif
c       last iteration minimized offset - calculate breakpoint SSE
        if(ibdebug .gt. 1) 
     *    print *,' est gt toofloop: ', est
        goto 190          
      endif      

c     step change only, 0 sloped segments
  130 if(lsqopt .eq. 1) then
        cmodel = 'FITTPR0'
        if(mknt1 .ge. mintest .and. mknt2 .ge. mintest) then
          call ttest(kmo, yseg, end1, end2, tchgpt, mintest)
          ttpr0 = critval(mknt-2, 3)
          SSEful = SSEflat1 + SSEflat2
          call bayes(mknt, SSEful, 3, qtpr, rsq1, rsq2)
          est = ymean1 - ymean2
          if(ibdebug .gt. 1) 
     *      write(6,'(a ,5f8.2,2(" -------"),3f8.2,2i5)')
     *        cmodel, qtpr, rsq1, rsq2, ymean1, ymean2, tchgpt, ttpr0,
     *        est, mknt1, mknt2
        else
          if(ibdebug .gt. 1) 
     *      print *,cmodel,' - Unable to fit model - skipping'
        endif  
      endif
      
c     use the kendall-theil method (with 0 sloped segments)
      cmodel = 'KTHTPR0'
      if(mkntme1 .ge. mintest .and. mkntme2 .ge. mintest) then
        call ttest(kmo, yseg, end1, end2, tchgpt, mintest)
        ttpr0 = critval(mknt-2, 3)
        SSEful = SSEflatme1 + SSEflatme2
        call bayes(mknt, SSEful, 3, qtpr, rsq1, rsq2)
        est = ymed1 - ymed2
        if(ibdebug .gt. 1) 
     *    write(6,'(a,5f8.2,2(" -------"),3f8.2,2i5)')
     *      cmodel, qtpr, rsq1, rsq2, ymed1, ymed2, tchgpt, ttpr0, est,
     *      mkntme1, mkntme2
        if(qtpr.lt.qmin) then
          qmin = qtpr
          iqtype = 3
          rmuq(1) = ymed1
          rmuq(2) = ymed2
          slpq(1) = 0.0
          slpq(2) = 0.0
          sseq = sqrt(SSEful/(mkntme1+mkntme2))
          mkfinq(1) = mkntme1
          mkfinq(2) = mkntme2
          qstat = tchgpt
          qcrit = ttpr0
          qoff = est
          if(ictest(iopt) .eq. 0 .or. tchgpt .ge. ttpr0) then
            tmin = qtpr
            ittype = 3
            rmut(1) = ymed1
            rmut(2) = ymed2
            slpt(1) = 0.0
            slpt(2) = 0.0
            sset = sqrt(SSEful/(mkntme1+mkntme2))
            mkfint(1) = mkntme1
            mkfint(2) = mkntme2
            tstat = tchgpt
            tcrit = ttpr0
            toff = est
          endif
        endif
      else
        if(ibdebug .gt. 1) 
     *    print *,cmodel,' - Unable to fit model - skipping'
      endif  
      
c     for iopt = 6 return only TPR0 solution
      if(iopt .eq. 6) goto 200

c     skip TPR1 when estamt needs consistency between slope segments
c      if(iopt .ge. 4) goto 150
c     Removed - 31oct06

c     step change with equal (constant) sloped segments
      if(mknt1 + mknt2 .ge. minslptest) then
        if(lsqopt .eq. 1) then
          cmodel = 'FITTPR1'
          call fittpr1(nobs, yseg, ichgpt, rmutpr, 
     *      slptpr1, SSEful, mknt1, mknt2, ierr)
          if(ierr .eq. 0) then
            mKnt = mknt1 + mknt2
            F1n3 = ((SSEred - SSEful)/1.)/(SSEful/(mknt-3))
            ftpr = critval(mknt-3, 4)
            call bayes(mknt, SSEful, 4, qtpr, rsq1, rsq2)
            y1 = rmutpr(1) + slptpr1 * xseg(ichgpt)
            y2 = rmutpr(2) + slptpr1 * xseg(beg2)
            est = y1 - y2
            if(ibdebug .gt. 1) 
     *        write(6,'(a,5f8.2,f8.3," -------",3f8.2,2i5)') 
     *          cmodel, qtpr, rsq1, rsq2, rmutpr, slptpr1, F1n3, 
     *          ftpr, est, mknt1, mknt2
          else
            if(ibdebug .gt. 1) 
     *        print *,cmodel,' - Unable to fit model - skipping'
          endif  
        endif
      endif    

      if(mkntme1 + mkntme2 .ge. minslptest) then
        cmodel = 'KTHTPR1'
        call kthtpr1(nobs, yseg, ichgpt, rmutpr, 
     *    slptpr1, SSEful, mkntme1, mkntme2, ierr)
        if(ierr .eq. 0) then
          mKnt = mkntme1 + mkntme2
          F1n3 = ((SSEredmed - SSEful)/1.)/(SSEful/(mknt-3))
          ftpr = critval(mknt-3, 4)
          call bayes(mknt, SSEful, 4, qtpr, rsq1, rsq2)
          y1 = rmutpr(1) + slptpr1 * xseg(ichgpt)
          y2 = rmutpr(2) + slptpr1 * xseg(beg2)
          est = y1 - y2
          if(ibdebug .gt. 1) 
     *      write(6,'(a,5f8.2,f8.3," -------",3f8.2,2i5)') 
     *        cmodel, qtpr, rsq1, rsq2, rmutpr, slptpr1, F1n3, 
     *        ftpr, est, mkntme1, mkntme2
          if(qtpr .lt. qmin) then
            qmin = qtpr
            iqtype = 4
            rmuq(1) = rmutpr(1)
            rmuq(2) = rmutpr(2)
            slpq(1) = slptpr1
            slpq(2) = slpq(1)
            sseq = sqrt(SSEful/(mkntme1+mkntme2))
            mkfinq(1) = mkntme1
            mkfinq(2) = mkntme2
            qstat = F1n3
            qcrit = ftpr
            qoff = est
            if((istest(iopt).eq.0 .or. abs(slptpr1).lt.slpmthmax) .and.
     *        (ictest(iopt) .eq. 0 .or. F1n3 .ge. ftpr)) then
              tmin = qtpr
              ittype = 4
              rmut(1) = rmutpr(1)
              rmut(2) = rmutpr(2)
              slpt(1) = slptpr1
              slpt(2) = slpq(1)
              sset = sqrt(SSEful/(mkntme1+mkntme2))
              mkfint(1) = mkntme1
              mkfint(2) = mkntme2
              tstat = F1n3
              tcrit = ftpr
              toff = est
            endif
          endif
        else
          if(ibdebug .gt. 1) 
     *      print *,cmodel,' - Unable to fit model - skipping'
        endif
      else
        if(ibdebug .gt. 1) 
     *    print *,cmodel,' - Unable to fit model - skipping'
      endif  
      
c     step change with any sloped segments (full two phase regression)
  150 if(lsqopt .eq. 1) then
        cmodel = 'FITTPR2'
        if(mknt1 .ge. minslptest .and. mknt2 .ge. minslptest) then
          y1 = yint1 + slope1 * xseg(ichgpt)
          y2 = yint2 + slope2 * xseg(beg2)
          est = y1 - y2
          mKnt = mknt1 + mknt2
          SSEful = SSEslope1 + SSEslope2
          F2n4 = ((SSEred - SSEful)/2.)/(SSEful/(mknt-4))
          ftpr = critval(mknt-4, 5)
          call bayes(mknt, SSEful, 5, qtpr, rsq1, rsq2)
          if(ibdebug .gt. 1) 
     *      write(6,'(a,5f8.2,2f8.3,3f8.2,2i5)') cmodel, qtpr, rsq1,
     *      rsq2, yint1, yint2, slope1, slope2, F2n4, ftpr, est,
     *      mknt1, mknt2
        else
          if(ibdebug .gt. 1) 
     *      print *,cmodel,' - Unable to fit model - skipping'
        endif  
      endif  

      cmodel = 'KTHTPR2'
      if(mkntme1 .ge. minslptest .and. mkntme2 .ge. minslptest) then
        y1 = yintme1 + slpme1 * xseg(ichgpt)
        y2 = yintme2 + slpme2 * xseg(beg2)
        est = y1 - y2
        mKnt = mkntme1 + mkntme2
        SSEful = SSEslpme1 + SSEslpme2
        F2n4 = ((SSEredmed - SSEful)/2.)/(SSEful/(mknt-4))
        ftpr = critval(mknt-4, 5)
        call bayes(mknt, SSEful, 5, qtpr, rsq1, rsq2)
        if(ibdebug .gt. 1) 
     *    write(6,'(a,5f8.2,2f8.3,3f8.2,2i5)') cmodel, qtpr, rsq1, rsq2,
     *    yintme1, yintme2, slpme1, slpme2, F2n4, ftpr, est, mkntme1, 
     *    mkntme2
        if(qtpr .lt. qmin) then
          qmin = qtpr
          iqtype = 5
          rmuq(1) = yintme1
          rmuq(2) = yintme2
          slpq(1) = slpme1
          slpq(2) = slpme2
          sseq = sqrt(SSEful/(mkntme1+mkntme2))
          mkfinq(1) = mkntme1
          mkfinq(2) = mkntme2
          qstat = F2n4
          qcrit = ftpr
          qoff = est
          idfree = 5
          if((istest(iopt) .eq. 0 .or.(abs(slpme1) .lt. slpmthmax
     *      .and. abs(slpme2) .lt. slpmthmax)) .and.
     *      (ictest(iopt) .eq. 0 .or. F2n4 .ge. ftpr)) then
            tmin = qtpr
            ittype = 5
            rmut(1) = yintme1
            rmut(2) = yintme2
            slpt(1) = slpme1
            slpt(2) = slpme2
            sset = sqrt(SSEful/(mkntme1+mkntme2))
            mkfint(1) = mkntme1
            mkfint(2) = mkntme2
            tstat = F2n4
            tcrit = ftpr
            toff = est
          endif
        endif
      else
        if(ibdebug .gt. 1) 
     *    print *,cmodel,' - Unable to fit model - skipping'
      endif  

c     step change with flat-to-sloped segments
  160 if(lsqopt .eq. 1) then
        cmodel = 'FITTPR3'
        if(mknt1 .ge. mintest .and. mknt2 .ge. minslptest) then
          y1 = ymean1
          y2 = yint2 + slope2 * xseg(beg2)
          est = y1 - y2
          mKnt = mknt1 + mknt2
          SSEful = SSEflat1 + SSEslope2
          F1n3 = ((SSEred - SSEful)/1.)/(SSEful/(mknt-3))
          ftpr = critval(mknt, 4)
          call bayes(mKnt, SSEful, 4, qtpr, rsq1, rsq2)
          if(ibdebug .gt. 1) 
     *      write(6,'(a,5f8.2," -------",f8.3,3f8.2,2i5)') 
     *        cmodel, qtpr, rsq1, rsq2, ymean1, yint2, slope2, F1n3,
     *        ftpr, est, mknt1, mknt2
        else
          if(ibdebug .gt. 1) 
     *      print *,cmodel,' - Unable to fit model - skipping'
        endif  
      endif

      cmodel = 'KTHTPR3'
      if(mkntme1 .ge. mintest .and. mkntme2 .ge. minslptest) then
        y1 = ymed1
        y2 = yintme2 + slpme2 * xseg(beg2)
        est = y1 - y2
        mKnt = mkntme1 + mkntme2
        SSEful = SSEflatme1 + SSEslpme2
        F1n3 = ((SSEredmed - SSEful)/1.)/(SSEful/(mknt-3))
        ftpr = critval(mknt, 4)
        call bayes(mKnt, SSEful, 4, qtpr, rsq1, rsq2)
        if(ibdebug .gt. 1) 
     *    write(6,'(a,5f8.2," -------",f8.3,3f8.2,2i5)') 
     *      cmodel, qtpr, rsq1, rsq2, ymed1, yintme2, slpme2, F1n3,
     *      ftpr, est, mkntme1, mkntme2
        if(qtpr .lt. qmin) then
          qmin = qtpr
          iqtype = 6
          rmuq(1) = ymed1
          rmuq(2) = yintme2
          slpq(1) = 0.0
          slpq(2) = slpme2
          sseq = sqrt(SSEful/(mkntme1+mkntme2))
          mkfinq(1) = mkntme1
          mkfinq(2) = mkntme2
          qstat = F1n3
          qcrit = ftpr
          qoff = est
          idfree = 4
          if((istest(iopt).eq.0 .or. abs(slpme2).lt.slpmthmax) .and.
     *      (ictest(iopt) .eq. 0 .or. F1n3 .ge. ftpr)) then
            tmin = qtpr
            ittype = 6
            rmut(1) = ymed1
            rmut(2) = yintme2
            slpt(1) = 0.0
            slpt(2) = slpme2
            sset = sqrt(SSEful/(mkntme1+mkntme2))
            mkfint(1) = mkntme1
            mkfint(2) = mkntme2
            tstat = F1n3
            tcrit = ftpr
            toff = est
          endif
        endif
      else
        if(ibdebug .gt. 1) 
     *    print *,cmodel,' - Unable to fit model - skipping'
      endif  

c     step change with sloped-to-flat segments
 170  if(lsqopt .eq. 1) then
        cmodel = 'FITTPR4'
        if(mknt1 .ge. minslptest .and. mknt2 .ge. mintest) then
          y1 = yint1 + slope1 * xseg(ichgpt)
          y2 = ymean2
          est = y1 - y2
          mKnt = mknt1 + mknt2
          SSEful = SSEslope1 + SSEflat2
          F1n3 = ((SSEred - SSEful)/1.)/(SSEful/(mknt-3))
          ftpr = critval(mknt, 4)
          call bayes(mKnt, SSEful, 4, qtpr, rsq1, rsq2)
          if(ibdebug .gt. 1) 
     *      write(6,'(a,5f8.2,f8.3," -------",3f8.2,2i5)') 
     *        cmodel, qtpr, rsq1, rsq2, yint1, ymean2, slope1, F1n3,
     *        ftpr, est, mknt1, mknt2
        else
          if(ibdebug .gt. 1) 
     *      print *,cmodel,' - Unable to fit model - skipping'
        endif  
      endif  

      cmodel = 'KTHTPR4'
      if(mkntme1 .ge. minslptest .and. mkntme2 .ge. mintest) then
        y1 = yintme1 + slpme1 * xseg(ichgpt)
        y2 = ymed2
        est = y1 - y2
        mKnt = mkntme1 + mkntme2
        SSEful = SSEslpme1 + SSEflatme2
        F1n3 = ((SSEredmed - SSEful)/1.)/(SSEful/(mknt-3))
        ftpr = critval(mknt, 4)
        call bayes(mKnt, SSEful, 4, qtpr, rsq1, rsq2)
        if(ibdebug .gt. 1) 
     *    write(6,'(a,5f8.2,f8.3," -------",3f8.2,2i5)') 
     *      cmodel, qtpr, rsq1, rsq2, yintme1, ymed2, slpme1, F1n3, 
     *      ftpr, est, mkntme1, mkntme2
        if(qtpr .lt. qmin) then
          qmin = qtpr
          iqtype = 7
          rmuq(1) = yintme1
          rmuq(2) = ymed2
          slpq(1) = slpme1
          slpq(2) = 0.0
          sseq = sqrt(SSEful/(mkntme1+mkntme2))
          mkfinq(1) = mkntme1
          mkfinq(2) = mkntme2
          qstat = F1n3
          qcrit = ftpr
          qoff = est
          idfree = 4
          if((istest(iopt).eq.0 .or. abs(slpme1).lt.slpmthmax) .and.
     *      (ictest(iopt) .eq. 0 .or. F1n3 .ge. ftpr)) then
            tmin = qtpr
            ittype = 7
            rmut(1) = yintme1
            rmut(2) = ymed2
            slpt(1) = slpme1
            slpt(2) = 0.0
            sset = sqrt(SSEful/(mkntme1+mkntme2))
            mkfint(1) = mkntme1
            mkfint(2) = mkntme2
            tstat = F1n3
            tcrit = ftpr
            toff = est
          endif
        endif
      else
        if(ibdebug .gt. 1) 
     *    print *,cmodel,' - Unable to fit model - skipping'
      endif
      
c     if the segment intercept test is on: recalc chgpt loc for loop-back
  180 if(iopt .ne. 3) then
c       no segment intercept test requested - finish 
        goto 200
      else
c       must be correct type
        if(iqtype .lt. 5) goto 200
        if(isloop .eq. 0) then
c         save min Qval parms on first pass
          ichgloop = ichgpt
          rmuloop1 = rmuq(1)
          rmuloop2 = rmuq(2)
          slploop1 = slpq(1)
          slploop2 = slpq(2)
          mkloop1 = mkntme1
          mkloop2 = mkntme2
          SSEfloop = SSEslpme1 + SSEslpme2
          idloop = idfree
        endif
c       recalculate x intercept for next iteration from this one
        i0 = -1 * int(((rmuloop2 - rmuloop1)/(slploop2-slploop1))+.5)
        ix0 = i0 - imo1 + 1
        if(ibdebug .gt. 1) 
     *    print *,' est,estloop,imo1,ix0: ', qoff,estloop,imo1,ix0
c       if orig. is close enough to brkpt est - use as brkpt
        if(abs(ix0 - ichgpt) .le. 2 .and. isloop .eq. 0) goto 190
c       zero intercept must be within series range and
c       no need to recalculate if the new chgpt == old chgpt
        if(ix0 .lt. mintest .or. ix0 .gt. nobs-mintest .or.
     *    ix0 .eq. ichgpt)then
c              orig. lin est outside of range - loop out
          if(isloop .eq. 0) goto 200
c         loop lin est outside of range or == last - test last loop est
          goto 190
        endif  
c       changepoint must be on a non-missing value 
        do i = ix0, 1, -1
          if(yseg(i) .ne. amiss) goto 185
        enddo
c       the following should never happen
        print *,' missing lower segment?!?!?'
        goto 200  
c       save best Qval offset and new changepoint location for looping
  185   isloop = 1
        estloop = qoff
        ichgpt = i
c       loop back and recalculate
        goto 100
      endif      

c     In the loop-back test and the breakpoint position (i.e. minimum offset)
c       has been found. Test against the original minimum Qval case
  190 if(mkloop1 .ge. mintest .and. mkloop2 .ge. mintest) then
        mKnt = mkloop1 + mkloop2
        Fxnx=((SSEredmed - SSEfloop)/(idloop-2))/
     *    (SSEfloop/(mknt-idloop+1))
        ftpr = critval(mknt, idloop)
        call bayes(mknt, SSEfloop, idloop, qtpr, rsq1, rsq2)
        if(ibdebug .gt. 1) 
     *    write(6,'(a,5f8.2,2f8.3,3f8.2,2i5)') cmodel, 
     *    qtpr, rsq1, rsq2,rmuloop1,rmuloop2,slploop1,slploop2, Fxnx,
     *    ftpr, est, mkloop1, mkloop2
        if(qtpr .lt. qmin .and. abs(slploop1) .lt. slpmthmax .and.
     *    abs(slploop2) .lt. slpmthmax)then
          iqtype = 8
          rmuq(1) = rmuloop1
          rmuq(2) = rmuloop2
          slpq(1) = slploop1
          slpq(2) = slploop2
          sseq = sqrt(SSEfloop/(mkloop1+mkloop2-1))
          mkfinq(1) = mkloop1
          mkfinq(2) = mkloop2
          qstat = Fxnx
          qcrit = ftpr
          qoff = est
          
c         the breakpoint Qval is less than the original - replace chgpt
c         back out the NEW iend1 from the compressed array
          ichgpt = ichgloop
          lend1 = end1
          do i = 1, end2
            if(i .gt. ichgpt) goto 195
            if(yseg(i) .ne. amiss) end1 = i
          enddo
  195     if(ibdebug .gt. 1)           
     *      print *,' Replaced ',ichgpt0,lend1,' with SIT ',ichgpt,end1
        else 
c         if adjusted intercept NOT successful in loop back - leave chgpt
          ichgpt = ichgpt0
        endif  
      else 
        if(ibdebug .gt. 1) 
     *    print *,cmodel,' - Unable to fit model - skipping'
      endif

c     calculate the offset for the best model      
c     if istraight comes out 1 then use the min SLR model
c        else use the incoming model and calc offset
  200 if(ibdebug .ge. 3) then
        write(6,'(" Post: ", 2i2, 9f9.2, 2i5)') iopt,
     *    iqtype, qmin, sseq, rmuq, slpq, qstat, qcrit, qoff, mkfinq
c       Print out "coincident" minbic output ONLY if different from "post" 
c        if(ittype .ne. iqtype) 
c     *   write(6,'(" Coin: ", 2i2, 9f9.2, 2i5)') iopt,
c     *    ittype, tmin, sset, rmut, slpt, tstat, tcrit, toff, mkfint
      endif  
      istraight = 0
      if( iqtype .lt. 3) then
c       if iqtype == 1 or 2, there is NO MOVE!
        istraight = 1
      else if(iopt .eq. 2 .or. iopt .eq. 6) then
c       send back the TPR0 fit for any NON-SLR best model      
        rmuq(1) = ymed1
        rmuq(2) = ymed2
        slpq(1) = 0.0
        slpq(2) = 0.0
      else if(iopt .eq. 4) then
c       Used for offset amplitude estimate - check for realistic slope
        if(abs(slpq(1)) .ge. slpmthmax .or.
     *    abs(slpq(2)) .ge. slpmthmax) then
          if(ibdebug .gt. 1) 
     *      write(6,'("Iopt:", i2, " Fail slpchk:",2f7.3,2i5)') 
     *        iopt, slpq(1), slpq(2), mkfinq(1), mkfinq(2)
          if(abs(slpq(1)) .ge. slpmthmax .and.
     *      abs(slpq(2)) .lt. slpmthmax) then
            ifail = 1
          else if(abs(slpq(1)) .lt. slpmthmax .and.
     *      abs(slpq(2)) .ge. slpmthmax) then
            ifail = 2
          else if(abs(slpq(1)) .ge. slpmthmax .and.
     *      abs(slpq(2)) .ge. slpmthmax) then
            ifail = 3
          endif       
        endif  
      else if(iopt .eq. 5) then
c       used to test unsupported metadata chgpts - check both slope and critvar
        if(qstat .lt. qcrit) then
          if(ibdebug .gt. 1) 
     *      write(6,'("Iopt:", i2, " Fail critvar: ",2f7.2,2i5)') 
     *        iopt, qstat, qcrit, mkfinq(1), mkfinq(2)
          istraight = 1
        endif    
        if( abs(slpq(1)) .ge. slpmthmax .or.
     *    abs(slpq(2)) .ge. slpmthmax) then
          if(ibdebug .gt. 1) 
     *      write(6,'("Iopt:", i2, " Fail slpchk:",2f7.3,2i5)') 
     *        iopt, slpq(1), slpq(2), mkfinq(1), mkfinq(2)
        endif  
      else
c       used for model vs stat test comparison and interim chgpt tests
c        check the critical variance threshold (t- or f-test) of proper model
        if(qstat .lt. qcrit) then
          if(ibdebug .gt. 1) 
     *      write(6,'("Iopt:", i2, " Fail critvar: ",2f7.2,2i5)') 
     *        iopt, qstat, qcrit, mkfinq(1), mkfinq(2)
          istraight = 1
        endif    
      endif
      if(istraight .eq. 1) then
c       either already a straight line or failed the iopt defined test
c       use the best SLR model
        qoff = 0.0
        mkfinq(1) = mknt
        mkfinq(2) = 0
        qstat = 0.0
        qcrit = 99.0
        if(qslr0 .lt. qslr1) then
          qmin = qslr0
          iqtype = 1
          rmuq(1) = ymed
          rmuq(2) = ymed
          slpq(1) = 0.0
          slpq(2) = 0.0
          sseq = sqrt(sseflat/mknt)
        else
          qmin = qslr1
          iqtype = 2
          rmuq(1) = Yintmed
          rmuq(2) = rmuq(1)
          slpq(1) = slpmed
          slpq(2) = slpq(1)
          sseq = sqrt(SSEredmed/mknt)
        endif
      else
c       did not fail whichever test for the iopt input -
c         keep incoming model and calculate offset
        y1 = rmuq(1) + slpq(1) * xseg(ichgpt)
        y2 = rmuq(2) + slpq(2) * xseg(beg2)
        qoff = y1 - y2        
      endif
      
c     when all is said and done, which test result is passed back?
c     icoin = 0 : return the post threshold results
c     icoin = 1 : return the coincident threshold results
      if(icoin .eq. 1) then
        iqtype = ittype
        rmuq(1) = rmut(1)
        rmuq(2) = rmut(2)
        slpq(1) = slpt(1)
        slpq(2) = slpt(2)
        sseq = sset
        mkfinq(1) = mkfint(1)
        mkfinq(2) = mkfint(2)
        qstat = tstat
        qcrit = tcrit
        qoff = toff
      endif        

c     redefine slpq to the kthline solution no matter what MINBIC model             
      slpq(1) = slpme1
      slpq(2) = slpme2

      return
      end

c ===================================================================
c      THE FOLLOWING SET OF ROUTINES FIT A GIVEN MODEL (WITH A CHANGEPOINT
C        AS NEEDED) AND RETURN A FIT AND RESIDUAL SERIES FOR EACH. ALSO
C        RETURNED ARE THE INTERCEPT <MU> AND SLOPE <ALPHA> FOR THE SEGMENT(S)
c ===================================================================

c ===================================================================

      subroutine fittpr1(nObs, rSeries, iChgPnt,
     &  rMu, rAlpha, rSumSqrT, mKnt1, mknt2, ierr)
!-----------------------------------------------------------------------
! Subroutine to calculate two phase regression for series with same
!   slope in both segments
!
!                          -on input-
!  nObs - sample size
!  rSeries - dependent data (predictand)
!  iChgPnt - last position in first segment
!  rMiss - Missing value code (assumed to be a large negative number) 
!
!                          -on output-
!  rResid - array of least squares residuals
!  rMu - y-intercept of regression line for each segment
!  rAlpha - slopes of each segment
!  rT - t-statistic for the slope of the straight line fit
!  rSumSqrE - Sum standard error of residual series
!  rFIt - least squares fit of time series rYdata
!
!-----------------------------------------------------------------------      
      INCLUDE 'inhomog.parm.mthly.incl' 
      INCLUDE 'inhomog.MDparm.mthly.incl'

      real rSeries(kmo)
      real rTmpResid(kmo)
      real rResid(kmo)
      real rMu(2)
      real rSumSqrE(2)
      
      rmiss = amiss
      ierr = 0
      
!c  calculate the residuals for the "full model" with constant slope
        rMn1 = 0.0
        rMn2 = 0.0
        rMn = 0.0
        rTmn1 = 0.0
        rTmn2 = 0.0
        iN1 = 0
        iN2 = 0
        iN = 0
        rAlpha = 0.0
        rMu(1) = 0.0
        rMu(2) = 0.0
        rNum = 0.0
        rDen = 0.0
        rSumSqrT = 0.0
        rSumSqrE(1) = 0
        rSumSqrE(2) = 0
        rSumSqrX = 0
        rSb = 0.0
        rT = 0.0
        
        do j = 1, nobs
          rResid(j) = rMiss
        enddo   

        do j = 1,iChgPnt
          if(rSeries(j) .gt. rMiss + 1) then
            rMn1 = rMn1 + rSeries(j)
            rTmn1 = rTmn1 + j
            iN1 = iN1 + 1
          endif
        enddo
        
        if(iN1 .ge. mintest) then     
          rMn1 = rMn1/iN1
          rTmnt = rTmn1
          rTMn1 = rTMn1/iN1
        else
          ierr = 1
          return
        endif 
        
        do j=1,iChgPnt
          if(rSeries(j) .gt. rMiss + 1) then
             rTmpResid(j) = rSeries(j) - rMn1
          else
             rTmpResid(j) = rMiss
          endif
        enddo
        
!c      Resids for first segment ready
        
!c      Calculate mean for second segment
        do j = iChgPnt+1, nObs
          if(rSeries(j) .gt. rMiss) then 
            rMn2 = rMn2 + rSeries(j)
            rTmn2 = rTmn2 + j
            iN2 = iN2 + 1
          endif
        enddo
        
        if(iN2 .ge. mintest) then
          rMn2 = rMn2/iN2
          rTmnt = rTmnt + rTmn2
          rTmn2 = rTmn2/iN2
        else
          ierr = 1
          return
        endif
        
c       generate the center of both segments
        iNt = iN1 + iN2
        rTmnt = rTmnt/iNt
        
!c      Resids for second segment ready 
        do j=iChgPnt+1,nObs
          if(rSeries(j) .gt. rMiss) then 
            rTmpResid(j) = rSeries(j) - rMn2
          else
            rTmpResid(j) = rMiss
          endif
        enddo
        
!       now calculate common slope for two segments
        do j=1,iChgPnt
          if(rTmpResid(j) .gt. rMiss) then
            rNum = rNum + (j-rTmn1)*(rTmpResid(j))
            rDen = rDen + (j-rTmn1)*(j-rTmn1)
          endif
        enddo
        
        do j=iChgPnt+1,nObs
          if(rTmpResid(j) .gt. rMiss) then
            rNum = rNum + (j-rTmn2)*(rTmpResid(j))
            rDen = rDen + (j-rTmn2)*(j-rTmn2)
          endif
        enddo
        
        rAlpha = rNum / rDen
      
!       now calculate the residuals
        rMu(1) = rMn1 - rAlpha * rTmn1
        rMu(2) = rMn2 - rAlpha * rTmn2

        do j = 1,iChgPnt
          if(rSeries(j) .gt. rMiss + 1) then        
            rResid(j) = rSeries(j) - rMu(1) - rAlpha*j
            rSumSqrE(1) = rSumSqrE(1) + rResid(j)**2
            rSumSqrX = rSumSqrX + (real(j) - rTmnt)**2
           endif
        enddo
          
        do j = iChgPnt+1,nObs
          if(rSeries(j) .gt. rMiss + 1) then
            rResid(j) = rSeries(j) - rMu(2) - rAlpha*j
            rSumSqrE(2) = rSumSqrE(2) + rResid(j)**2
            rSumSqrX = rSumSqrX + (real(j) - rTmnt)**2
          endif
        enddo
        
        rSumSqrT = rSumSqrE(1) + rSumSqrE(2)

        rAmpStep = rMu(2)-rMu(1)
        mKnt1 = iN1
        mKnt2 = iN2
      return
      end 

c ===================================================================

      subroutine kthtpr1(nObs, rSeries, iChgPnt,
     &  rMu, rAlpha, rSumSqrT, mknt1, mKnt2, ierr)
!-----------------------------------------------------------------------
! Subroutine to calculate two phase regression for series with single
!   slope in both segments
!
!                          -on input-
!  nObs - sample size
!  rSeries - dependent data (predictand)
!  iChgPnt - last position in first segment
!  rMiss - Missing value code (assumed to be a large negative number) 
!
!                          -on output-
!  rResid - array of least squares residuals
!  rMu - y-intercept of regression line for each segment
!  rAlpha - slopes of each segment
!  rT - t-statistic for the slope of the straight line fit
!  rSumSqrE - Sum standard error of residual series
!  rFIt - least squares fit of time series rYdata
!
!-----------------------------------------------------------------------      
      INCLUDE 'inhomog.parm.mthly.incl' 
      INCLUDE 'inhomog.MDparm.mthly.incl'

      real rSeries(kmo), rSeg(kmo)
      real rResid(kmo)
      real rMu(2)
      real rSumSqrE(2)
      
      rmiss = amiss
      ierr = 0
      
!c  calculate the residuals for the "full model" with constant slope
        iN1 = 0
        iN2 = 0
        iNt = 0
        rAlpha = 0.0
        rMu(1) = 0.0
        rMu(2) = 0.0
        rNum = 0.0
        rDen = 0.0
        rSumSqrT = 0.0
        rSumSqrE(1) = 0
        rSumSqrE(2) = 0
        rSumSqrX = 0
        
        do j = 1, nobs
          rResid(j) = rMiss
        enddo   

!c      Calculate mean for first and second segment together
        sX = 0.0
        do j = 1,iChgPnt
          if(rSeries(j) .gt. rMiss+1) then
            iN1 = iN1 + 1
            rSeg(iN1) = rSeries(j)
            sX = sX + j
          endif
        enddo
        
        do j = iChgPnt+1, nObs
          if(rSeries(j) .gt. rMiss+1) then 
            iN2 = iN2 + 1
            rSeg(iN1 + iN2) = rSeries(j)
            sX = sX + j
          endif
        enddo
        
        iNt = iN1 + iN2
        rTmnt = sX / iNt
        
c       use kendall-theil method with internal breakpoint and single slope
        call kendalltheil(iNt, rSeg, rMu, rAlpha, 
     *    SSEred, mknt, iN1+1, ierr)
        if(ierr .eq. 1) return

        do j = 1,iChgPnt
          if(rSeries(j) .gt. rMiss+1) then        
            rResid(j) = rSeries(j) - rMu(1) - rAlpha*j
            rSumSqrE(1) = rSumSqrE(1) + rResid(j)**2
            rSumSqrX = rSumSqrX + (real(j) - rTmnt)**2
           endif
        enddo
          
        do j = iChgPnt+1,nObs
          if(rSeries(j) .gt. rMiss+1) then
            rResid(j) = rSeries(j) - rMu(2) - rAlpha*j
            rSumSqrE(2) = rSumSqrE(2) + rResid(j)**2
            rSumSqrX = rSumSqrX + (real(j) - rTmnt)**2
          endif
        enddo
        
        rSumSqrT = rSumSqrE(1) + rSumSqrE(2)

        rAmpStep = rMu(2)-rMu(1)
        mKnt1 = iN1
        mKnt2 = iN2
      return
      end 

c ===================================================================

      real function critval(nobs, iqtype)
      real alpha
      integer nobs, iqtype
c
c     set critical values for t and f tests 

c     this set is at 95% confidence
c      real tval(34)/ 6.31, 2.92, 2.35, 2.13, 2.02, 1.94, 1.90,
c     *   1.86, 1.83, 1.81, 1.80, 1.78, 1.77, 1.76, 1.75, 1.75, 1.74,
c     *   1.73, 1.73, 1.73, 1.72, 1.72, 1.71, 1.71, 1.71, 1.71, 1.70,
c     *   1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.70/
c     this set is at 95% confidence (two-tailed)
      real tval(34)/ 12.71, 4.30, 3.18, 2.78, 2.57, 2.45, 2.37,2.31,
     *   2.26, 2.23, 2.20, 2.18, 2.16, 2.15, 2.13, 2.12, 2.11, 2.10,
     *   2.09, 2.09, 2.08, 2.07, 2.07, 2.06, 2.06, 2.06, 2.05, 2.05,
     *   2.05, 2.04, 2.02, 2.00, 1.98, 1.96/
c     this is the F1 column
      real f1val(34)/ 161.4, 18.50, 10.1,  7.7,  6.6,  6.0,  5.6,
     *    5.3,  5.1,  5.0,  4.8,  4.7,  4.7,  4.6,  4.5,  4.5,  4.5,
     *    4.4,  4.4,  4.4,  4.3,  4.3,  4.3,  4.3,  4.2,  4.2,  4.2,
     *    4.2,  4.2,  4.2,  4.1,  4.0,  3.9,  3.8/
c     this is the F2 column
      real f2val(34)/ 199.5, 19.00, 9.55, 6.94, 5.79, 5.14, 4.74,
     *   4.46, 4.26, 4.10, 3.98, 3.89, 3.81, 3.74, 3.68, 3.63, 3.59, 
     *   3.55, 3.52, 3.49, 3.47, 3.44, 3.42, 3.40, 3.39, 3.37, 3.35, 
     *   3.34, 3.33, 3.32, 3.23, 3.15, 3.07, 3.00/
c     this is the F3 column
      real f3val(34)/ 215.7, 19.16, 9.28, 6.59, 5.41, 4.76, 4.35,
     *   4.07, 3.86, 3.71, 3.59, 3.49, 3.41, 3.34, 3.29, 3.24, 3.20, 
     *   3.16, 3.13, 3.10, 3.07, 3.05, 3.03, 3.01, 2.99, 2.98, 2.96, 
     *   2.95, 2.93, 2.92, 2.84, 2.76, 2.68, 2.60/


c     threshold lengths for t & f values
      integer length(34) / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
     *   14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
     *   29, 30, 40, 60, 120, 9999/ 

      do n = 33,1,-1
        if(nobs .ge. length(n)) go to 10
      end do 

c     for iqtype = 3 (TPR0) return the t-test stat
   10 if(iqtype .eq. 3) then
        val1 = tval(n)
        val2 = tval(n+1)
      else if(iqtype .eq. 4) then
c       for the iqtype = 4 (TPR1 & TPR3 & TPR4) return the F1 Stat
        val1 = f1val(n)
        val2 = f1val(n+1)
      else if(iqtype .eq. 5) then
c       for the iqtype = 5 (TPR2) return the F2 stat
        val1 = f2val(n)
        val2 = f2val(n+1)
      else
        print *,' Unknown iqtype in alpha'
        stop
      endif        
        
c     one-to-one correspondense at low N
      if(n .le. 30) then
        critval = val1
      else  
c       use ratio for larger N
        lseg = length(n+1) - length(n)
        nseg = nobs - length(n)
        critval = ((val2-val1)/lseg)*nseg + val1
      endif  
   
      return
      end

c     =======================================================================

      subroutine bayes(iN,rSumSqr,irq,rsq,rsq1,rsq2)
!-----------------------------------------------------------------------
!
!  This function calculates the Schwarz Bayesian Information Criterian
!  statistic for a series and its model fit.  The calculation is based  
!  on Schwarz (1978), Ann. Stat., 6, 461-464 and Seidel and Lanzante
!  (2004), JGR-Atmos, 109, D14108.
! 
!                      -on input-
!  iN - number of non-missing values
!  rSumSqr - Sum Sqaure of residuals of model fit to rSeries
!  rq - dimensions (degrees of freedom) of model required to fit rSeries
! 
!                     -on output-
!  rSq - Schwarz statistic (the lower the value, the better the model
!        fit
!-----------------------------------------------------------------------    
      real rq
      integer iN

c      log10 - Schwarz
      rsq1 = iN * log10(rSumSqr/real(iN))
      rsq2 = irq * log10(real(iN))
c     natural log - Sawa
c      rsq1 = iN * log(rSumSqr/real(iN))
c      rsq2 = irq * log(real(iN))
      rSq = rsq1 + rsq2
      
      return
      end

c     =======================================================================
c     
c     =======================================================================

      subroutine snits (n,z,ts,mcnt)
c**********************************************************************
c
c     Subroutine snits computes the standard normal inhomogeneity test
c     statistic for a step change or shift in the standardized
c     difference data between a candidate and the average of its
c     neighbors.  The algorithm to calculate ts is Eq.(4) in Alexander-
c     sson and Moberg, Int'l Jour. Climat., 17, 25-34 (1997).  
c
c     23 September, 1999; 18 October, 1999; 9 December, 1999
c     6 July, 2001; November 2001; April 2002--to include date intervals
c
c
c                               -on input-
c
c       n       The number of months for z & ts arrays
c
c       mcnt    the number of non-missing months or years
c               in the time series (n and mcnt can be combined
c               if no missing months are passed to the routine)
c
c       z       The n by 1 z-series of standardized monthly
c               temperature differences produced by sbqzoft
c
c                               -on output-
c
c       ts      The mcnt series of the snit test statistic
c               ts for a shift or step-change in monthly mean
c               temperature at the candidate station
c       
c***********************************************************************
      INCLUDE 'inhomog.parm.mthly.incl' 

      dimension z(n),ts(n)

      do i = 1, n
        ts(i) = amiss
      enddo

c     Calculate the likelihood ratio test statistic ts
      iCount = 1
      do i = 1, mcnt - 1    !use this if window is used
        if(z(i) .gt. amiss + 1.0) then
          zmn1 = 0.0
          zmn2 = 0.0
          n1 = 0
          n2 = 0
          do j = 1,i
            if(z(j) .gt. amiss + 1.0) then
              zmn1 = zmn1 + z(j)
              n1 = n1 + 1
            endif  
          enddo
          if(n1 .ne. 0) then     
            zmn1 = zmn1/n1
          else
            goto 10
          endif
          do j = i+1, mcnt
            if(z(j) .gt. amiss + 1.0) then
              zmn2 = zmn2 + z(j)
              n2 = n2 + 1
            endif  
          enddo
          if(n2 .ne. 0) then
            zmn2 = zmn2/n2
          else
            goto 10
          endif
          ts(i) = n1*zmn1*zmn1 + n2*zmn2*zmn2
          iCount=iCount+1
        endif    
   10 enddo
      return
      end

c     =======================================================================
      
      subroutine lookup (num,t90,t95,t975,iopt)
c**********************************************************************
c
c     Subroutine lookup computes the critical values for the
c     standard normal inhomogeneity test at the 90%, 95%, and 97.5%
c     confidence levels (one-sided).  The critical values are
c     dependent on the number of data and are the same regardless
c     of whether a step-change (shift) or trend is being tested.
c     19 October, 1999    6 July, 2001
c
c                               -on input-
c
c       num     The number of months in the z-series
c
c                               -on output-
c
c       t90     The critical values of t at the 90%, 95%, and 97.5%
c       t95     one-sided confidence levels, respectively.
c       t975
c
c***********************************************************************
      INCLUDE 'inhomog.parm.mthly.incl' 
      INCLUDE 'inhomog.comm.mthly.incl'
      parameter (numsig = 14)
      parameter (numtpr = 96)

      dimension sig90(numsig),sig95(numsig),sig975(numsig)
      integer nsig(numsig)
c     two phase regression - used with iopt = 2
      real nTPR(numtpr),rTPR95(numtpr),rTPR975(numtpr),rTPR99(numtpr)

c     Set look-up tables of critical values of ts for one-sided
c     90%, 95%, and 97.5% significance levels.
c     The values   are taken from the 1986 and 1997 papers by Alexandersson 
c     and Alexandersson and Moberg, respectively.
      data nTPR/   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
     *  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
     *  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
     *  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,
     *  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,
     *  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,
     *  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,
     *  94,  95,  96,  97,  98,  99, 100/
      data rTPR95/   792.38, 55.93, 24.17, 16.59, 13.50, 11.69, 10.52,
     *   9.76,  9.20,  8.80,  8.51,  8.25,  8.09,  7.97,  7.81,  7.69,
     *   7.65,  7.58,  7.49,  7.38,  7.33,  7.32,  7.26,  7.23,  7.15,
     *   7.17,  7.20,  7.17,  7.14,  7.04,  7.07,  7.06,  7.04,  7.01,
     *   7.03,  7.04,  7.04,  6.97,  6.99,  7.00,  6.93,  6.99,  6.96,
     *   6.95,  6.94,  6.93,  6.92,  6.90,  6.91,  6.89,  6.91,  6.90,
     *   6.90,  6.91,  6.89,  6.93,  6.91,  6.93,  6.87,  6.90,  6.93,
     *   6.85,  6.86,  6.88,  6.88,  6.86,  6.93,  6.87,  6.86,  6.88,
     *   6.88,  6.87,  6.89,  6.88,  6.86,  6.88,  6.90,  6.91,  6.89,
     *   6.86,  6.90,  6.86,  6.90,  6.89,  6.90,  6.89,  6.89,  6.90,
     *   6.87,  6.91,  6.89,  6.88,  6.91,  6.91,  6.91,  6.88/
      data rTPR975/ 3221.94,115.97, 39.70, 24.65, 18.92, 15.63, 13.76,
     *  12.59, 11.56, 11.01, 10.52, 10.15,  9.88,  9.68,  9.37,  9.22,
     *   9.11,  9.02,  8.93,  8.71,  8.60,  8.61,  8.52,  8.46,  8.42,
     *   8.37,  8.39,  8.31,  8.30,  8.16,  8.22,  8.20,  8.13,  8.16,
     *   8.12,  8.15,  8.10,  8.04,  8.03,  8.07,  7.99,  8.05,  8.02,
     *   7.97,  7.99,  7.93,  7.94,  7.93,  7.92,  7.91,  7.91,  7.90,
     *   7.87,  7.92,  7.86,  7.88,  7.90,  7.91,  7.84,  7.86,  7.86,
     *   7.81,  7.79,  7.81,  7.85,  7.80,  7.89,  7.78,  7.79,  7.81,
     *   7.81,  7.81,  7.83,  7.82,  7.78,  7.81,  7.86,  7.83,  7.81,
     *   7.81,  7.80,  7.77,  7.82,  7.79,  7.81,  7.78,  7.82,  7.82,
     *   7.77,  7.81,  7.82,  7.77,  7.81,  7.83,  7.79,  7.77/
      data rTPR99/ 21411.19,304.57, 75.56, 42.31, 28.44, 22.64, 19.20,
     *  17.22, 15.32, 14.42, 13.57, 12.99, 12.52, 12.13, 11.67, 11.44,
     *  11.17, 11.00, 10.85, 10.59, 10.48, 10.43, 10.31, 10.20, 10.09,
     *   9.97, 10.03,  9.98,  9.84,  9.69,  9.76,  9.74,  9.65,  9.62,
     *   9.68,  9.60,  9.55,  9.46,  9.50,  9.54,  9.40,  9.53,  9.41,
     *   9.32,  9.34,  9.37,  9.32,  9.30,  9.28,  9.28,  9.23,  9.21,
     *   9.12,  9.24,  9.15,  9.15,  9.23,  9.28,  9.17,  9.11,  9.12,
     *   9.05,  9.08,  9.06,  9.08,  9.09,  9.18,  9.00,  9.04,  9.05,
     *   9.03,  9.04,  9.07,  9.01,  8.97,  9.05,  9.08,  9.01,  9.02,
     *   9.03,  9.05,  8.99,  9.01,  8.96,  8.99,  8.96,  8.96,  9.02,
     *   8.93,  9.03,  9.05,  8.97,  9.00,  8.95,  8.97,  8.96/

c     Set look-up tables of critical values of ts for one-sided
c     90%, 95%, and 97.5% significance levels.
c     The values   are taken from the 1986 and 1997 papers by Alexandersson 
c     and Alexandersson and Moberg, respectively.
      data  sig90/4.27, 5.05,6.10,6.65,7.00,7.25,7.40,7.55,7.70,7.80,
     &  7.85,8.05,8.20,8.35/
      data  sig95/4.54, 5.70,6.95,7.65,8.10,8.45,8.65,8.80,8.95,9.05,
     &  9.15,9.35,9.55,9.70/
      data sig975/4.71, 6.25,7.80,8.65,9.25,9.65,9.85,10.1,10.2,10.3,
     &  10.4,10.8,11.0,11.2/
      data nsig/5,  10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 
     &  200, 250/

c     SNIT one-sided critical value estimations
      if(iopt .eq. iTstat) then
c       no good if less than 10....
        if (num .lt. nsig(1)) then
          t90 = 99999.
          t95 = 99999.
          t975 = 99999.
c       if num is greater than 250 set sig-lev = sig at 250
        else if(num .ge. nsig(numsig)) then
          t90 = sig90(numsig)
          t95 = sig95(numsig)
          t975 = sig975(numsig)
        else
c         Select the critical value of t for the 3 significance
c         levels for the input number of months n
          do m = 1,numsig
            if(num .lt. nsig(m)) goto 10
          enddo
      
c         Estimate critical values of ts using linear interpolation
c         between table values -- 
   10     n21 = nsig(m) - nsig(m-1)
          nn1 = num - nsig(m-1)
          t90 = (( sig90(m)- sig90(m-1))/n21)*nn1 +  sig90(m-1)
          t95 = (( sig95(m)- sig95(m-1))/n21)*nn1 +  sig95(m-1)
          t975 = ((sig975(m)-sig975(m-1))/n21)*nn1 + sig975(m-1)
        endif  
      
      else if(iopt .eq. iFstat) then
      
c       Critical value estimation using the L&R Two Phase Regression tables
        if(num .lt. nTPR(1)) then
          t90 = 99999.
          t95 = 99999.
          t975 = 99999.
        else if(num .ge. nTPR(numTPR)) then
          t90 = 99999.
          t95 = rTPR95(numTPR)
          t975 = rTPR975(numTPR)
        else
c         Select the critical value of two-phase
          do m = 1,numtpr
            if(num .lt. nTPR(m)) goto 20
          enddo
   20     t90 = 99999.
          t95 = rTPR95(m)
          t975 = rTPR975(m)
        endif
        
      else
        write(6,*)' Unknown critical value option: ', iopt
        stop 999  
        
      endif    

      return
      end

c ===================================================================

      subroutine ttest(kmo, z, end1, n2, tt, minseg)
      include 'inhomog.parm.mthly.incl'
      
c     routine to calculate student's t test
c     from Statistical Analysis for Climate Research pp. 111-118
c          section 6.6.5 Unequal Variances

c     z - input array of realizations (assume no missing)
c     n1 - begin of the first realization
c     n2 - end of the second realization
c     end1 - breakpoint between populations
c     tt - t-test statistic

c     xsq1, xbar1, xn1 - sx^2/nx, avg x, num x 
c     xsq2, xbar2, xn2 - sy^2/ny, avg y, num y 
c     

      real z(kmo),xsq1,xsq2,xsq3,xn1,xn2,xbar1,xbar2,tt
      integer k,end1,n1,n2

      n1 = 1

c     The first population
      xsq1 = 0.0
      xn1 = 0.0
      xbar1 = 0.0
      xbar2 = 0.0
      do k = n1, end1
        if(z(k) .gt. amiss + 1.) then
          xbar1 = xbar1 + z(k)
          xn1 = xn1 + 1.0
        endif  
      enddo
      
      if (xn1 .ge. minseg) then
        xbar1 = xbar1 / xn1
      else 
        tt = 0.0
        return
      end if

      do k = n1, end1
        if(z(k) .gt. amiss + 1.) then
          xsq1 = xsq1 + (z(k) - xbar1)**2.
        endif  
      enddo  
      xsq1 = xsq1 / (xn1-1.0) / xn1
        
c     The second population 
      xn2 = 0.0
      xsq2 = 0.0
      do k = end1+1, n2
        if(z(k) .gt. amiss + 1.) then
          xbar2 = xbar2 + z(k)
          xn2 = xn2 + 1.0
        endif  
      enddo
      
      if (xn2 .ge. minseg) then
        xbar2 = xbar2 / xn2
      else
        tt=0.0
        return
      end if
      
      do k = end1+1, n2
        if(z(k) .gt. amiss + 1.) then
          xsq2 = xsq2 + (z(k) - xbar2)**2.
        endif  
      enddo  
      xsq2 = xsq2 / (xn2-1.0) / xn2
      
c     now, put them all together
      xsq3 = xsq1 + xsq2
      if(xsq3 .lt. eps) xsq3 = eps
      tt = abs(xbar1-xbar2) / sqrt(xsq3)

      return
      end

c ===================================================================

      subroutine tpr1table (num,F95)
c**********************************************************************
c
c     Subroutine lookup computes the critical values for the TPR1
c     changepoint model (step change with a constant slope).  
c     The critical values are dependent on the sample size -- num.
c
c                               -on input-
c
c       num     The number of values in a series
c
c                               -on output-
c
c       F95     one-sided 95% confidence level
c
c     The values are taken from the 2003 paper by Wang
c     J,Climate, 16, 3383-3385 (Table 1) except for n = 5
c     which was obtained via simulations of the null hypothesis
c     M.J. Menne -- 3 Sept 2005.  
c
c***********************************************************************
      dimension sig95(17)
      integer nsig(17)

c     Set look-up tables of critical values of Fc for one-sided
c     95% significance levels.

      data  sig95/63.67,15.56,11.34,11.15,11.07,11.06,11.07,11.06,
     &11.07,11.08,11.085,11.127,11.208,11.310,11.537,11.749,12.064/
      data nsig/5,10,20,30,40,50,60,70,80,90,100,150,200,250,500,1000,
     &2500/

c     no good if less than 5....
      if (num.lt.5) then
        F95 = 9999
        return
      endif
      
c     if num is greater than 2500 set sig-lev = sig at 2500
      if(num .ge. 2500) then
        F95 = sig95(17)
        return
      endif

c     Select the critical value of t for the 3 significance
c     levels for the input number of months n
      do m = 1,17
        if(num .lt. nsig(m)) goto 10
      enddo
      
c     Estimate critical values of ts using linear interpolation
c     between table values -- 
   10 n21 = nsig(m) - nsig(m-1)
      nn1 = num - nsig(m-1)
      F95 = (( sig95(m)- sig95(m-1))/n21)*nn1 +  sig95(m-1)
      return
      end
c

c ===================================================================

      subroutine twophreg(n,z,Fc,mcnt)

c     v2 07sep05 includes the impact of missing data which is especially
c         significant with sloped segment

c**********************************************************************
c
c     Subroutine snits computes the standard normal inhomogeneity test
c     statistic for a step change or shift in the standardized
c     difference data between a candidate and the average of its
c     neighbors.  The algorithm to calculate ts is Eq.(4) in Alexander-
c     sson and Moberg, Int'l Jour. Climat., 17, 25-34 (1997).  
c
c     23 September, 1999; 18 October, 1999; 9 December, 1999
c     6 July, 2001; November 2001; April 2002--to include date intervals
c
c
c                               -on input-
c
c       n       The number months or years for z & ts arrays
c
c       mcnt    the number of non-missing months or years
c               in the time series (n and mcnt can be combined
c               if no missing months are passed to the routine)
c
c       z       The n by 1 z-series of standardized monthly
c               temperature differences produced by sbqzoft
c
c                               -on output-
c
c       Fc      The mcnt series of the Fmax test statistic
c       
c***********************************************************************
      INCLUDE 'inhomog.parm.mthly.incl' 

      dimension z(n),Fc(n)

c     Calculate the Fmax test statistic assuming no change in slope
c     between the two time series

c      calculate SSE for reduced model

      do i = 1, n
        Fc(i) = amiss
      enddo  

      SSEred = 0.0
      zmnR = 0.0
      rMnR = 0.0
      alphaR = 0.0
      tmnR = 0.0
      nR = 0
      rNumR = 0.0
      rDenR = 0.0
      
c  Get the least squares estimate of the mean and the slope (muhat and
c  alphahat)

      do i=1,mcnt
        if(z(i) .gt. amiss + 1.0) then
          zmnR = zmnR + z(i)
          tmnR = tmnR + i
          nR = nR + 1
        endif  
      end do
      if(nR .ne. 0) then
        zmnR = zmnR / nR
        tmnR = tmnR / nR
      else
        return
      endif
      
      do i=1,mcnt
        if(z(i) .gt. amiss + 1.0) then
          rNumR = rNumR + (i - tmnR)*(z(i) - zmnR)
          rDenR = rDenR + (i - tmnR)*(i - tmnR)
        endif  
      enddo
      
      alphaR = rNumR / rDenR
      rMuR = zmnR - alphaR * tmnR
      
      do i=1,mcnt
        if(z(i) .gt. amiss + 1.0)
     *    SSEred = SSEred + (z(i) - rMuR - alphaR*i)**2
      enddo
      
c
c  calculate the SSEs for the "full model"
c
      
      do i = 2, mcnt - 1
c       the output stat series(i) == missing if z(i) == missing
        if(z(i) .lt. amiss + 1.0) goto 10
 
        SSEful = 0.0
        zmn1 = 0.0
        tmn1 = 0.0
        n1 = 0
        do j = 1,i
          if(z(j) .gt. amiss + 1.0) then
            zmn1 = zmn1 + z(j)
            tmn1 = tmn1 + j
            n1 = n1 + 1
          endif  
        enddo
        
        if(n1 .ne. 0) then     
          zmn1 = zmn1/n1
          tmn1 = tmn1/n1
        else
          goto 10
        endif 
        
        rNum1 = 0.0
        rDen1 = 0.0
        alpha1 = 0.0
        rMu1 = 0.0
        SSE1 = 0.0
        do j=1,i
          if(z(j) .gt. amiss + 1.0) then
            rNum1 = rNum1 + (j - tmn1)*(z(j) - zmn1)      
            rDen1 = rDen1 + (j - tmn1)*(j - tmn1)
          endif  
        enddo
          
      
        alpha1 = rNum1 / rDen1
        rMu1 = zmn1 - alpha1 * tmn1
        
c       Mean for first segment ready calculate SSE
        
        do j = 1,i
          if(z(j) .gt. amiss + 1.0) 
     *      SSE1 = SSE1 + (z(j) - rMu1 - alpha1*j)**2
        enddo
        
        zmn2 = 0.0
        tmn2 = 0.0
        n2 = 0
        alpha2 = 0.0
        rMu2 = 0.0
        rNum2 = 0.0
        rDen2 = 0.0
        SSE2 = 0.0
        if(i+1 .lt. mcnt) then
          do j = i+1, mcnt
            if(z(j) .gt. amiss + 1.0) then
              zmn2 = zmn2 + z(j)
              tmn2 = tmn2 + j
              n2 = n2 + 1
            endif  
          enddo
          if(n2 .ne. 0) then
            zmn2 = zmn2/n2
            tmn2 = tmn2/n2
          else
            goto 10
          endif
        
          do j=i+1,mcnt
            if(z(j) .gt. amiss + 1.0) then
              rNum2 = rNum2 + (j - tmn2)*(z(j) - zmn2)
              rDen2 = rDen2 + (j - tmn2)*(j - tmn2)
            endif  
            enddo
        
          alpha2 = rNum2 / rDen2
          rMu2 = zmn2 - alpha2 * tmn2
        
          if (i+1 .lt. mcnt) then
            do j = i+1,mcnt
              if(z(j) .gt. amiss + 1.0)
     *          SSE2 = SSE2 + (z(j) - rMu2 - alpha2*j)**2
            enddo
          endif
        endif

        SSEful = SSE1 + SSE2
        Fc(i) = ((SSEred - SSEful)/2.)/(SSEful/(mcnt-4))
   10 enddo

      return
      end
      
c ===================================================================
      subroutine twophreg1(n,z,F,mcnt)
!c**********************************************************************
!c
!c     Subroutine twophreg computes the two-phase test
!c     statistic for a step change and/or trend shift in a series
!c     using the method described by Lund and Reeves (2002)
!c
!c                               -on input-
!c
!c       n       The number months or years in the z-series
!c
!c       mcnt    the number of non-missing months or years
!c               in the time series (n and mcnt can be combined
!c               if no missing months are passed to the routine)
!c
!c       z       The n by 1 z-series of standardized monthly
!c               temperature differences produced by sbqzoft
!c
!c         rMiss   Missing value indicator, assumed to be a large negative
!c               value
!c
!c                               -on output-
!c
!c       Fc      The mcnt series of the snit test statistic
!c               ts for a shift or step-change in monthly mean
!c               temperature at the candidate station
!c       
!c***********************************************************************
      INCLUDE 'inhomog.parm.mthly.incl' 

      real z(n)
      real F(n)
      real resid(n)
      
      rmiss = amiss

!c      Initialize Fc statistic
      do iFc=1,n
        F(iFc) = amiss
      enddo
      
!c
!c     Calculate the Fmax test statistic assuming no change in slope
!c     between the two time series
!c
!c
!c      calculate SSE for reduced model

      SSEred = 0.0
      zmnR = 0.0
      rMnR = 0.0
      alphaR = 0.0
      tmnR = 0.0
      nR = 0
      rNumR = 0.0
      rDenR = 0.0
      
!c    Get the least squares estimate of the mean and the slope (muhat and
!c    alphahat) for "reduced" model

      do i=1,mcnt
        if(z(i) .gt. rMiss + 1) then
          zmnR = zmnR + z(i)
          tmnR = tmnR + i
          nR = nR + 1
        endif
      end do
      
      if(nR .ne. 0) then
        zmnR = zmnR / nR
        tmnR = tmnR / nR
      else
        return
      endif
      
      do i=1,mcnt
        if(z(i) .gt. rMiss + 1) then
          rNumR = rNumR + (i - tmnR)*(z(i) - zmnR)
          rDenR = rDenR + (i - tmnR)*(i - tmnR)
        endif
      enddo
      
      alphaR = rNumR / rDenR
      rMuR = zmnR - alphaR * tmnR
      
      do i=1,mcnt
        if(z(i) .gt. rMiss + 1)
     *    SSEred = SSEred + (z(i) - rMuR - alphaR*i)**2
      enddo
      
!c
!c    calculate the SSE for the "full model" with constant slope
!c
      
      do i = 2,mcnt - 2
        if(z(i) .lt. rMiss + 1) goto 10
       
        zmn1 = 0.0
        zmn2 = 0.0
        zmn = 0.0
        tmn1 = 0.0
        tmn2 = 0.0
        n1 = 0
        n2 = 0
        nF = 0
        alpha = 0.0
        rMu1 = 0.0
        rMu2 = 0.0
        rNum = 0.0
        rDen = 0.0
        SSE1 = 0.0
        SSE2 = 0.0
        SSEful = 0.0
        
        do j = 1,i
          if(z(j) .gt. rMiss + 1) then
            zmn1 = zmn1 + z(j)
            tmn1 = tmn1 + j
            n1 = n1 + 1
          endif
        enddo
        
        if(n1 .ne. 0) then     
          zmn1 = zmn1/n1
          tmn1 = tmn1/n1
        else
          goto 10
        endif 
        
        do j=1,i
          if(z(j) .gt. rMiss + 1) then
            resid(j) = z(j) - zmn1
          else
            resid(j) = rMiss
          endif
        enddo
        
!c      Resids for first segment ready
        
!c      Calculate mean for second segment

        do j = i+1, mcnt
          if(z(j) .gt. rMiss + 1) then
            zmn2 = zmn2 + z(j)
            tmn2 = tmn2 + j
            n2 = n2 + 1
          endif
        enddo
        
        if(n2 .ne. 0) then
          zmn2 = zmn2/n2
          tmn2 = tmn2/n2
        else
          goto 10
        endif
        
!c      Resids for second segment ready 
        
        do j=i+1,mcnt
          if(z(j) .gt. rMiss + 1) then
            resid(j) = z(j) - zmn2
          else
            resid(j) = rMiss
          endif
        enddo
        
!       now calculate common slope for two segments
        do j=1,mcnt
          if(resid(j) .gt. rMiss + 1) then
            zmn = zmn + resid(j)
            nF = nF + 1
          endif
        enddo
        
        if (nF .ne. 0) then
          zmn = zmn / nF
        else
          stop 'zmn twophreg1'
        endif

        do j=1,i
          if(resid(j) .gt. rMiss + 1) then
            rNum = rNum + (j-tmn1)*(resid(j))
            rDen = rDen + (j-tmn1)*(j-tmn1)
          endif
        enddo
        
        do j=i+1,mcnt
          if(resid(j) .gt. rMiss + 1) then
            rNum = rNum + (j-tmn2)*(resid(j))
            rDen = rDen + (j-tmn2)*(j-tmn2)
          endif
        enddo
      
        alpha = rNum / rDen
      
!       now calculate SSE1 and SSE2
        rMu1 = zmn1 - alpha * tmn1
        rMu2 = zmn2 - alpha * tmn2

        do j = 1,i        
          if(z(j) .gt. rMiss + 1) then
            SSE1 = SSE1 + (z(j) - rMu1 - alpha*j)**2
          endif
        enddo
          
        if(i+1 .le. mcnt-2) then
        
          do j = i+1,mcnt
            if(z(j) .gt. rMiss + 1)
     *        SSE2 = SSE2 + (z(j) - rMu2 - alpha*j)**2
          enddo
        
        endif
        
        SSEful = SSE1 + SSE2
        F(i) = ((SSEred - SSEful)/1.)/(SSEful/(float (mcnt-3)))
        
   10 enddo
      return
      end 
      
c ===================================================================

      subroutine kendalltheil(nObs, rYdata, rYint, rSlope,
     *  rSumSqrE, iN, islope, ierr )
!-----------------------------------------------------------------------
! Subroutine to calculate the kendall_theil slope, intercept and
!     probability for the input series
!      Based upon ranking all the paired combinations of the observed
!        points.
!   NOTE: for this version, Y is assumed to be a monotonically increasing
!        array (months) WITH MISSING DATA
!
!                          -on input-
!  nObs - sample size (including missing)
!  rYdata - dependent data (predictand)
!  islope - the option to force the slope to 0 for certain models
!            0 = force slope to zero; else calculate slope
!            additionally - if islope is less than nObs then
!                         separate the segments at that point
!                         (for TPR1 constant slope emulation)
!
!                          -on output-
!  rYint - y - intercept of regression line for each segment (see islope)
!  rResid - residual of the fit
!  rFit - fitted regression line
!  rSlope - slope of linear regression line
!  rT - slope error
!  rSumSqrE - sum square of residuals
!  iN - total number of non-missing values used
!  ierr - algorithm unable to compute regression 
!
!-----------------------------------------------------------------------      
      INCLUDE 'inhomog.parm.mthly.incl' 
      INCLUDE 'inhomog.MDparm.mthly.incl'

      parameter(minsbin=-1000,maxsbin=1000, sinc=.0002)
       integer nObs,iN 
      real rSlope,rYint(2),rMiss
      real rYdata(kmo),rResid(kmo),rFit(kmo),rTemp(60*60)
      real rX(kmo),rY(kmo)
      integer islpbin(minsbin:maxsbin)
      
      ierr = 0
      rMiss = amiss
      iN = 0
      lastgood = 0
      iend1 = 0
      ibeg2 = 0
      rT = -9999.0
      rSumSqrE = 0.0
      ktdebug = 0
      hsinc = sinc * 0.5
      
c     remove missing data from incoming arrays
      do i = 1, nObs
        if(rYdata(i) .ne. rMiss) then
          iN = iN + 1
          rX(iN) = real(iN)
          rY(iN) = rYdata(iN)
          if(ibeg2.eq.0 .and. (islope.le.i .and. islope.ne.0))then
            iend1 = lastgood  
            ibeg2 = iN
          endif  
        endif
        rResid(i) = amiss
        rFit(i) = amiss
        lastgood = iN
      enddo
      if(iend1 .eq. 0 .or. islope .eq. 0) iend1 = iN
      if(iN .lt. 5) then
        ierr = 1
        return
      endif  

      if(ktdebug .gt. 0) 
     *  write(6,'("iN, iend1, ibeg2:",3i5)') iN, iend1, ibeg2

      if(iN .lt. 60) then
c       original sort ascending method
        nslp = 0
c       Generate paired slopes for the first segment
        do i = 1, iend1 - 1
          do j = i+1, iend1
c           if x(j) != x(i) then calc & save paired slope
            if((rX(j) - rX(i)) .ne. 0.0) then
              nslp = nslp + 1
              rTemp(nslp) = (rY(j)-rY(i))/(rX(j)-rX(i))
            endif
          enddo
        enddo
c       Generate paired slopes for the second segment (if any)
        do i = ibeg2, iN - 1
          do j = ibeg2+1, iN
c           if x(j) != x(i) then calc & save paired slope
            if((rX(j) - rX(i)) .ne. 0.0) then
              nslp = nslp + 1
              rTemp(nslp) = (rY(j)-rY(i))/(rX(j)-rX(i))
            endif
          enddo
        enddo
      
c       if islope option is set to 0, force slope to zero
        if(islope .eq. 0) then
          rSlope = 0.0
        else  
c         kt slope is the median of all values
          call SortAscending(nslp, rTemp)
          imed = nslp/2
          if(mod(nslp,2) .eq. 1) imed = imed + 1        
          rSlope = rTemp(imed)
          if(ktdebug .gt. 0) write(6,'("slope, ic, imet:",f7.2,2I5)')
     *      rSlope,nslp,imed
        endif
      
      else
      
c      if(iN .gt. 60) then
c       new "bin"ed slope method when number of values are large
        do ind = minsbin,maxsbin
          islpbin(ind) = 0
        enddo
        
        nslp = 0
c       Generate paired slopes for the first segment
        do i = 1, iend1 - 1
          do j = i+1, iend1
c           if x(j) != x(i) then calc & save paired slope
            if((rX(j) - rX(i)) .ne. 0.0) then
              nslp = nslp + 1
              rT = (rY(j)-rY(i))/(rX(j)-rX(i))
              if(rt .gt. -1*hsinc) then
                ind = (rt + hsinc)/sinc
                if(ind .gt. maxsbin) ind = maxsbin
              else
                ind = (rt - hsinc)/sinc
                if(ind .lt. minsbin) ind = minsbin
              endif
              islpbin(ind) = islpbin(ind) + 1    
            endif
          enddo
        enddo
c       Generate paired slopes for the second segment (if any)
        do i = ibeg2, iN - 1
          do j = ibeg2+1, iN
c           if x(j) != x(i) then calc & save paired slope
            if((rX(j) - rX(i)) .ne. 0.0) then
              nslp = nslp + 1
              rT = (rY(j)-rY(i))/(rX(j)-rX(i))
              if(rt .gt. -1*hsinc) then
                ind = (rt + hsinc)/sinc
                if(ind .gt. maxsbin) ind = maxsbin
              else
                ind = (rt - hsinc)/sinc
                if(ind .lt. minsbin) ind = minsbin
              endif
              islpbin(ind) = islpbin(ind) + 1    
            endif
          enddo
        enddo
      
        if(mod(nslp,2) .eq. 1) then
          imed = (nslp + 1)/2
        else
          imed = nslp / 2
        endif
c        print *,' nslp, imed', nslp, imed
        
        isum = 0
        do ind = minsbin, maxsbin
          isum = isum + islpbin(ind)
c          if(ktdebug .gt. 0 .and. islpbin(ind) .gt. 0) 
c     *      print *,' ind, bin, hist, cum ', ind, 
c     *        real(ind)*sinc + hsinc,islpbin(ind), isum
          if(imed .le. isum) goto 100
        enddo
        
c       if islope option is set to 0, force slope to zero
  100   if(islope .eq. 0) then
          rSlope = 0.0
        else  
c         kt slope is the median of all values
          sdif = real(isum-imed)/real(islpbin(ind))*sinc
          sind = real(ind)*sinc + hsinc
          rslope=sind - sdif
        endif
        
      endif
      
c     Always make the first segment intercept from 1 to iend1
c     intercept is y-median - slope * x-median
      nval1 = iend1
      imed = nval1/2
      if(mod(nval1,2) .eq. 1) imed = imed + 1        
c     X is already sorted... (see assumption @ start of subroutine)
      rXmed = rX(imed)
      call SortAscending(nval1, rY)
      rYmed = rY(iMed)
      rYint(1) = rYmed - rslope * rXmed
      if(ktdebug .gt. 0) 
     *  write(6,'("Seg1 - Xmed, Ymed, slope, Yint:",3f7.2,f7.3)')
     *   rXmed, rYmed, rslope, rYint(1)
      
c     if there is a segment after the break - do it
      if(ibeg2 .lt. iN) then
        nval2 = iN - ibeg2 + 1
        do i = ibeg2, iN
          rX(i - ibeg2 + 1) = rX(i)
          rY(i - ibeg2 + 1) = rY(i)
        enddo
        imed = nval2/2
        if(mod(nval2,2) .eq. 1) imed = imed + 1        
        rXmed = rX(imed)
        call SortAscending(nval2, rY)
        rYmed = rY(iMed)
        rYint(2) = rYmed - rslope * rXmed
        if(ktdebug .gt. 0) 
     *    write(6,'("Seg2 - Xmed, Ymed, slope, Yint:",3f7.2,f7.3)')
     *     rXmed, rYmed, rslope, rYint(2)
      endif  

c     setup arrays to determine root mean square
      iN = 0
      rSumSqrX = 0.0
      do i = 1, nObs
        if(rYdata(i) .ne. rMiss) then
          iN = iN + 1
          if(iN .le. iend1) then
            Yint = rYint(1)
          else
            Yint = rYint(2)
          endif    
          rResid(i) = (Yint + rSlope*i) - rYdata(i)
          rFit(i) = Yint + rSlope*i
          rSumSqrE = rSumSqrE + rResid(i)**2
          rSumSqrX = rSumSqrX + (real(i) - rXmed)**2
        endif
      enddo 
      
      if (iN .ge. 5) then
        rSeSqr = rSumSqrE / (iN -2)
        rSb = sqrt(rSeSqr / rSumSqrX)
        rT = rSlope / rSb
        if(ktdebug .gt. 0) write(6,'("sse, rT:",2f7.2)')rSumSqrE,rT
      endif

      return
      end 

c ===================================================================
c     this set of routines are basic average, median, linear slope
      
c ---------------------------------------------------------------------

      subroutine lsqline(n, rMiss, xData, yData, xmean, ymean, slope, 
     *  yInt, sseflat, sseslope, count)
!-----------------------------------------------------------------------
! Subroutine to calculate least squares simple linear regression line
! given an array of dependent (yData) and independent (xData) data 
! values.
!                          -on input-
!  n - sample size
!  rMiss - Missing value code (assumed to be a large negative number) 
!  xData - independent data (predictor)
!  yData - dependent data (predictand)
!
!                          -on output-
!  xmean - arithmetic average of x
!  ymean - arithmetic average of y
!  slope - slope of linear regression line
!  yInt - y - intercept of regression line
!  resid - array of least squares residuals
!
!-----------------------------------------------------------------------      
       integer n, count
      real sumx,sumy,xmean,ymean,slope,yint,rMiss
      real xdata(n), ydata(n)

      count = 0
      sumx = 0
      sumy = 0
      xmean = rmiss
      ymean = rmiss
      slope = rmiss
      yint = rmiss
      sseflat = 0.0
      sseslope = 0.0
      
      do i=1,n
        if ((xdata(i).gt.rMiss + 1).and.(ydata(i).gt.rMiss + 1)) then
          count = count + 1.0
          sumx = sumx + xdata(i)
          sumy = sumy + ydata(i)
        endif
      enddo        

      if(count .ge. 2) then
        xmean = sumx / count
        ymean = sumy / count
      else
        return
      endif  
          
!      calculate slope and y-intercept
      rnum1 = 0.0
      rden1 = 0.0
      do i = 1, n
        if ((xdata(i).gt.rMiss + 1).and.(ydata(i).gt.rMiss + 1)) then
          rNum1 = rNum1 + (xdata(i)-xmean)*(ydata(i)-ymean)      
          rDen1 = rDen1 + (xdata(i)-xmean)*(xdata(i)-xmean)
        endif
      enddo
        
      slope = rNum1 / rDen1
      yint = ymean - slope * xmean

!     calculate residuals
      do i=1,n
        if ((xdata(i).gt.rMiss + 1).and.(ydata(i).gt.rMiss + 1)) then
          sseflat = sseflat + (ymean - ydata(i)) ** 2
          resid = (yint + slope*xdata(i)) - ydata(i)
          sseslope = sseslope + resid ** 2
        endif  
      enddo 
      
      return
      end 

c ===================================================================


      subroutine kthline(inval, rMiss, xData, yData, xmed, ymed, slope,
     *  yInt, sseflat, sseslope, nval)
!-----------------------------------------------------------------------
! Subroutine to calculate least squares simple linear regression line
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
!  resid - array of least squares residuals
!
!-----------------------------------------------------------------------      
      INCLUDE 'inhomog.parm.mthly.incl'
      parameter(minsbin=-1000,maxsbin=1000, sinc=.0002)
       integer inval, nval, nslp
      real xmed,ymed,slope,yint,rMiss
      real xdata(nmo), ydata(nmo)
      real rX(nmo), rY(nmo), rTemp(60*60)
      integer islpbin(minsbin:maxsbin)

      nval = 0
      xmed = rmiss
      ymed = rmiss
      slope = rmiss
      yint = rmiss
      sseflat = 0.0
      sseslope = 0.0
      hsinc = sinc * 0.5
      
      do i=1,inval
        if ((xdata(i).gt.rMiss + 1).and.(ydata(i).gt.rMiss + 1)) then
          nval = nval + 1
          rx(nval) = xdata(i)
          ry(nval) = ydata(i)
        endif
      enddo        
      if(nval .lt. 2) return

c      if(nval .lt. 5) then
c        print *,' Low count'
c      endif  

!     calculate slope
      if(nval .lt. 60) then
c       use the sort-ascending for low number of values      
        nslp = 0
        do i = 1, nval-1
          do j = i+1, nval
            if(rX(j) .ne. rX(i)) then
              nslp = nslp + 1
              rTemp(nslp) = (rY(j)-rY(i))/(rX(j)-rX(i))
            endif
          enddo    
        enddo

        call SortAscending(nslp, rTemp)
        if(mod(nslp,2) .eq. 1) then
          imed = (nslp + 1)/2
          slope = rTemp(imed)
        else
          imed = nslp / 2
          slope = (rTemp(imed) + rTemp(imed+1))/2
        endif
        
      else
        
c      if(nval .gt. 60) then
c       use bined distribution when number of values are large
        do ind = minsbin,maxsbin
          islpbin(ind) = 0
        enddo
        
        nslp = 0
        do i = 1, nval-1
          do j = i+1, nval
            if(rX(j) .ne. rX(i)) then
              nslp = nslp + 1
              rT = (rY(j)-rY(i))/(rX(j)-rX(i))
              if(rt .gt. -1*hsinc) then
                ind = (rt + hsinc)/sinc
                if(ind .gt. maxsbin) ind = maxsbin
              else
                ind = (rt - hsinc)/sinc
                if(ind .lt. minsbin) ind = minsbin
              endif
              islpbin(ind) = islpbin(ind) + 1    
            endif
          enddo
        enddo

        if(mod(nslp,2) .eq. 1) then
          imed = (nslp + 1)/2
        else
          imed = nslp / 2
        endif
        
        isum = 0
        do ind = minsbin, maxsbin
          isum = isum + islpbin(ind)
          if(imed .le. isum) goto 100
        enddo
        
  100   sdif = real(isum-imed)/real(islpbin(ind))*sinc
        sind = real(ind)*sinc + hsinc
        slope=sind - sdif
        
      endif
      
!     calculate y-intercept
      call SortAscending(nval, rx)
      call SortAscending(nval, ry)
      if(mod(nval,2) .eq. 1) then
        imed = (nval + 1) / 2
        xmed = rx(imed)
        ymed = ry(imed)
      else
        imed = nval / 2
        xmed = (rx(imed) + rx(imed+1)) / 2  
        ymed = (ry(imed) + ry(imed+1)) / 2  
      endif  
                  
      yint = ymed - slope * xmed

!     calculate residuals
      do i=1,inval
        if ((xdata(i).gt.rMiss + 1).and.(ydata(i).gt.rMiss + 1)) then
          sseflat = sseflat + (ymed - ydata(i)) ** 2
          resid = (yint + slope*xdata(i)) - ydata(i)
          sseslope = sseslope + resid ** 2
        endif  
      enddo 
      
      return
      end 

C     ************************************************************

      subroutine sort1(n,arr)

c     Shell's method.

c     n is the number of elements to sort in array ARR, 
c     ndx is an associated array that is sorted along with ARR
      real arr(n)

      aln2i=1./0.69314718
      tiny=1.E-5
      lognb2=int(alog(float(n))*aln2i+tiny)
      m=n
      do nn=1,lognb2
        m=m/2
        k=n-m
        do j=1,k
          i=j
 553      continue
          l=i+m
c         descending sort....
c          if(arr(l).gt.arr(i))then
c         ascending sort....
          if(arr(l).lt.arr(i))then
            t=arr(i)
            arr(i)=arr(l)
            arr(l)=t
            i=i-m
            if(i.ge.1)go to 553
          endif
        enddo
      enddo

 900  continue
      return

      end

!******************************************************************************

      subroutine SortAscending(N,X)

  ! PURPOSE: Sort an array X into ascending order.

  ! AUTHOR(S): Public Domain - www.netlib.org/napack/sort.f 

      INTEGER I,J,K,L,M,N
      REAL X(N),Y(N),S,T
      I = 1
  10  K = I
  20  J = I
      I = I + 1
      IF ( J .EQ. N ) GOTO 30
      IF ( X(I) .GE. X(J) ) GOTO 20
      Y(K) = REAL(I)
      GOTO 10
  30  IF ( K .EQ. 1 ) RETURN
      Y(K) = REAL(N + 1)
  40  M = 1
      L = 1
  50  I = L
      IF ( I .GT. N ) GOTO 120
      S = X(I)
      J = INT(Y(I))
      K = J
      IF ( J .GT. N ) GOTO 100
      T = X(J)
      L = INT(Y(J))
      X(I) = L
  60  IF ( S .GT. T ) GOTO 70
      Y(M) = S
      M = M + 1
      I = I + 1
      IF ( I .EQ. K ) GOTO 80
      S = X(I)
      GOTO 60
  70  Y(M)= T
      M = M + 1
      J = J + 1
      IF ( J .EQ. L ) GOTO 110
      T = X(J)
      GOTO 60
  80  Y(M) = T
      K = M + L - J
      I = J - M
  90  M = M + 1
      IF ( M .EQ. K ) GOTO 50
      Y(M) = X(M+I)
      GOTO 90
  100 X(I) = J
      L = J
  110 Y(M) = S
      K = M + K - I
      I = I - M
      GOTO 90
  120 I = 1
  130 K = I
      J = X(I)
  140 X(I) = Y(I)
      I = I + 1
      IF ( I .LT. J ) GOTO 140
      Y(K) = REAL(I)
      IF ( I .LE. N ) GOTO 130
      IF ( K .EQ. 1 ) RETURN
      GOTO 40

      end subroutine SortAscending 

