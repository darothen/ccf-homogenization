 
      subroutine splitmerge(idunit, rrtemp, subnet, iopt, istnpt,
     *  isigy1, isigy2, mtype, asigx, azscr, nsig, iwrite, itarg, 
     *  it2pair, ipair, ip2targ, rTraw, imob, imoe)
     
c Vers   Date    Description
c
c  21f 19dec06  allow any tpr model to pass tpr0 model consistency test in testseg 
c
C     General note: During the splitting calls to TESTSEG there is a set of
c       work arrays used (with "ex" in their names) to accumulate the changes
c       to the changepoint arrays. All internal results in TESTSEG must be
c       accounted for in these arrays. 
c       During the merging calls to TESTSEG, changes are made only to the
c       incoming arrays (with "in" in their names)

c     subroutine parameters
c       idunit = output unit number for writing station data
c       rrtemp = cand & net temperature data, 12 monthly + annual
c       subnet = cand & net coop numbers
c       iopt definitions come from inhomog.comm.mthly.incl...
c       isigy1 = year-month changepoints (end of before (earlier) window)
c       isigy2 = year-month changepoints (begin of after (more current) window)
c       asigx = associated offsets for changepoint years
c       azscr = offset z-score wrt MINBIC segmented sum-square error
c       mtype = associated data models for changepoint years
c       nsig = number of changepoints
c       iwrite = toggle for writing out station data
c           :=0 for no output, else :=index of candidate for output
c       itarg = index of the first (target) station
c       ipair = index of the second station
c

      INCLUDE 'inhomog.parm.mthly.incl' 
      INCLUDE 'inhomog.comm.mthly.incl'
      INCLUDE 'inhomog.restart.mthly.incl'
      INCLUDE 'inhomog.MDparm.mthly.incl'
c
      character*13 stnstr
      character*2 otag
      character*6 subnet(nstns)
      integer ounit
      integer job

c     paired difference temperature series
      dimension rTraw(nmo)
      dimension qx(nmo),qy(nmo)
      
c     inhomogeneous breakpoint array (including begin/end of series)
      integer inhyr(ninh), inhmod(ninh)
      integer lsplitin2(ninh), lsplitout(ninh), lmergein2(ninh),
     *  lmergeout(ninh)
      real offinh(ninh),segslp(2,ninh),sseseg(ninh)
c     keep track of significant year/months (isigy), models (mtype), 
c       and offset (asigx)
      integer isigy1(ninh), isigy2(ninh), mtype(ninh), istnpt(ninh)
      real asigx(ninh), azscr(ninh)
c     inhomog hits delineated by none = 0, stn1 = 1, stn2 = 2, both = 3
      integer inhstns(ninh)
c     current lowest imo for each type of inhomog hit (1,2,3)
      integer inhmo(3)
      integer beg1,end1,beg2,end2
      real rrtemp(2,nmo)
      integer mknt2(2)
      real rslpq(2), rmuq(2)
      otag = c2tech(iopt)
      
      mintest = 5
c      mintest = minlen

c     >>>>>>>>>>>>> Initialize incoming inhomog arrays <<<<<<<<<<<<

c ------  Initialize data and reference series  ------
c     set up the temp arrays and the yr position array
c     imob & imoe are the begin and end indices of paired station's data
      imob = 0
      imoe = nmo
      do imo = 1, nmo
        if(rrtemp(1, imo) .gt. amiss+1. .and. 
     *    rrtemp(2, imo) .gt. amiss+1. ) then
          if(imob .eq. 0) imob = imo
          imoe = imo
          rTraw(imo) = rrtemp(1,imo) - rrtemp(2,imo)
        else
c          print *,' One of pair miss: ', imo
          rTraw(imo) = amiss  
        endif  
      enddo
      if(imob .eq. 0) then
        write(6,*)subnet(1),' Station has no data - skipping'
        goto 999
      endif
      call imo2iym(iyr,imth,imob)
      iBegYr = iyr
      iBegMth = imth
      call imo2iym(iyr,imth,imoe)
      iEndYr = iyr  
      iEndMth = imth
      
c     Initialize arguments to test segments:
c      iFirst requires the series to be split the first pass
c      numchgpt initializes number of change points found
C     Mindy determines if ifirst gets set to 1
      iFirst = 0
      ichange = 1
      ipass = 1

c     ARRAYS FOR CHGPTS FROM BIC and TESTSEG
c     inhyr is the month index of the chgpts
c     inhyr(1) is the begin of series, inhyr(mindy) is the end of series
c     inhmod is the model
c     inhstns is the station of the pair with the SHF "hit" (1,2,or 3=both)
c     offinh is the offset
c     BEFORE/AFTER ADJOINING SEGMENTS:
c     segslp are the slopes
c     sseseg are the sum-square error

      do inh = 1, ninh
        inhyr(inh) = 0
        inhstns(inh) = 0
      enddo  
          
c     number obs in current segment wrt inhomog date array (from station history)
      nobs = 0
c     first inhomog date array postion is beginning of data
      mindy = 1
      inhyr(mindy) = imob
c     initialize last obs value when nhit > 0
      lnhobs = inhyr(mindy)

c     initialize last yr/month with inhomog
      call imo2iym(lnhyr,lnhmth,imob)
c     set number of separate paired inhomog changepoints in series
      inhsum = 0

c      print *,' At MLRTest - nobs,yr,mth,rTraw1&2'
c     go through the current hits and the data for both stations
      inhadd = 0
c     inhmo holds the earliest chgpt from each and both stations
c       for aligning with the full period chgpts
      do inh = 1,3
        inhmo(inh) = 9999
      enddo
        
      lhitmo = imob
      do imo = imob, imoe
        call imo2iym(iyr,im,imo)
c       if there are any station inhomog, set toggle and test with last
c         Stn 1 only inhomog: inhadd := 1 
c         Stn 2 only inhomog: inhadd := 2
c         Stn 2 & 3  inhomog: inhadd := 3 
    9   iredo = 0
        if(nhits(itarg,imo).gt.0) then
          if(idebug .ge. 3)  print *,' Stn1 - nhits ', subnet(1),
     *      itarg, iyr, im, nhits(itarg,imo)
          lhitmo = imo
          if(imo .lt. inhmo(1)) inhmo(1) = imo
          if(inhadd .ne. 1 .and. inhadd .ne. 3) inhadd = inhadd + 1
        endif  
        if(nhits(ipair,imo).gt.0) then
          if(idebug .ge. 3) print *,' Stn2 - nhits ', subnet(2),
     *      ipair, iyr, im, nhits(ipair,imo)
          lhitmo = imo
          if(imo .lt. inhmo(2)) inhmo(2) = imo
          if(inhadd .ne. 2 .and. inhadd .ne. 3) inhadd = inhadd + 2
        endif  
        if(inhadd .eq. 3 .and. imo .lt. inhmo(3)) inhmo(3) = imo
        if(imo .eq. imoe) then
          inhadd = 3
        endif  
        
        if(rTraw(imo) .gt. amiss+1.) then
c         increment number in segment
          nobs = nobs+1
c         if there has been a station inhomog
          if(inhadd .gt. 0 .and. nobs .gt. 1) then
c           adjust if the earliest chgpt is in the previous compressed
c             location
            if(inhmo(inhadd) .lt. imo) then
              nobs = nobs - 1
              iredo = 1
            endif  
c           check for more than minlen data since the last inhomog
            if(nobs .lt. mintest) then
c             if not, reset to last and merge stat inhomog type
              if(idebug .ge. 3) 
     *          print *,' Not enough data from: ', lnhyr, lnhmth,
     *            ' to ', iyr, im
              if(lnhadd.eq.1 .and. inhadd.eq.2) inhadd = 3
              if(lnhadd.eq.2 .and. inhadd.eq.1) inhadd = 3
              if(lnhadd.eq.3) inhadd = 3
c             v13c - removed not enough data segments from testseg input
c             backup and set data for testseg to missing
              lobs = lnhobs + 1
              if(mindy .eq. 1) lobs = lnhobs
              print *,' rTraw to miss: ', lobs, imo
              do iobs = lobs, imo
                rTraw(iobs) = amiss
              enddo  
            else
c             keep sum of metadata added
              if(idebug .ge. 3) 
     *          print *,' Enough data from: ', lnhyr, lnhmth,
     *            ' to ', iyr, im
              inhsum = inhsum + 1
c             compile inhomog array for segtest - check against the 
c               current earliest chgpt
              mindy = mindy + 1
c              inhyr(mindy) = imo
              inhyr(mindy) = lhitmo
              lnhobs = inhyr(mindy)
              lnhyr = iyr
              lnhmth = im
            endif 
c           reset the number of obs for current segment
            nobs = 0
c           Keep inhomog type for adjustment allocation
            inhstns(mindy) = inhadd
            lnhadd = inhadd
            inhadd = 0
            do inh = 1, 3
              inhmo(inh) = 9999
            enddo  
c           if current earliest is not current month, check current month
            if(iredo .eq. 1) goto 9
          endif  
        endif
      enddo

      if(idebug .ge. 3) then
        print *,' Inhomog data: mindy,inhyr,yr,mth,inhstns'
        do i = 1, mindy
          indx = inhyr(i)
          call imo2iym(iy,im,indx)
          write(6,'(5i5)') i,indx,iy,im,inhstns(i)
        enddo  
      endif  

c     The yr/mth index arrays have been set for all windows with enough
c       data for the stat tests, now set last yr/mth 
      do imo = imob, imoe
        if(rTraw(imo).gt.amiss+1.) lastmo = imo
      end do  ! end candidate/neighbor loop

c     if needed, readjust the last inhomog at the last paired data value
c      if(lastmo .lt. imoe) then
      if(inhyr(mindy) .lt. lastmo) then
        if(idebug .ge. 3) 
     *    print *,' Last month chg from: ',inhyr(mindy),' to ',lastmo
        inhyr(mindy) = lastmo
        imoe = lastmo
      endif

c     >>>>>>>>>>>>> Split/merge section <<<<<<<<<<<<
      
c    ------ Test segment until there are no changes ------
c        12oct05 - added looping through the BIC after all of 
c        the changes have been found by the Undoc split/merge

c     if there was no changepoints - break first series anyway
      if(mindy .eq. 2) iFirst = 1
      
c     set the station output string
      stnstr = subnet(1) // '-' // subnet(2)

c     initialize number of chgpts in last split/merge testseg calls
      nsplitin2 = 0
      nsplitout = 0
      nmergein2 = 0
      nmergeout = 0

      do inh = 1, ninh
        lsplitin2(inh) = 0
        lmergein2(inh) = 0
        lsplitout(inh) = 0
        lmergeout(inh) = 0
        inhmod(inh) = 0
        offinh(inh) = 0.0
        sseseg(inh) = 0.0
        segslp(1,inh) = 0.0
        segslp(2,inh) = 0.0
      enddo

c     if this is the metadata inclusion method(=3) and neither the UCP
c      nor the histories contain any changepoints, or there where no
c      additional suspect changepoints in the metadata then skip segtest
      if((iFirst .eq. 1 .or. inhsum .eq. 0) .and. method .eq. 3) then
        print *,' skip segtest ifirst, inhsum, method: ',ifirst,inhsum,
     *    method
        iFirst = 0
        ichange = 0
        goto 20
      endif  

c     if this is the undocmented stats call, do merge/split and
c       start with the undocmented stats (not BIC)
      if(method .ne. 3) then
        iMerge = 0
        iBic = 0
        ioptseg = iopt
c       set inqtype for the duration of the looping
        if(iopt .eq. indMD) then
          inqtype = 3
        else if(iopt .eq. indMM) then
          inqtype = 5
        else if(iopt .eq. indXL) then
          inqtype = 4
        else
          print *,' Unknown MLRTest stat: ', iopt
          stop
        endif        
      else  
c       if this is the metadata inclusion method then only merge
c         using the BIC
        iMerge = 1
        iBic = 1
        ioptseg = iopt
        inqtype = 5
      endif  

C      For MINBIC inclusion in method=1
  10  do while (ichange .ne. 0 .or. ibic .eq. 0)
        ichange = 0
        imchange = 0
        ipchange = 0
        ierr = 0
        if(idebug .ge. 3)
     *    write(6,*)' Segtest pass: ',ipass, ' ibic: ',ibic, ' opt: ',
     *      ioptseg, ' imerge: ', iMerge


        if(iMerge .eq. 0) then
          if(idebug .ge. 2)
     *      print *,' Parse segments (isplit = 1), ipass:', ipass
          call testseg(rTraw, ioptseg, inqtype, inhmod, inhstns,
     *      iFirst, mindy, inhyr, nsplitin2, lsplitin2,
     *      offinh, segslp, sseseg, 1, ipchange, subnet, ierr)
c          print *,' ipchange ', ipchange
          if(ierr .ne. 0) goto 20
c         move the prior last split out to the next split input
c         keep this output as the last split out
          do inh = 1, ninh
            lsplitin2(inh) = lsplitout(inh)
            lsplitout(inh) = inhyr(inh)
          enddo
          nsplitin2 = nsplitout
          nsplitout = mindy
        endif  

c       If ifirst pass complete, loop back and test shorter segments
        if(iFirst .eq. 1) then
          iFirst = 0
          goto 10
        endif  
        
        if(idebug .ge. 2)
     *    print *,' Merge segments (isplit = 0), ipass:', ipass
        call testseg(rTraw, ioptseg, inqtype, inhmod, inhstns, 
     *    iFirst, mindy, inhyr, nmergein2, lmergein2,
     *    offinh, segslp, sseseg, 0, imchange, subnet, ierr)
c        print *,' imchange ', imchange
        if(ierr .ne. 0) goto 20
c       move the prior last merge out to the next merge input
c       keep this output as the last merge out
        do inh = 1, ninh
          lmergein2(inh) = lmergeout(inh)
          lmergeout(inh) = inhyr(inh)
        enddo
        nmergein2 = nmergeout
        nmergeout = mindy
        ichange = ipchange + imchange
        ipass = ipass + 1

c       test to see if the undoc stats are finished with the split/merge
        if(ibic .eq. 0 .and. (ichange .eq. 0 .or. ipass .gt. 10)) then
c         reset for the BIC runs
          ichange = 1
          ibic = 1
          ipass = 1
          iMerge = 1
          ioptseg = indKW
        else if(ibic .eq. 1 .and. ipass .le. 10) then
          iochange = 0
c         skip Complementary chgpts model at this time        24 May 06
          goto 15

c         Bringing "complementary chgpts" model earlier into the process
c           to incorporate short divergent sloped and multiple segments. 
c         Due to the ordering of the statistical tests, effectively 
c           making TPR2 the main workhorse

c         initialize divergence trigger (number of interim chgpts)
          intr = 0

          ichg = 2
          do while (ichg .le. mindy-1)
c           count the number of observations in the upcoming segment
            ilen = 0
            do imo = inhyr(ichg)+1, inhyr(ichg+1)
              if(rTraw(imo) .gt. amiss + 1.0) ilen = ilen + 1
            enddo
     
c           Complementary chgpt test triggered by upcoming seglen < minlenshf
c             and continues until upcoming seglen >= minlenshf
            if(ilen .lt. minlenshf) then
              if(intr .eq. 0) then
c               begin of comp chgpt test fragment
                indydiv1 = ichg
                imodiv1 = inhyr(ichg)+1
c               segment immediately before comp chgpt test
                imopre1 = inhyr(ichg-1)+1
                imopre2 = inhyr(ichg)
                idomain = 1
                qsign = offinh(ichg)
              endif
c             current end of comp chgpt test fragment
              imodiv2 = inhyr(ichg+1)
              intr = intr + 1
              if(qsign * offinh(ichg) .lt. 0.0) idomain = 0
            else if(intr .gt. 0) then
c             check last estimate
              if(qsign * offinh(ichg) .lt. 0.0) idomain = 0
              if(idomain .eq. 0) then
c               to get here, short fragments have been accumulated (intr > 0) 
c               and the est adj are on both sides of zero (idomain = 0)
c               and at the end of the short segments (2nd seg len >= minlenshf)
c               set end of the complementary chgpt test
                imodiv2 = inhyr(ichg)
c               segment immediately after comp chgpt test
                imopost1 = inhyr(ichg) + 1
                imopost2 = inhyr(ichg+1)
c               generate series and indices for minbic and erase comp chgpt data
                ix = 0
                mknt = 0
c               process from begin of pre segment to end of post segment
                do imo = imopre1, imopost2
                  ix = ix + 1
                  qx(ix) = imo
c                 catch and keep end of pre segment
                  if(imo .eq. imopre2) iqob1 = ix
c                 keep data outside of the comp chgpt test
                  if(imo .lt. imodiv1 .or. imo .gt. imodiv2) then
                    qy(ix) = rTraw(imo)
                  else
c                   erase data within the comp chgpt test
                    qy(ix) = amiss
                    ndelete(itarg,it2pair,imo) = 'D'
                    rTraw(imo) = amiss
                  endif  
                enddo
                numx = imopost2 - imopre1 + 1
c               go back through MINBIC for finals
                call minbic(2, qx, qy, iqob1, numx, critval, curstat, 
     *            qmin, toff, rmuq, rslpq, rsseq, inqtype, iqtype, 
     *            mknt2, ifail)
                call imo2iym(iyr1, imth1, imopre1)
                call imo2iym(iyrb, imthb, imopre2)
                call imo2iym(iyrb2, imthb2, imopost1)
                call imo2iym(iyr2, imth2, imopost2)
                if(idebug .ge. 1) write(6,'( a, "-", a, " COMP: ", 
     *            2(2(i5, i2.2, i5), " : "), 3f7.2, 2f7.3, 3i5)')
     *            subnet(1),subnet(2), iyr1, imth1, imopre1,
     *            iyrb, imthb, imopre2, iyrb2, imthb2, imopost1, 
     *            iyr2, imth2, imopost2, curstat,  
     *            critval, toff, rslpq, mknt2, iqtype

c               go back to the first comp chgpt index (indydiv1)
                indyback = ichg - indydiv1
c               collapse the breakpoint array
                do jchg = ichg+1, mindy
                  inhyr(jchg - indyback) = inhyr(jchg)
                  inhmod(jchg - indyback) = inhmod(jchg)
                  inhstns(jchg - indyback) = inhstns(jchg)
                  offinh(jchg - indyback) = offinh(jchg)
                  sseseg(jchg - indyback) = sseseg(jchg)
                  segslp(1,jchg - indyback) = segslp(1,jchg)
                  segslp(2,jchg - indyback) = segslp(2,jchg)
                enddo
                mindy = mindy - indyback
                ichg = ichg - indyback

                iochange = iochange + 1
              endif
              intr = 0
            endif
c            print *,' iochange, ichg ',iochange, ichg
            ichg = ichg + 1
          enddo  

   15     ichange = ichange + iochange
        endif  
c        print *,' ichange ', ichange
c       if MINBIC has been turned ON and
c         either no changes have occured or there have been 10 passes
c         then leave testseg loop
        if(ibic .eq. 1 .and. (ichange .eq. 0 .or. ipass .gt. 10))
     *     goto 20  
        
      enddo             

c     >>>>>>>>>>> BIC model vs Split/merge stat section <<<<<<<<<<
   20 nsig = 0

c     initialize output arrays
      do inh = 1, ninh
        isigy1(inh) = 0
        isigy2(inh) = 0
        asigx(inh) = 0.0
        azscr(inh) = 0.0
        mtype(inh) = 0
        istnpt(inh) = 0
      enddo  
      do ichg = 2, mindy-1
        modtype = inhmod(ichg)
        wadj = offinh(ichg)
        end1 = inhyr(ichg)
        do indy = inhyr(ichg) + 1, inhyr(ichg+1)
          if(rTraw(indy) .gt. amiss+1.) then
            beg2 = indy
            goto 30
          endif
        enddo  
   30   continue
        rslp1 = segslp(1,ichg)
        rslp2 = segslp(2,ichg)
        
c       if the segment is "unconfounded" accumulate adjustments
c         and the model type has not been changed to slr0 or slr1
        iprocess = 0
        if(method .eq. 3) then
          iprocess = 1
        else if(inhstns(ichg) .ne. 3 .and. modtype .ge. 3) then
          iprocess = 1
c         test the incoming stattest (inqtype) against the best minbic 
c           data model (modtype). The following compatibilities are assumed:
c            inqtype    modtype
c            3 TPR0     3 TPR0
c            4 TPR1     3 TPR0 & 4 TPR1
c            5 TPR2     3 TPR0, 4 TPR1, 5 TPR2, 6 TPR3 & 7 TPR4
c           any other inqtype-modtype combinations are removed
c          if(inqtype .eq. 3 .and. modtype .ne. 3) then

c         V21F modification - ALLOW ALL TPR MODELS FOR TPR0
          if(inqtype .eq. 3 .and. modtype .lt. 3) then
            iprocess = 0
          else if(inqtype .eq. 4) then
            if(modtype .ne. 3 .and. modtype .ne. 4) iprocess = 0
          else if(inqtype .eq. 5 .and. modtype .lt. 3) then
            iprocess = 0
          endif  
        endif  
        if(iprocess .eq. 0) then
          if(idebug .gt. 2)
     *      print *,'Removed incompatible in/out qtype:',inqtype,modtype
          if(idebug .gt. 1) then
            call imo2iym(iye1,ime1,end1)
            call imo2iym(iyb2,imb2,beg2)
            write(6,'(a,2i5,1x,a," TESTSEG SKIP:",f7.2,2(2i5,i3),2i3)')
     *        stnstr, itarg, ipair, otag, wadj,
     *        end1, iye1, ime1,
     *        beg2, iyb2, imb2, modtype, inhstns(ichg)
          endif
        else
          nsig = nsig + 1
c         year/month of sig changepoint (beg-end of span)
          isigy1(nsig) = end1
          isigy2(nsig) = beg2-1
          asigx(nsig) = wadj
          azscr(nsig) = wadj / sseseg(ichg)
          mtype(nsig) = modtype
          istnpt(nsig) = inhstns(ichg)
          if(asigx(nsig) .gt. 10. .or. asigx(nsig) .lt. -10.) then
            print *, subnet(1),'-',subnet(2),' BIG ASIGX: ', 
     *        asigx(nsig),azscr(nsig),beg1,end1,beg2,end2
          endif  

          if(idebug .ge. 1) then
            call imo2iym(iye1,ime1,end1)
            call imo2iym(iyb2,imb2,beg2)
            write(6,'(a,2i5,1x,a," TESTSEG ADJ: ",2f7.2,2f8.4,2(2i5,i3),
     *       2i3)') stnstr, itarg, ipair, otag, asigx(nsig), 
     *       azscr(nsig), rslp1, rslp2, end1, iye1, ime1,
     *       beg2, iyb2, imb2, modtype, istnpt(nsig)
          endif
          call imo2iym(iy, im, end1)
        endif  
   80 enddo  
      
  100 continue
     
  999 return
      end
c         
c     =======================================================================
c
      subroutine testseg(rTraw, ioptseg, inqtype, inhmod, inhstns,
     *  inFirst, mindy, inhyr, lindy, lnhyr, offinh, segslp, sseseg,
     *  iSplit, ichange, subnet, ierr)
c     loop through the segments - splitting at definitely inhomog,
c       keep when definitely homog, and restraining when questionable
      INCLUDE 'inhomog.parm.mthly.incl' 
      INCLUDE 'inhomog.comm.mthly.incl'
      INCLUDE 'inhomog.restart.mthly.incl'
      INCLUDE 'inhomog.MDparm.mthly.incl'

c      = cand & net temperature data, 12 monthly + annual
c     ioptseg is the Technique number (ref ser + stat test) 
c            See: inhomog.comm.mthly.incl for definitions
c     inqtype is the BIC model Q type for the inhomog splitting algorithm
c            (set as indMD=3,indMM=5,indXL=4) for minbic
c     inhyr is the inhomogeneous breakpoint array 
c            (including begin/end of series)
c     mindy is the total number of date indices in the chgpt array
c     iFirst indicates whether this is the Initial test of the series
c            currently this is used to force a chgpt at the peak whether
c            significant or not... to test the shorter series
c     iSplit indicates whether this is
c            Annual Merging call = 0
c            Annual Parsing call = 1
c            Monthly Merging call = 2
c            Monthly Parsing call = 3
c     ichange (output) indicates whether any change occured in the
c            segments
        
      integer inhyr(ninh), lnhyr(ninh), inhmod(ninh), inhstns(ninh),
     *   iexyr(ninh), iexmod(ninh), iexstns(ninh)
      real offinh(ninh), offiex(ninh),segslp(2,ninh),segsex(2,ninh),
     *  sseseg(ninh), ssesex(ninh)

c     raw input data with different array dimensions
      dimension rTraw(nmo)
      
c     difference between cand & composite neighbors (q) and 
c        standardized (z) series
      dimension qx(nmo),qy(nmo),z(nmo)
      
c     T-statistic and 9-point Binomial Average series derived from z
      dimension ts(nmo)

c      dimension tstat(kmo)
      integer mknt, mknt2(2)
      
      integer ounit

c     temp arrays for data transfer
      dimension rDatOut(nmo,nStns)

c     inhyr is the "inhomogeneous year" (chgpt) array, that is, the
c            segment array (including endpoints)
c     mindy is the total number of date indices in the chgpt array
      character*2 otag
      integer job
      character*6 subnet(nstns)
      character*9 limsense
      
      real acan(numyr), aref(numyr)
      integer ipyr(numyr), ifind(numyr)
      integer StatTest, end1
      
      real rmuq(2), rslpq(2)
      mintest = 5
c      mintest = minlen

      limsense = ' limit>: '

      job = 2
c     ioptseg can be negative at the end of the split/merge technique
c       to ONLY estimate the amplitudes
      iopt = abs(ioptseg)

c     set the outtag and out unit for this iopt version
      otag = c2tech(iopt)
      ounit = bunit + iopt

c     always start out with the first segment
      indy1 = 1
      if(isplit .eq. 1) then
        indy2 = 2
      else
        indy2 = 3
      endif    

      ichange = 0
      
c     Need a temporary variable that does not change the value in 
c     the calling program
      iFirst = inFirst
            
c     initialize inhomog expansion series array and indices
c     these are working arrays which are copied back to the incoming
c     arrays at the end of the split process
      if(isplit .eq. 1) then
        mexpn = mindy
        do ichg = 1, mindy
          iexyr(ichg) = inhyr(ichg)
          iexmod(ichg) = inhmod(ichg)
          iexstns(ichg) = inhstns(ichg)
          offiex(ichg) = offinh(ichg)
          ssesex(ichg) = sseseg(ichg)
          segsex(1,ichg) = segslp(1,ichg)
          segsex(2,ichg) = segslp(2,ichg)
        enddo 
      endif  
      iexy1 = indy1
      iexy2 = indy2
      
      lndy1 = 1
      lob1 = lnhyr(lndy1)

c ------    loop through all of the incoming segments (parse)  ------
c ------                or adjacent incoming segments (merge)  ------
      do while (indy2 .le. mindy)
      
        StatTest = inhomog
        if(indy1 .eq. 1) then
          iob1 = inhyr(indy1)
        else  
          iob1 = inhyr(indy1)+1
        endif  
        iob2 = inhyr(indy2)
c       if merge call - save interior chgpt
        if(isplit .eq. 0) iyrstat = inhyr(indy1 + 1)

c       generate yr/mth at beg, end, and brkpt of segments 
        call imo2iym(iyr1, imth1, iob1)
        call imo2iym(iyr2, imth2, iob2)
        if(ioptseg .lt. 0 .or. iscrit(iopt) .eq. iBstat) goto 10

c       is this segment stable (that is, is the segment the same as the 
c        last time testseg was called)
        do lndy1 = 1, lindy-1
          if(inhyr(indy1) .eq. lnhyr(lndy1)) then
            do inc = 1, indy2-indy1
              if(inhyr(indy1+inc) .ne. lnhyr(lndy1+inc)) goto 10
            enddo
            if(idebug .ge. 2)
     *        print *,'Stable segment: ',iob1,iyr1,imth1,iob2,iyr2,imth2
            goto 800
          endif
        enddo

c       initialize final stat series
   10   do im = 1, nmo
          qx(im) = amiss
          qy(im) = amiss
          z(im) = amiss
          ts(im) = 0.0
        enddo

c ------ Calculate/Recalculate ref series for current segment --------
c       estimate variance of each monthly time series
c       - also works for annual
c       mknt = number of non-missing months between iob1 & iob2
c       numx = number of serial months between iob1 & iob2
        mknt = 0
        ix = 0
        do nx = iob1,iob2
        ix = ix + 1
        qx(ix) = nx
        qy(ix) = rTraw(nx)
          if(rTraw(nx).gt.amiss+1.) mknt = mknt + 1
        enddo
        numx = iob2 - iob1 + 1      
        
        if(mknt .lt. mintest) then
          if(idebug .ge. 2)
     *      write(6,'(a, " ", a, " Skip 0 Segment too Short ", i4,
     *        i3, " to ",i4,i3,2i5)') subnet(1), otag, iyr1, imth1,
     *        iyr2, imth2, mknt, igood
c         Do nothing - goto next segment
          goto 800
        endif  
     
c       standardize series for split/merge stats
        if(iscrit(iopt) .ne. iBstat) then
          call standard(qy,z,1,numx,nmo)
        else
c         for estamt use full temperature series
          do ix = 1, numx
            z(ix) = qy(ix)
          enddo
        endif
              
c  ------ Generate test statistic for shift or step-change ----------
c       pass the window interval instead of kmo
c       see inhomog.comm.mthly.incl for ioption definitions
c          MLRTest.v17g modification - make minbic offset estimation call 
c          (iopt = 2) common for all of the split/merge calls (not indKW)
        if(ioptseg .gt. 0 .and. iscrit(iopt) .lt. iBstat) then
          if(iscrit(iopt) .eq. iTstat) then
c           generate statistic (and series for chgpt) for max-likely
            call snits(kmo,z,ts,mknt)
            iqtype = 3  
          else if(iscrit(iopt) .eq. iFstat) then
c           generate statistic (and series for chgpt) for 2-phase
            call twophreg(kmo,z,ts,mknt)
            iqtype = 5
          else if(iscrit(iopt) .eq. iXstat) then
c           generate statistic (and series for chgpt) for 2-phase w/const slope
            call twophreg1(kmo,z,ts,mknt)
            iqtype = 4
          endif

c    -------- Evaluate test statistic wrt threshold ---------------------
          iPeak = 0
          rPeak = 0.0
          do iknt = 2, numx - 1
            if(ts(iknt) .gt. rPeak) then
              iPeak = iknt
              rPeak = ts(iknt)
            endif  
          enddo  
                 
          if(iscrit(iopt) .eq. iXstat) then
c           Two Phase w/const slope critical value
            call tpr1table(mknt, critval)
          else  
c           Call lookup to get the 90%, 95%, and 97.5% critical values for
c            Tmax or Fmax test for a full segment
            call lookup(mknt,crit90,crit95,crit975,iscrit(iopt))
            critval = crit95
          endif  

          if(iPeak .le. 0) then
            write(6,'(a," ",a,"            No found peaks ",i4,i3,
     *        " to ",i4,i3)')subnet(1),otag,iyr1,imth1,iyr2,imth2
            StatTest = homog
            if(isplit .eq. 0) then
              if(idebug .ge. 2) 
     *         write(6,'(a," ",a," Compress 1 out peak at ",2f7.2," ")')
     *          subnet(1), otag, iyrb, imthb
c             for merging, this is a change (collapse)
              goto 200
            else
c             for parsing, this is NOT a change
              goto 800
            endif  
          endif  
     
c         test the peak stat against the critical value
          if(rPeak.lt.critval) StatTest = homog
          curstat = rPeak
          end1 = qx(iPeak)
          call imo2iym(iyrb, imthb, end1)

c         Fragment First if either homog or inhomog
          if(iFirst .eq. 1) then
c           force first to split if homog, to test shorter segments
            if(StatTest .eq. homog) StatTest = inhomog
            if(idebug .ge. 2) write(6,'(a, "-", a, " ", a, 
     *        "             FIRST series ", i4, i3, " to ", i4, i3,
     *        " | at ", i4, i3, " ts: ", f7.2, a9, f6.2)')
     *        subnet(1), subnet(2), otag, iyr1, imth1, iyr2, 
     *        imth2, iyrb, imthb, curstat, limsense, critval
            goto 300
          else if(isplit .eq. 0) then
            if(StatTest .eq. homog)then
c             in a merge pass, collapse a homog chgpt
              if(idebug .ge. 2) write(6,'(a, "-", a, " ", a,
     *          " Compress 2 out peak at  ", i4, i3, " | ",
     *          f7.2, a, f7.2, " ")') subnet(1), subnet(2),
     *          otag, iyrb, imthb, curstat, limsense, critval
              goto 200
            else  
c             in a merge pass, leave an inhomog segment alone
              if(idebug .ge. 2) write(6,'(a, "-", a, " ",a,
     *          " Peak kept in merge at   ", i4, i3, " | ", 
     *          " ts: ", f7.2, a9, f6.2)')subnet(1), subnet(2),
     *          otag, iyrb, imthb, curstat, limsense, critval
              goto 800
            endif
          else
            if(StatTest .eq. homog) then
c             in a split pass, leave a homog segment alone
              if(idebug .ge. 2) write(6,'(a, "-", a, " ", a, 
     *          "       Homogeneous series ", i4, i3, " to ", i4, i3,
     *          " | at ", i4, i3, " ts: ", f7.2, a9, f6.2)')
     *          subnet(1), subnet(2), otag, iyr1, imth1, iyr2,
     *          imth2, iyrb, imthb, curstat, limsense, critval
              goto 800
            else
c             in a split pass, fragment an inhomog segment
              if(idebug .ge. 2) write(6,'(a, "-", a, " ", a,
     *          " Inhomogeneity for series ", i4, i3, " to ", i4, i3, 
     *          " | at ", i4, i3, " ts: ", f7.2, a9, f6.2)') 
     *          subnet(1), subnet(2), otag, iyr1, imth1, iyr2, 
     *          imth2, iyrb, imthb, curstat, limsense, critval
              goto 300
            endif  
          endif  

C         This is the ONLY parsing change point! These must NOT be
c         added back into the orginal chgpt array, but kept in an
c         array of their own
  300     iFirst = 0
          ichange = ichange + 1

c         update expansion year series array and indices
          mexpn = mexpn + 1
          do ichg = mexpn, iexy2+1, -1
            iexyr(ichg) = iexyr(ichg-1)
            iexmod(ichg) = iexmod(ichg-1)
            iexstns(ichg) = iexstns(ichg-1)
          enddo 
c         set the end of the series to the highest peak
          iexyr(iexy2) = end1
          iexmod(iexy2) = iqtype
          iexstns(iexy2) = 0
          iexy2 = iexy2 + 1
          iexy1 = iexy1 + 1
          goto 800
        
        else if(ioptseg .lt. 0 .or. iscrit(iopt) .eq. iBstat) then
        
c     -------------------------  run the BIC  -----------------------
          end1 = inhyr(indy1+1)
          call imo2iym(iyrb, imthb, end1)
          do ix = 1, numx
            if(end1 .lt. qx(ix)) then
              iqob1 = ix - 1
              goto 30
            endif
          enddo
        
c         determine the optimum Bayesian Info Criteria for found station/chgpt
   30     iyrstat = iqob1
          if(idebug .ge. 2)
     *      print *,' Entering MINBIC:',iyr1,imth1,iyrb,imthb,iyr2,imth2
          call minbic(2, qx, qy, iqob1, numx, critval, curstat,
     *      qmin, toff, rmuq, rslpq, rsseq, inqtype, iqtype, mknt2, 
     *      ifail)

c ------- tests for "break slope" model inconclusive, see versions before
c            v20e for brkpt model (=8)

c         save BIC data
          inhmod(indy1+1) = iqtype
          offinh(indy1+1) = toff
          sseseg(indy1+1) = rsseq
          segslp(1,indy1+1) = rslpq(1)
          segslp(2,indy1+1) = rslpq(2)
          call imo2iym(iyrb2, imthb2, end1 + 1)
          if(idebug .ge. 2) write(6,'( a, "-", a, " BIC: ", 
     *      2(2(i5, i2.2, i5), " : "), 3f7.2, 2f7.3, 3i5)')
     *      subnet(1),subnet(2), iyr1, imth1, iob1, iyrb, imthb, end1, 
     *      iyrb2, imthb2, end1 + 1, iyr2, imth2, iob2, curstat,  
     *      critval, toff, rslpq, mknt2, iqtype
            
          if(iqtype .lt. 3) then
            StatTest = homog
c           compress MINBIC segment here ONLY IF an SLR model is returned
c            goto 200
c           DO NOT collapse chgpt if SLR model - v21b
            goto 800
          else
            StatTest = inhomog
            goto 800
          endif

        else
          write(6,*)' Unknown statistical option in MLRTest, stopping'
          stop  
        endif    
     
c       collapse the breakpoint array
  200   do ichg = indy1+1, mindy-1
          inhyr(ichg) = inhyr(ichg+1)
          inhmod(ichg) = inhmod(ichg+1)
          inhstns(ichg) = inhstns(ichg+1)
          offinh(ichg) = offinh(ichg+1)
          sseseg(ichg) = sseseg(ichg+1)
          segslp(1,ichg) = segslp(1,ichg+1)
          segslp(2,ichg) = segslp(2,ichg+1)
        enddo
        mindy = mindy - 1
        indy2 = indy2 - 1
        indy1 = indy1 - 1

        ichange = ichange + 1
          
c ----- Update first and last chgpt date pointers for next segment -------
  800   indy2 = indy2 + 1
        indy1 = indy1 + 1
        iexy2 = iexy2 + 1
        iexy1 = iexy1 + 1
  810 enddo  ! end of series segmenting loop

  998 if(isplit .eq. 1 .and. ichange .gt. 0) then
c       repopulate the incoming array if this is a split pass with changes
        mindy = mexpn
        do ichg = 1, mexpn
          inhyr(ichg) = iexyr(ichg)
          inhmod(ichg) = iexmod(ichg)
          inhstns(ichg) = iexstns(ichg)
          offinh(ichg) = offiex(ichg)
          sseseg(ichg) = ssesex(ichg)
          segslp(1,ichg) = segsex(1,ichg)
          segslp(2,ichg) = segsex(2,ichg)
        enddo 
      endif          

  999 return
      end

c     =======================================================================
      subroutine standard(rRaw,rStd,nx1,nx2,iTotlObs)
c-----------------------------------------------------------------------
c     Standardize rRaw from nx1 to nx2, output results in rStd
c-----------------------------------------------------------------------      
      INCLUDE 'inhomog.parm.mthly.incl'  
        
      dimension rRaw(iTotlObs), rStd(iTotlObs)

      rSum = 0.0
      rNum = 0.0
      do i=nx1,nx2
        if(rRaw(i) .ne. amiss) then
          rSum = rRaw(i) + rSum
          rNum = rNum + 1
        endif  
      enddo
      
      rMean = rSum / rNum
      rVarSumm = 0.0
      do i = nx1,nx2
        if(rRaw(i) .ne. amiss) rVarSumm = rVarSumm + (rRaw(i)-rMean)**2
      enddo
      
      rsqVar = sqrt(rVarSumm / (rNum-2))
      
      do i = nx1,nx2
        if(rRaw(i) .ne. amiss) then
          rStd(i) = (rRaw(i) - rMean) / rsqVar
        else
          rStd(i) = amiss
        endif    
      enddo
                
      return
      end

C     ******************************************************************
C     *  Subroutine acresid.
C     *  This subroutine obtains the autocorrelations of the residuals
C     ******************************************************************

      subroutine acresid (nr,resid,dwstat,alag1)
      INCLUDE 'inhomog.parm.mthly.incl'  

      dimension ac(2)
      real resid(nmo)
      real*8 z(nmo)

      do 100 i=1,nr
        z(i) = resid(i)
  100 continue
  
c     Use Claude Duchon's lag autocorrelation routine
      call uacovf(z, nr, amean, var, ac, 2)
      alag1 = ac(1) / var
 
      return
      end
 
