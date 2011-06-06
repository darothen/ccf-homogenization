c***********************************************************************
c                          program ushcn_corr_2004.f
c                      from the ushcn program hcn_corr.f
c 
c version     modification                                     date
c -------     ------------                                  ----------
c
c            minor change - added toggle to turn off minimum
c              month checking (icheck)                        09jan09
c
c             minor change station number to character        24sep08
c
c  2004.v3    three limiting parameters specific the network 
c               from input through output
c             ndist - number nearest-neighbor in input dist file
c             minpair (minsta*2) - min pair per year/month
c             nstns - number of stations in output corr file
c                     (also min stn limit for neighborhood)
c                                                           15 jul 2005
c
c  2004       catch the 100 nearest neighbors               6 apr 2005
c               with monthly anomaly correlations 
c
c  2000.v2    major change to sort by correlation of first-difference
c               so that step changes do not dramatically reduce the
c               correlation of the stations                 19 sep 03
c
c  2000.v1c   resurrected seasonal correlation calculations
c                  for areal edits                          21 apr 2003 cw
c             minsea upped to 20 due to canada urban task   14 aug 2003 

c  2000.v1b   numcor fixed and parameters globalized with
c                  rest of pre-homogeneity codes            15 jan 2003 cw

c  2000       added correlations to annual output files

c   1.1       sun workstation version                       10/24/1995
c 
c   1.0       original
c 
c description:  this program calculates the 10/20 highest correlated
c               stations of the 40 nearest neighbors.  the highest 20 
c               are calculated for temperature while either the 
c               highest 10 or 20 are calculated for precipitation.
c               the major difference between this version and the 
c               hcn version are 
c             1) this routine uses separate lists for the candidates 
c                   and the network neighbors
c             2) reads from lat,lon lists instead of the ushcn station
c                   history file
c             3) ingests the data from separate files instead of one
c                   big, merged file

      program hcncor
      include 'prehomog.parm.incl'
      include 'prehomog.comm.incl'  
      include 'prehomog.corr.incl'

c     parameter : description
c      ibegyr   : begin year of correlation
c      ilstyr   : last year of correlation
c      nstns    : number of neighbors in final subnet 
c      nstns-1   : number of correlations
c      maxnets  : maximum number of potential subnetworks (cand+neighs)
c      minmths  : minimum number of cadidate w/network station pairs for
c                 correlation

      integer id(maxnets), netid(ndist)
      real bad
      real c1(numyr*12), n1(numyr*12)
      real c1diff(numyr*12), n1diff(numyr*12)

      integer iargc, narg
      character*132 argv
     
      character*132 ofile/''/
      character*6 cstn,cid(maxnets)
      real corrs(ndist),var(ndist),corrn(ndist)
      real corsq(maxnets,12),stdmn(maxnets,12)
      integer nim(12)
      integer pairs(ndist),ptr2(ndist)
      integer ptr(maxnets,ndist),netptr(ndist)
      integer igood/0/, iskip/0/,inone/0/

c     file names for the network stations and the candidate
      character*80 cfname
      
c     output arrays
      integer idout(ndist),ptrout(ndist)
      real corrout(ndist)
      
c     work arrays for minimizing # neighbors
      integer ksum(ibegyr:ilstyr, 12), lowtoo(ibegyr:ilstyr, 12),
     *  iremove(ndist), jsum(ibegyr:ilstyr, 12)
      
      integer netid0/0/
      real corrs0/0.0/

c     initialize variables
c     number of stations to correlate comes from 

      bad = amiss
      iwrite = 0
      idebug = 0
      
c     default to HCN input data type
      indatyp = 0
      nc = 0
            
c     first two parameters are required and are decoded as follows:
c       four digit year (>= ibegyr+minmths)
c       elements to process (1=max, 2=min, 3=avg, 4=pcp)
      print *,' From prehomog include headers:'
      print *,' Minmths: ', minmths
      print *,' Maxnets: ', maxnets
      print *,' Minstns: ', minstns
      print *,' Nstns:   ', nstns

c     begin command line code      
c     get the number of command line arguments
      narg = iargc()
      if(narg .lt. 3) then
        print *,' not all required parameters given'
        goto 10
      endif
      
      metafile = ''
      outfile = ''
      distfile = ''
      canddir = ''
      netdir = ''
      cdtype = ''
      
c     Default minimum month checking (<8) is ON
      icheck = 1
      
      iarg = 1

c     get last year of data to process
      call getarg(iarg,argv)
      iarg = iarg + 1
      read(argv,fmt='(i4)') endyr
      if(endyr .lt. ibegyr + mincorr) then
        print *, 'command line input endyr:', endyr, 
     *    ' too low, must be at least:', ibegyr+mincorr
        stop
      endif
      print *, ' highest year to process: ', endyr

c     get element(s) to process
      call getarg(iarg,argv)
      iarg = iarg + 1
      read(argv,fmt='(i1)') inel
      if(inel .lt. 1 .or. inel .gt. maxelem) then
        print *, ' element must be 1 to ', maxelem
        stop
      endif  
      print *,' processing element: ', celem(inel)

c     go thru remaining command line arguments - keep the valid ones
c        and ignoring the undefined ones.
      do while(iarg .le. narg)
        call getarg(iarg, argv)
        iarg = iarg + 1
        if(argv .eq. '-d') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          distfile = argv
          print *, ' candidate neighbor distance input file :',distfile
        else if(argv .eq. '-m') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          metafile = argv
          print *, ' candidate neighbor meta input file :', metafile
        else if(argv .eq. '-o') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          outfile = argv
          print *, ' cand-neigh correlation output files :', outfile
        else if(argv .eq. '-C') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          canddir = argv
          print *, ' base candidate monthly data directory:', canddir
        else if(argv .eq. '-N') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          netdir = argv
          print *, ' base network monthly data directory:', netdir
        else if(argv .eq. '-p') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          cdtype = argv
          print *, ' process stage for correlation data :', cdtype
        else if(argv .eq. '-u') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i4)') nc
          print *, ' candidate-neighbor merge - nc:',nc
        else if(argv .eq. '-t') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i1)') indatyp
          print *,' Input data formal type (0=hcn;1=3220): ',indatyp
        else if(argv .eq. '-X') then
          icheck = 0
          print *,' Minimum months checking is turned off'
        else
          print *, ' unknown argument :', argv, ': skipping'
        endif
      enddo
      
      if(metafile.ne.'' .and. outfile.ne.'' 
     *   .and. distfile.ne.'' .and. cdtype.ne.'' .and. canddir.ne.'' 
     *   .and. netdir.ne.'') goto 11
      
   10 print *,' normals correlation network generation routine'
      print *, ' usage: normals_corr_2000 year elem ',
     *  '-d distfile -m metafile -C canddir -O netdir -o outfile ',
     *  '-p process'
      print *,'   required parameters'
      print *,'      year         four digit year'
      print *,'      elem         1=max, 2=min, 3=avg, 4=pcp'
      print *,'     -u nc         number of candidate stations'
      print *,'     -p cdtype     data process level for corel'
      print *,'    input files/directories ---'
      print *,'     -d distfile   candidate neighbor distance file'
      print *,'     -m metafile   candidate neighbor stations file'
      print *,'     -o outfile    output filename (base)'
      print *,''
      print *,'     -C canddir    base candidate data directory'      
      print *,'     -N netdir     base network data directory'
      print *,' note: raw data files assumed to be parsed into '
      print *,'        sudirectories by state (01,02...,50)'
      stop

   11 continue
      if(nc .eq. 0) print *, ' candidate-neighbor merge disabled'
      if(icheck .eq. 1) print *,' Minimum month checking in ON'

c     set the last non-blank character in the base file directory
      icdlen = lnblnk(canddir)
      ictlen = lnblnk(cdtype)
      if(canddir(icdlen:icdlen).ne.'/') then
        canddir(icdlen+1:icdlen+1)='/'
        icdlen = icdlen + 1
      endif  
      print *,icdlen, ':', canddir(1:icdlen), ':'
      indlen = lnblnk(netdir)
      if(netdir(indlen:indlen).ne.'/') then
        netdir(indlen+1:indlen+1)='/'
        indlen = indlen + 1
      endif  
      print *, indlen, ':', netdir(1:indlen), ':'
      
c     open, read and close merged cand-network station metafile
      open(10,file=metafile,status='old',err=150)

      do inet = 1,maxnets
        read(10,20,err=170,end=35) id(inet)
   20   format(i6)
        write(cid(inet), '(i6.6)') id(inet)

c       fill out the pathname with the appropriate directory early
c       nc indicates the number of candidates in the metafile from
c       the command line
        if(indatyp .eq. 0) then
c         HCN default input
          if(inet .le. nc) then
            nfname(inet) = canddir(1:icdlen) // cdtype(1:ictlen) // 
     *        '/' // cid(inet) // celem(inel) // '.' // cdtype
          else
            nfname(inet) = netdir(1:indlen) // cdtype(1:ictlen) //
     *        '/' // cid(inet) // celem(inel) // '.' // cdtype
          endif
        else  
c         HOFN input - all the files are in one directory as stn".3220"
          nfname(inet) = canddir(1:icdlen) // cid(inet) // '.3220'
        endif
      enddo

   35 close(10)
      nnets = inet - 1
      print *,' Number of networks: ', nnets

c     open annual output file
      ofile = outfile(1:lnblnk(outfile))
      open(20,file=ofile,err=210)

c     get data for all network reference stations (hcn) from 
c       begin year to present for current season (english units)
      print *,' reading in data'
      do j = 1,nnets
        call getdata(j, indatyp)
      enddo

c     open candidate station input file
      open(10,file=distfile,status='old',err=176)

c     get correlations for nearest neighbors, sort them, and then
c     write along with station ids
      do j = 1, nnets
c       reading the hcn-cdmp merged station list
        read(10,*,err=174,end=174) icid
c       generate the candidate data file name
        write(cstn, '(i6.6)') icid
        cfname = nfname(j)
        if(mod(j,100) .eq. 1) print *,j,' Correlated: ',cstn

c       initialize stat output variables
        do im = 1, 12
          corsq(j,im) = 0.
          stdmn(j,im) = 0.
        enddo  

c       now get the pointers to the network for the candidate 

   24   read(10,25,err=174,end=174) (ptr(j,k),k=1,ndist)
   25   format(<ndist>(i6,1x),/)
        nid = ptr(j,1)
        
        if(id(nid) .ne. icid) then
          print *,' lost synchronization between dist and meta indices'
          goto 130
        endif  

c       populate the current candidates data arrays
        do iy = ibegyr, ilstyr
          do im = 1, 12
            cdata(iy,im) = rdata(nid,iy,im)
          enddo
        enddo    

c       have all the necessary data - compute correlation with networ
c       initialize first cell of correlation and pointer arrays
        corrs(1) = 1.
        ptr2(1) = ptr(j,1)
        corrn(1) = 1.
          
        do k = 2,ndist
c         initialize correlation and pointer arrays
          corrs(k) = 0.
          ptr2(k) = ptr(j,k)
          pairs(k) = 0
          
c         align non-missing years between candidate and neighbor
          do im = 1, 12
            nim(im) = 0
          enddo
          nal = 0
          ifirst = 0
          ilast = 0  
          do iy = ibegyr, ilstyr
            do im = 1, 12
              if(cdata(iy,im) .ne. amiss .and. 
     *          rdata(ptr(j,k),iy,im) .ne. amiss) then
                if(ifirst .eq. 0) ifirst = iy
                ilast = iy
                nim(im) = nim(im) + 1
                nal = nal + 1
                c1(nal) = cdata(iy,im)
                n1(nal) = rdata(ptr(j,k),iy,im)
              endif
            enddo  
          enddo    
            
          do im = 1, 12
            if(nim(im) .lt. minmths) then
              corrs(k) = 0.
              corrn(k) = 0.
              ptr2(k) = 0
              if(idebug .eq. 1) then
                print *,nid,' skipping pair: ',cid(nid),' - ', 
     *            cid(ptr(j,k)),' im= ',im,' nim= ', nim(im)
              endif
              goto 70
            endif
          enddo
            
c         make the first difference of the candidate station
          call frstdif(nal, c1, c1diff, bad)

c         make the first difference of each network station series
          call frstdif(nal, n1, n1diff, bad)

c         calculate correlation
          call correl(nal, c1diff, n1diff, amiss, corrs(k), 
     *      var(k), var(1), pairs(k))
          if(pairs(k) .lt. minmths) then
            print *,cid(nid),'-',cid(ptr(j,k)),' too few corr-pairs ',
     *        nid,ptr(j,k),pairs(k),corrs(k)
          endif  

          corrn(k) = corrs(k)
   70   enddo

c       sort 10/20 best correlated stations
        call sortco(corrn,ndist,ndist,ptr2)

c       attempt to equalize number of neighbors over period of record

c       summ all of the neighbor year/months
        minpair = minstns * 2
        kstns = 0
        do iy = ibegyr, ilstyr
          do im = 1, 12
            ksum(iy,im) = 0
            jsum(iy,im) = 0
            lowtoo(iy,im) = 0
            if(cdata(iy,im) .ne. amiss) then
              do k = 2, ndist
                if(corrn(k) .gt. corrlim) then
                  if(rdata(ptr2(k),iy,im) .ne. amiss)
     *              ksum(iy,im) = ksum(iy,im) + 1
                else
c                 do not use neighbors below correlation threshold
                  goto 75
                endif
              enddo
   75         kstns = k - 1
              jsum(iy,im) = ksum(iy,im)
c             detect dates when not enough neighbors
              if(ksum(iy,im) .lt. minpair) then
                if(idebug .eq. 1)
     *            print *,' Total less than minpair: ',cid(nid),iy,im
                lowtoo(iy,im) = 1
              endif
            endif
          enddo
        enddo

c       if total number of neighbors is greater than number desired
c       keep data in low info periods as we reduce the number of neighbors
c       from least to highest corr
        jstns = kstns
        if(kstns .le. nstns) goto 90
        do k = jstns, 2, -1
          kptr = ptr2(k)
          iremove(k) = 1
          npair = 0
          do iy = ibegyr, ilstyr
            do im = 1, 12
              if(cdata(iy,im) .ne. amiss .and. 
     *          rdata(kptr,iy,im) .ne. amiss) then
                npair = npair + 1
c               if current neighbor contributes data in a low
c                 information time period - it cannot be removed
                if(ksum(iy,im) .le. minpair) then
                  if(idebug .eq. 1)
     *              print *,' Cannot remove: ',cid(nid), '-',  
     *                cid(kptr), iy, im, ksum(iy,im), lowtoo(iy,im)
                  iremove(k) = 0
                  goto 85
                endif  
              endif  
            enddo
          enddo  
c         remove the station if possible
   85     if(iremove(k) .eq. 1) then
c           when the number of neighbors at ANY year/month drops below
c           the maximum - discontinue reduction
            if(kstns-1 .lt. nstns) goto 90
            if(idebug .eq. 1)
     *        print *,' Remove: ',cid(nid),'-',cid(kptr),npair,corrn(k)
            corrn(k) = 0.0
            kstns = kstns - 1
c           update neighbor data count for period of record
            do iy = ibegyr, ilstyr
              do im = 1, 12
                if(cdata(iy,im) .ne. amiss .and. 
     *            rdata(kptr,iy,im) .ne. amiss) then
                  ksum(iy,im) = ksum(iy,im) - 1
                endif
              enddo
            enddo
          endif
        enddo
        
   90   if(idebug .eq. 1) then
          do iy = ibegyr, ilstyr
            do im = 1, 12
              if(jsum(iy,im) .gt. 0) 
     *          print *,'Original-Final: ', cid(nid), iy, im, 
     *            jsum(iy,im), ksum(iy,im)
            enddo
          enddo  
        endif  

        print *,'Orig-Final Number stns: ', cid(nid), ' ', jstns, kstns

c       sort 10/20 best correlated stations
c        must go through the entire ndist elements to get the nstns network!
        call sortco(corrn,ndist,ndist,ptr2)

c       write out best network
        do k = 1, nstns
          if(corrn(k) .gt. corrlim) then
            idout(k) = id(ptr2(k))
            ptrout(k) = ptr2(k)
            corrout(k) = corrn(k)
          else
c            if(ptr2(k) .ne. 0)
c     *        print *,cid(nid),'-',cid(ptr2(k)),' too low corr ',
c     *        nid,ptr2(k),corrn(k)
            idout(k) = 0
            ptrout(k) = 0
            corrout(k) = 0.0
          endif
        enddo      
            
        write(20,110) (idout(k),k=1,nstns),(ptrout(k),k=1,nstns),
     *      (corrout(k),k=1,nstns)
  110   format(<nstns>(i6.6,1x),/,<nstns>(i6,1x),/,<nstns>(f6.2,1x))
        iwrite = iwrite + 1

      enddo
  
c     close data file
  130 close(20)

c     skip error messages
      go to 310

c     print error message and abort program
  150 print *,'cannot open network neighborhood input file :', metafile
      go to 290

  170 print *,'cannot read network neighborhood file.', metafile
      go to 290

  174 print *,'cannot read cand-network file:', distfile
      go to 290
      
  176 print *, 'cannot open distance file:', distfile
      go to 290

  210 print *, 'cannot open output file.', ofile
      go to 290

  250 write(6,260)
  260 format('cannot close output file.')

  290 write(6,300)
  300 format('program aborting.')

  310 stop

      end

c***********************************************************************
c end of program hcncor.                                               *
c***********************************************************************


c***********************************************************************
c                          subroutine getdata                           *
c                                                                      *
c   2004.v2   linux workstation version returns monthly anomalies
c               for 1st diff correlation in main
c                                                        20 apr 2005 cw
c
c   1.1       sun workstation version                       12/27/1994 *
c                                                                      *
c   1.0       original                                                 *
c                                                                      *
c author:  claude williams                                             *
c                                                                      *
c modifiers:  pam hughes                                               *
c             david bowman                                             *
c                                                                      *
c usage:  getdata(integer endyr, integer stnind)                        *
c                                                                      *
c description:  this subroutine fetches monthly for one hcn station.   *
c                                                                      *
c***********************************************************************

c ---- fix so that if the stnid < nc (number of candidates) then
c ----  read from candidate data directory
c ---- else
c ----  read from neighbor data directory

      subroutine getdata(stnind, indatyp)
      include 'prehomog.parm.incl'  
      include 'prehomog.comm.incl'  
      include 'prehomog.corr.incl'  

c     indatyp indicates the input data type and units
c       0 = HCN data format and English
c       1 = 3220 data format and English
c       2 = GHCN data format and Metric

c     local variables
      integer stnind, idata(12)
      integer dyear
      integer nummth(12)
      integer imiss /-9999/
      
c     3220 input variables
      character*6 dstn
      character*3 cRecType
      character*4 cElmType
      integer iMonth(12),iDay(12),iValue(12)
      
      real ddata(ibegyr:ilstyr,12),tdata(12)
      real summth(12)
      real adiv(2) /10.,100./

c     initialize variables
      dstn = "      "
      iflag = 0

c     initialize data array for current station
      do iy = ibegyr,ilstyr
        do im = 1,12
          ddata(iy,im) = amiss
          rdata(stnind,iy,im) = amiss
        enddo  
      enddo

c     open station data file
      print *, nfname(stnind)(1:lnblnk(nfname(stnind)))
      open(17, file=nfname(stnind), err=110)
    
      do im = 1,12
        summth(im) = 0.0
        nummth(im) = 0
      enddo  

      do while(1 .eq. 1)
c       get data values and day flags
        if(indatyp .eq. 0) then
c         candidate station data in normals format          
          read(17,1000,end=100,err=110) dstn,dyear,(idata(i),i=1,12)
 1000     format(a6,1x,i4,12(i6,1x))
          if(dyear .lt. ibegyr .or. dyear .gt. ilstyr) go to 95
          if(inel .ne. 4) then
            rdiv = adiv(1)
          else
            rdiv = adiv(2)
          endif    
          do im = 1, 12
            if(idata(im) .eq. imiss) then
              ddata(dyear,im) = amiss
            else
              if(inel .eq. 4 .and. idata(im) .lt. 0) then
                print *,' precip less than 0!!!!!', dstn, dyear, i
              endif    
              ddata(dyear,im) = float(idata(im))/rdiv
              summth(im) = summth(im) + ddata(dyear,im)
              nummth(im) = nummth(im) + 1
            endif  
          end do
        else if(indatyp .eq. 1) then
c         incoming data in TD3220 format (originally for HOFN)
          read(17,1500,end=95) cRecType,dstn,cElmType,dYear,
     &        iTotlVal,(iMonth(iNVal),iDay(iNVal),iValue(iNVal),
     &        iNVal=1,iTotlVal)
 1500     format (a3,a6,2x,a4,2x,i4,6x,i3,12(2i2,i6,2x))
          if(cRecType .ne. 'MLY') then
            print *, nfname, ' is NOT monthly data'
            goto 95
          endif  

          if(dyear .lt. ibegyr .or. dyear .gt. ilstyr) go to 95
          rdiv = adiv(1)
          if((cElmType .eq. 'MMXT' .and. inel .eq. 1 ) .or. 
     *      (cElmType .eq. 'MMNT' .and. inel .eq. 2)) then 
            if (iTotlVal .gt. 12) iTotlVal = 12
            do iNVal=1,iTotlVal
              im = iMonth(iNVal)
c             !!!!!! Assumed NO Estimated input !!!!!!
              ddata(dyear,im) = float(iValue(iNVal))/rdiv
              summth(im) = summth(im) + ddata(dyear,im)
              nummth(im) = nummth(im) + 1
            end do
          endif  

        else
c         reference (hcn) style format (metric to english conv)
          read(17,2000,end=100,err=110) dstn,dyear,(tdata(i),i=1,12)
 2000     format(a6,3x,i4,12f8.2)
          if(dyear .lt. ibegyr .or. dyear .gt. ilstyr) go to 95
          do im = 1,12
            if(tdata(im) .lt. -99.) then
              ddata(dyear,im) = amiss
            else  
              if(inel .ne. 4) then
                tdata(im) = (tdata(im) * 1.8) + 32.0
              else
                if(tdata(im) .lt. 0) then
                  print *,' precip less than 0!!!!!', dstn, dyear, i
                endif  
                ddata(dyear,im) = tdata(im) / 2.54
                summth(im) = summth(im) + ddata(dyear,im)
                nummth(im) = nummth(im) + 1
              endif
            endif  
          end do
        endif      
   95 end do
c     check to ensure station may be of value
  100 do im = 1,12
        if(nummth(im) .lt. minmths .and. icheck .eq. 1) then
          print *,stnind, ' skipping: ',dstn,' insufficient months:',
     *       nummth
          goto 102
        endif
      enddo
c     save station data      
      do im = 1,12  
        avgmth = summth(im) / nummth(im)
        irflag(stnind,im) = nummth(im)
        do iy = ibegyr, ilstyr
          if(ddata(iy,im) .ne. amiss)
     *      rdata(stnind, iy, im) = ddata(iy,im) - avgmth
c          print *, stnind, iy, rdata(iy, stnind)
        enddo
      end do

  102 close(17, err=105)        
      return

  105 print *, 'cannot close: ', nfname(stnind)
      go to 150
      
c     print error message and abort program
  110 print *, 'cannot open or read: ',
     *  nfname(stnind)(1:lnblnk(nfname(stnind))),' - skip'
      return

  150 print *, 'program is aborting.'
      stop

      end

c***********************************************************************
c end of subroutine getdata.                                            *
c***********************************************************************


c***********************************************************************
c                          subroutine sortco                           *
c                                                                      *
c date of last change:  24 march 1994                                  *
c                                                                      *
c modifier:  david bowman                                              *
c                                                                      *
c version     modification                                     date    *
c -------     ------------                                  ---------- *
c   1.1       sun workstation version                       03/24/1994 *
c                                                                      *
c   1.0       original                                                 *
c                                                                      *
c author:  claude williams                                             *
c                                                                      *
c modifiers:  pam hughes                                               *
c             david bowman                                             *
c                                                                      *
c usage:  sortco(double precision corrs(nstns-1), integer nstns,       *
c                integer numnbr, integer ptr2(nstns-1))                 *
c                                                                      *
c description:  this subroutine determines and sorts the 10/20 most    *
c               correlated of a hcn station's 40 nearest neighbors.    *
c                                                                      *
c notes:  changed to sort only in descending order, to sort            *
c         correlation array as well as pointer array, and to sort only *
c         for the 10/20 best correlated neighbors.                     *
c         - d. bowman:  01/21/94                                       *
c                                                                      *
c results:  sorts the first 10/20 of the correlation and pointer       *
c           arrays for the 40 nearest neighbors for each candidate     *
c           station.                                                   *
c                                                                      *
c variables:                                                           *
c                                                                      *
c   corrs      correlations between the candidate station and the      *
c   (nstns-1)   station's 40 nearest neighbors                          *
c                                                                      *
c   cortmp     holds correlation value when sorting                    *
c                                                                      *
c   i          loop counter                                            *
c                                                                      *
c   j          loop counter                                            *
c                                                                      *
c   max        index of maximum value                                  *
c                                                                      *
c   nstns     number of values to sort                                *
c                                                                      *
c   numnbr     number of hcn nearest neighbors                         *
c                                                                      *
c   ptr2       pointer array that sorts the 10/20 best correlated of a *
c   (nstns-1)   hcn station's 40 nearest neighbors                      *
c                                                                      *
c   ptrtmp     holds pointer value when sorting                        *
c                                                                      *
c***********************************************************************

      subroutine sortco(corrs,nstns,numnbr,ptr2)

      real corrs(*)
      integer ptr2(*)

      real cortmp
      integer ptrtmp

c     sort in descending order
      do 20 i = 2,nstns
        max = i
        do 10 j = i+1,numnbr
          if(corrs(j) .gt. corrs(max)) max = j
   10   continue

c       order correlation and pointer arrays
        if(max .ne. i) then
          cortmp = corrs(i)
          corrs(i) = corrs(max)
          corrs(max) = cortmp

          ptrtmp = ptr2(i)
          ptr2(i) = ptr2(max)
          ptr2(max) = ptrtmp
        end if

   20 continue

      return

      end

c***********************************************************************
c end of subroutine sortco.                                            *
c***********************************************************************

c     =======================================================================

      subroutine frstdif (nx,x,y,amiss)
c*****************************************************************
c
c     subroutine frstdif computes a first difference time series of
c     monthly mean temperatures taking into account missing monthly
c     values.  7 september, 1999; 26 november, 1999
c
c                               -on input-
c
c       nx      the number (odd) of monthly mean temperatures
c
c       x       the 2-d (ns, nx) monthly mean temperatures
c
c                               -on output-
c
c       ny      the number (even) of first differences ny = nx -1
c
c       y       the 2-d (ns, ny) first differences
c
c******************************************************************
      real x(nx), amiss
      real y(nx)
      kount = 0

c     initialize
      do im = 1, nx
        y(im) = amiss
      enddo

c
c     find first good value
c
      do im1 = 1,nx
        if (x(im1).ne.amiss) goto 5
      enddo
      goto 99    
c
c     do first difference filter accounting for missing data
c
    5 kount = 1
c      print *, 'first', im1
      do im = im1+1,nx
        if (x(im) .eq. amiss) then
          kount = 0
        else
          kount = kount + 1
          if (kount.gt.1) then
            y(im-1) = (x(im) - x(im-1))/2
c            print *,im,kount,x(im),x(im-1),y(im-1)
          endif  
        endif
      enddo  
 99   return
      end     

c     =======================================================================

      subroutine correl(nx,x,y,amiss,ccoef,ysd,xsd,kk)
c********************************************************************
c
c     subroutine correl computes the pearson correlation coefficent
c     between a candidate station and each of its neighbors.  account
c     is taken of missing data in the candidate and each neighbor so
c     that the mean, stnd dev, and cross-products are computed only
c     for those months when both candidate and neighbor have good,
c     as opposed to missing, data.  each correlation coefficient applies
c     to all months of the year.
c     9 september, 1999; 29 november, 1999
c
c                               -on input-
c
c       ns      number of stations (including candidate)
c
c       nm      number of months of good and missing data in each 
c               time series of monthly departures.
c
c       x       the (nm x ns) array of monthly departures for the candidate
c               and successive neighbor stations - each array sum is zero
c
c                               -on output-
c
c       ccoef   the (ns-1 x 1) array of candidate-neighbor pearson correlation
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
      real x(nx),y(nx),ccoef
      real xn(nx),yn(nx)
c
      eps = 0.000001

c     determine those months that both candidate and a neighbor have
c     good data.  they will vary with the particular candidate-neighbor
c     combination.  
c
      ccoef = 1.0
      kk = 0       !the counter for good candidate-neighbor data
      do im = 1,nx
        if (x(im).ne.amiss.and.y(im).ne.amiss) then
          kk = kk + 1
          xn(kk) = x(im)
          yn(kk) = y(im)
c          print *, kk, xn(kk), yn(kk)
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
          xmn = amiss
          ymn = amiss
          goto 11
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
c       calculate pearson correlation coefficient for each candidate-
c       neighbor combination
c
        sum = 0.0
        do 50 k = 1,kk
          sum = sum + (xn(k)-xmn)*(yn(k)-ymn)
 50     continue
        if(ysd .lt. eps) ysd = eps
        if(xsd .lt. eps) xsd = eps
        ccoef = sum/((kk-1)*xsd*ysd)
c        print *, kk, sum, xsd, ysd, ccoef
        goto 10

 11     ccoef = amiss

 10   continue
      return
      end


