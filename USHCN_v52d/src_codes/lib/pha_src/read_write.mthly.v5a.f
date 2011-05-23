c    ------------------------ Official Version USHCNv2 -------------------
c
c        5a   28may09   restucture for new data type subdir from state

      subroutine commline()
      
c     read the command line parameters given as pairs of "-option value"
      integer iargc, narg
      character*132 argv
      
c     This include file contains the system parameters
      INCLUDE 'inhomog.parm.mthly.incl'  
c     command line parameters
      include 'inhomog.comm.mthly.incl'
      
c     initialize command line parameters
      do it = 1, ntech
        itech(it) = 0
      enddo  
      inel = 0
      idebug = 0
      icorr = 0
      icoin = 1
      irecurse = 0
      iregion = 0
      cmetafmt = 0
      rmetafmt = 0
      ncand = 0
      immts = 0
      qscale = 0.0
      ihyear = 0
      ihtag = ''
      ctype = ''
      ntype = ''
      otype = ''
      incand = ''
      incoop = ''
      outcand = ''
      outcoop = ''
      odir = ''
      netfile = ''
      candfile = ''
      reffile = ''
      jrnlfile = ''
      mattdata = ''
      unique = ''
c     random gives number of generated random series network stations
      irandom = 0
      minslpknt = 0
      itimeres = -1
      
c     parameters defined in the confirm-filter study
      print *,' ---- Parameters defined in Confirm-Filter Study ----'
      iconfirm = 2
      print *, ' Confirmation number :', iconfirm
      nmrgyr = -2
      print *, ' Merge months depend on chgpt amp :', nmrgyr
      itech(1) = 1
      itech(2) = 0
      itech(3) = 0
c     NOTE: the SHAP output (indKW := ntech+1) is always turned ON!
      itech(indKW) = 1

      print *, ' T-test (+KW) Processing always enabled'

c     list of command line parms, in order
c     cC-d-e-F-g-H-j-lL-mM-nN-oO-pP-qQ-rR-sS-tT-u-W

      print *,' ---- Command line input ----'
c     get the number of command line arguments
      iarg = 1
      narg = iargc()
c     go thru command line arguments - keep the valid ones
c        and ignoring the undefined ones.
      do while(iarg .le. narg)
        call getarg(iarg, argv)
        iarg = iarg + 1
        if(argv .eq. '-n') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          netfile = argv
          print *, ' Cand-Ref Network Input file :', 
     *      netfile(1:lnblnk(netfile))
        else if(argv .eq. '-T') then
c          call getarg(iarg,argv)
c          print *,ntech, ':', argv(1:lnblnk(argv)), ':'
          iarg = iarg + 1
c          read(argv,fmt='(3i1)') (itech(i),i=1,ntech)
c          do i = 1,ntech
c            if(itech(i) .eq. 1) then
c              print *, ctech(i), ' processing option enabled'
c              iptech = iptech + 1
c            endif
c          enddo    
        else if(argv .eq. '-m') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          mattmeta = argv
          print *,' Matt meta Input file :',mattmeta(1:lnblnk(mattmeta))
         else if(argv .eq. '-r') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          mattdata = argv
          print *,' Matt data Input file :',mattdata(1:lnblnk(mattdata))
        else if(argv .eq. '-j') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          jrnlfile = argv
          print *, ' ----------- RESTART RUN --------------------'
          print *, ' Journal Input file :', jrnlfile
        else if(argv .eq. '-c') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i4)') ncand
          print *, ' Number of candidate stations in meta: ', ncand
        else if(argv .eq. '-t') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i1)') itimeres
          print *, ' Time resolution (0-ann;1-mthly): ', itimeres
        else if(argv .eq. '-o') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          otype = argv
          print *, ' Processed Stage for Homog Output :', 
     *      otype(1:lnblnk(otype))
        else if(argv .eq. '-p') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          ctype = argv
          print *, ' Processed Stage for Candidate Input :', 
     *      ctype(1:lnblnk(ctype))
        else if(argv .eq. '-q') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          ntype = argv
          print *, ' Processed Stage for Network Input :',
     *      ntype(1:lnblnk(ntype))
        else if(argv .eq. '-u') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          unique = argv
          print *, ' Unique descriptor for this run/rerun :', 
     *      unique(1:lnblnk(unique))
        else if(argv .eq. '-C') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          incand = argv
          ilen = lnblnk(incand)
          if(incand(ilen:ilen).ne.'/') incand(ilen+1:ilen+1)='/'
          print *, ' Base Candidate Input Directory ', 
     *      INCAND(1:lnblnk(INCAND))
          call getarg(iarg, argv)
          iarg = iarg + 1
          outcand = argv
          ilen = lnblnk(outcand)
          if(outcand(ilen:ilen).ne.'/') outcand(ilen+1:ilen+1)='/'
          print *, ' Base Candidate Output Directory: ', 
     *      OUTCAND(1:lnblnk(OUTCAND))
        else if(argv .eq. '-O') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          odir = argv
          ilen = lnblnk(odir)
          if(odir(ilen:ilen).ne.'/') odir(ilen+1:ilen+1)='/'
          print *, ' Directory for Test or HOFN output files :', 
     *      odir(1:lnblnk(odir))
        else if(argv .eq. '-N') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          incoop = argv
          ilen = lnblnk(incoop)
          if(incoop(ilen:ilen).ne.'/') incoop(ilen+1:ilen+1)='/'
          print *, ' Base Coop Input Directory: ', 
     *      INCOOP(1:lnblnk(INCOOP))
          call getarg(iarg, argv)
          iarg = iarg + 1
          outcoop = argv
          ilen = lnblnk(outcoop)
          if(outcoop(ilen:ilen).ne.'/') outcoop(ilen+1:ilen+1)='/'
          print *, ' Base Coop Output Directory: ', 
     *      outcoop(1:lnblnk(outcoop))
        else if(argv .eq. '-e') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i3)') inel
          if(inel .eq. 4) then
            print *, ' Error: Use v2 only for temperature'
            print *, '        Use v3 for precipitation'
            go to 10
          else if(inel .gt. maxelem .or. inel .lt. 1) then
            print *, ' Element parameter out of range'
            go to 10
          endif 
          icelem = celem(inel) 
          print *, ' Processing Meteorological Element:', inel, 
     *      ' : ', icelem          
        else if(argv .eq. '-d') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i1)') idebug
          print *, ' Debug Level:', idebug
        else if(argv .eq. '-H') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i4)') ihyear
          itemp = ihyear
          if(ihyear .lt. 0) then
c     -- CHANGES BY DANIEL ROTHENBERG <darothen@mit.edu> -- BEGIN
c        -- need to put parentheses around the -1 else the compiler
c        -- misinterpets what operation is wanted here
            itemp = ihyear * (-1)
c     -- CHANGES -- END
          else if(ihyear .eq. 0) then
            itemp = begyr  
          endif
          if(itemp .lt. begyr) then
            print *,' Error: abs(ihyear):',itemp,
     *        ' is less than begyr:', begyr
            goto 10
          endif  
          call getarg(iarg,argv)
          iarg = iarg + 1
          ihtag = argv
          print *, ' HOFN Input and Graph Output enabled start:', 
     *      ihyear, ' Out DTAG: ', ihtag(1:lnblnk(ihtag))
        else if(argv .eq. '-S') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i2)') minslpknt
          print *, ' Min Slope Seg:', minslpknt
        else if(argv .eq. '-s') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i3)') ndellim
          print *, ' Suspect ndellim:', ndellim
        else if(argv .eq. '-R') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i3)') irandom
          print *, ' Random Series Test neighbors:', 
     *      irandom 
        else if(argv .eq. '-F') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i6)') firstnet
          print *, ' First Network (by count) set to:', firstnet
        else if(argv .eq. '-L') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i6)') lastnet
          print *, ' Last Network (by count) set to:', lastnet
        else if(argv .eq. '-W') then
          irecurse = 1
          print *,' Recursion ENABLED with WMs series as input'
        else if(argv .eq. '-g') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i6)') iregion
          print *,' Output for US Region only:', iregion
        else if(argv .eq. '-l') then
          netloop = 1
          print *, ' Network Looping ENABLED'
c          print *, ' Reference metafile in GHCN format'
        else if(argv .eq. '-M') then
          immts = 1
          print *, ' USHCN v1 MMTS adjustment enabled'
        else if(argv .eq. '-P') then
          icoin = 0
          print *, ' Post Threshold test ENABLED'
        else if(argv .eq. '-Q') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(f5.2)') qscale
          print *, ' Inner-quartile filter scale: ', qscale
        else
          print *, ' Unknown argument :', argv, ': skipping'
        endif
      enddo
      
      if(iregion .eq. 0) then
        print *,' No Region specified - All stations will be output'
      endif  
      
      if(ncand .eq. 0) then
        print *,' NOTE: NO CANDIDATE STATIONS - ALL COOP DATA'
      endif  

      if(minslpknt .eq. 0) then
        minslpknt = minlen
        print *,' Minslpknt default (minlen):', minslpknt
      endif  

      if(icorr .eq. 0) then
        print *,' All pairs equally weighted in adj. est.'
      endif  
      
      if(icoin .eq. 1) then
        print *,' Co-incident Threshold testing ENABLED'
      endif  
      
      if(immts .eq. 0) then
        print *,' USHCN v1 MMTS not enabled'
      endif  
      
      if(netfile.ne.'' .and. incand .ne. '' .and. itimeres .ge. 0 .and.
     *   outcand .ne. '' .and. outcoop .ne. '' .and.
     *   incoop.ne. '' .and. ctype.ne.'' .and. odir.ne.''.and.
     *   otype.ne.'' .and. ntype.ne.'' .and. qscale .ne. 0. .and.
     *   ndellim .ge. 0 .and.
     *   (itech(1).ne.0.or.itech(2).ne.0.or.itech(3).ne.0)) goto 11
      
   10 print *,' Apply Inhomogeneity techniques to the Monthly Data'
      print *,'     TEMPERATURE ONLY'
      print *,'   REQUIRED PARAMETERS'
      print *,'     -T            Techniques to use (TPR0,TPR2,TPR1)'
      print *,'                   = 100 (version 52d)'
      print *,'     -l            Toggle network looping option ON'
      print *,'                   Default is CLASSIC no looping'
      print *,'     -c ncand      Num USHCN Stations in netfile'      
      print *,'     -t itimeres   Time resolution (0-ann;1-mthly)'
      print *,'     -Q qscale    Scale for inner-quartile filter'
      print *,'                   = 1.46 (version 52d)'
      print *,'     -s ndellim    Ndelete threshold for suspect data'
      print *,'     -g reg ne=1,se=2,ct=3,nc=4,so=5,np=6,sw=7,nw=8,we=9'
      print *,'    Input Files ---'
      print *,'     -n netfile    Cand-Ref Network Stations file'
      print *,'     -o otype      Output process level'
      print *,'     -p ctype      Candidate Input process level'
      print *,'     -q ntype      Network Input process level'
      print *,'     -e elem       Met elem(1=max,2=min,3=avg)'
      print *,'    Input Directories ---'
      print *,'     -C incand outcand  I/O Candidate Data Directories'
      print *,'     -N incoop outcoop  I/O Coop Data Directories'
      print *,'If neither unique nor journal entered - RESTART UNAVAIL'
      print *,'     -u unique     descriptor for this run/rerun'
      print *,'     -j journal   Input the journal file and restart'
      print *,'   MISC PARAMETERS'
      print *,'     -d           Debug level'
      print *,'     -S minslpknt Min seg length for sloped model'
      print *,'     -r mattdata  Matt data Input file'
      print *,'     -M           MMTS v1 Adjustment (Off is default)'
      print *,'     -m metafile  Matt meta Input file'
      print *,'     -H ihyear ihtag Generate output for HOFN Graphs'
      print *,'                    ihyear is the first year to output'
      print *,'                    if <= 0, used HCN Normal input'
      print *,'                    if > 0, use HOFN 3220 input'
      print *,'                    ihtag is the datetag for files'
      print *,'     -O odir      Directory for Test&Graph output files'
      print *,'     -R nneigh    Random test (with #neigh) Mmenne fmt'
      print *,'     -P           Toggle Post/Coincident Threshold test'
      print *,'     -W           Recursion ENABLED w/ WMs series input'
      print *,'     -B begstn    first station number to process ',
     *  '(default - first in file)' 
      print *,'     -E endstn    last station number to proces ',
     *  '(default - last in file)'
      print *,'     -F firstnet    first network to process ',
     *  '(default - first in file)' 
      print *,'     -L lastnet    last network to proces ',
     *  '(default - last in file)'
      print *,' NOTE: Raw data files assumed to be parsed into '
      print *,'        sudirectories by state (01,02...,50)'
      print *,''
      print *,' ------------ What is needed -------------'
      if(netfile.eq.'') print *,' Netfile missing'
      if(incand .eq. '') print *,' Input Candir missing'
      if(outcand .eq. '') print *,' Output Candir missing'
      if(itimeres .lt. 0) print *,' Timeres negative'
c      if(unique .eq. '' .and. jrnlfile .eq. '') 
c     *  print *,' Define Unique OR Journal'
      if(incoop.eq. '') print *,' Input Coopdir missing'
      if(outcoop.eq. '') print *,' Output Coopdir missing'
      if(ctype.eq.'') print *,' Ctype missing'
      if(odir.eq.'') print *,' Odir missing'
      if(otype.eq.'') print *,' Otype missing'
      if(ntype.eq.'') print *,' Ntype missing'
      if(qscale .eq. 0.) print *,' Qscale ne zero'
      if(ndellim .eq. 0.) print *,' Ndellim ne zero'
      if(itech(1).eq.0.and.itech(2).eq.0.and.itech(3).eq.0)
     *  print *,' At least one technique must be defined'
      
      stop

   11 continue

      return
      end

c     =======================================================================
      subroutine openunits(version)
      INCLUDE 'inhomog.comm.mthly.incl'
      
      character*132 version
      character*132 outfile
      
c     Open the output files to accumulate the data for each station
c       for all the iterations
      do intech = 1, ntech + 1
        if(itech(intech) .eq. 1) then
          outfile = odir(1:lnblnk(odir)) // version(1:lnblnk(version)) 
     *      // '.' // c2tech(intech)
          iounit = bunit + intech
          open(iounit,file=outfile,status='unknown')
        endif
      enddo

      return
      end
     
c     =======================================================================
      subroutine closeunits()
      INCLUDE 'inhomog.comm.mthly.incl'

c     Close all the debug output files that accumulate data for each
c       station for all the iterations
      do intech = 1, ntech + 1
        if(itech(intech) .eq. 1) then
          iounit = bunit + intech
          close(iounit)
        endif
      enddo

      return
      end

c     =======================================================================
      subroutine readnet(mmunit, nnunit,idunit,rnunit,orig,intr,ieof,
     *  ntstn)

c     WARNING: use only with inhomog.parm.incl & inhomog.comm.mthly.incl 
c         AFTER 10Feb04

c     version 3 breaks the 10 subnet barrier for MM's simulations
c                                                     26 Apr 04 cw
c
c     readnet version 2 is implemented from octpart.v3.f and newer.
c       BE AWARE: ONLY USE SKILL SCORE ROUTINES VERSION 5 AND ABOVE --- as of
c         this version the pseudo station numbers for testing now begin with 830000 
c         this removes the error of the **** in the last station, but also
c         offsets the station number by 10 for the older skill score programs!
c
c       The connection of the network-station read being one is removed
c       The target (candidate) station is the first subnet (neighbor) station
c       All of the other stations in the subnet must be a candidate in their
c         own subnet for the network to be closed (or circular) otherwise
c         the code won't work properly.                26sep03 cw
c
c     readnet reads in the candidate and its network of stations
c       nnunit is the unit number of the network definition file
c       idunit is the unit number for the data files
c       orig is the data array in stn,year,mth index order
c       ieof returns whether the end of file has been reached in the
c          network definition file (=1) or not (=0)

c     cstn is the desired station id
c     ista is the station's index

c     This include file contains the system parameters
      INCLUDE 'inhomog.parm.mthly.incl'  
      
c     command line parameters
      include 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'

      parameter(nmth = 13)
      parameter(iMaxVal = 100)

c     Raw monthly temperature data array (stn/yr/monann)
      real   orig(maxstns,begyr:endyr,13)
c     Indices of all the stations in each neighborhood
      integer indx(nstns)
c     Need to run, initially contains zeros and only when the candidate station 
c       with that index is read in, does it change to 1
      integer intr(maxstns)
c     for annual computations      
      integer num(maxstns),rnunit
      real sum(maxstns)
      
      integer indata(13), imiss /-9999/

c     station id and lat/lon/elev
      character*6 ntstn(maxstns),instn
      dimension rStnInfo(maxstns,3)
      
c     temporary 3220 input arrays
      integer iMonth(iMaxVal),iDay(iMaxVal),iValue(iMaxVal)
      character*3 cRecType
      character*4 cElmType

c      character*132 nfname, hcndir
      character*132 nfname
      character*1 inflag(13)
      character*1 ichar(21) /'0', '1', '2', '3', '4', '5', '6', '7', 
     * '8', '9', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'/

c     MM's missing value -999.999
      real ammis /-999./
      
c     variables for text output construction due to G77 inadequacies
      character*132 outstr, tmpstr

      inumyr = endyr - begyr + 1

      do istn = 1, maxstns
        ntstn(istn) = ''
        do iy = begyr, endyr
          do imth = 1, 13
            orig(istn,iy,imth) = amiss
          enddo  
        end do
      end do
      
      do istn = 1, maxstns
        do iy = begyr, endyr
          do imth = 1, 13
            tflg(istn,iy,imth) = ' '
          enddo  
        end do
      end do
      
c     initialize need-to-run for all subnetworks
      do is = 1, maxstns
        intr(is) = 0
      enddo
      
c     Set the element to the incoming parameter and input resolution
      rdiv = 10.

      ieof = 0
        
c     For Normals or HOFN type input, read in candidate & network indices
      if(ihyear .gt. 0 .or. irandom .eq. 0) then
c       read in all subnetworks from the correlation output
        itarg = 0
        print *,' nsnet: ', nsnet
        do while (1 .eq. 1)
          read(nnunit, 1000, end=80) instn
 1000     format(a6)
          read(nnunit, *, end = 80) (indx(i),i=1,nsnet)
          read(nnunit, *, end = 80) adum
          itarg = itarg + 1
          k = indx(1)
          if(k .gt. maxstns) then
            print *,' Number input stations > maxstns ',maxstns
            stop
          endif  
          ntstn(k) = instn 
          do i = 1,nsnet
            nindx(k,i) = indx(i)
          enddo  
          read(ntstn(k), fmt='(i6)') istn
        enddo

c       if no subnets have been read, end it here....
   80   if(itarg .eq. 0) goto 200
  
c       total number of stations (& subnets)
        numsubs = indx(1)
      endif  

c     HOFN input 
      if(ihyear .gt. 0) then
      
c        print *,' Begyr: ',Begyr,' Endyr: ',Endyr,' IHyear: ',ihyear
c        print *,' icelem: ', icelem
c       read 3220 datafiles
        do iStn = 1, numsubs
      
          nfname = incand(1:lnblnk(incand)) //  ntstn(istn)
     *      // '.3220'
          open(idunit,FILE=nfname,status='old',err=92)

          intr(iStn) = 1
          ntot = 0
          do while (1 .eq. 1) 
c         incoming data in TD3220 format (originally for HOFN)
            read(idunit,1500,end=95) cRecType,dstn,cElmType,iYear,
     &        iTotlVal,(iMonth(iNVal),iDay(iNVal),iValue(iNVal),
     &        iNVal=1,iTotlVal)
 1500       format (a3,i6,2x,a4,2x,i4,6x,i3,12(2i2,i6,2x))
c            print *,' cRectype, celmtype, iyear, itotlval: ', 
c     *        cRectype, celmtype, iyear, itotlval
            if(cRecType .ne. 'MLY') then
              print *, nfname, ' is NOT monthly data'
              goto 95
            endif  

            if(iYear .ge. BegYr .and. iYear .le. EndYr ) then
              aval = 0.0
              nval = 0
              if((cElmType .eq. 'MMXT' .and. icelem .eq. '_max' ) .or. 
     *          (cElmType .eq. 'MMNT' .and. icelem .eq. '_min')) then 
        
                if (iTotlVal .gt. 12) iTotlVal = 12
                do iNVal=1,iTotlVal

c                 !!!!!! Assumed NO Estimated input !!!!!!
                  orig(iStn, iYear, iMonth(iNVal)) = iValue(iNVal)/rdiv
                  aval = aval + orig(iStn, iYear, iMonth(iNVal))
                  nval = nval + 1
                  ntot = ntot + 1
                end do
c               Generate annual
                if(nval .eq. 12) then
                  orig(iStn, iyear, 13) = aval / nval
                else
                  orig(iStn, iyear, 13) = amiss  
                endif  
              endif !End condition "TMAX or TMIN"
            endif !End condition "Year"
            
      
   90     end do
          
   92     call perror(' Cannot open HOFN(3200) data: ' // nfname)

c   95     print *, ntstn(istn), ntot
   95     close(idunit)
        end do

  100   numsubs = iStn - 1
        print *, 'Number of stations read: ', numsubs

c     read in the Normals format
      else if(irandom .eq. 0) then

c       for candidate & each reference station read in its data
        DO K = 1,numsubs
          itarg = nindx(k,1)
          if(k .le. ncand) then
            if(irecurse .eq. 0) then
              nfname = incand(1:lnblnk(incand)) // 
     *         ctype(1:lnblnk(ctype)) // '/' // ntstn(itarg) // 
     *         icelem // "." // ctype(1:lnblnk(ctype))
            else
              nfname = incand(1:lnblnk(incand)) // 'WMs.' // 
     *          ctype(1:lnblnk(ctype)) // '/' // ntstn(itarg) // 
     *          icelem // '.WMs.' // ctype(1:lnblnk(ctype))
            endif
          else
            if(irecurse .eq. 0) then
              nfname = incoop(1:lnblnk(incoop)) //  
     *          ntype(1:lnblnk(ntype)) // '/' // ntstn(itarg) // 
     *          icelem // "." // ntype(1:lnblnk(ntype))
            else
              nfname = incoop(1:lnblnk(incoop)) //  'WMs.' // 
     *          ntype(1:lnblnk(ntype)) //  '/' // ntstn(itarg) // 
     *          icelem // '.WMs.' // ntype(1:lnblnk(ntype))
            endif
          endif
c         Open temporary i/o unit for data input
          open(idunit, FILE = nfname, err = 130)
          intr(k) = 1
          do while (1 .eq. 1)
              read(idunit, '(a6,i1,i4,12(i6,a1))', end = 140)instn, 
     *          isrc, iyear,(indata(imth),inflag(imth),imth=1,12)
              if(iyear .ge. begyr .and. iyear .le. endyr) then
                aval = 0.0
                nval = 0
                do imth = 1, 12
c                 convert external format and missing to internal values
                  if(indata(imth) .eq. IMISS) then
                    orig(itarg, iyear, imth) = AMISS
c                  else if(inflag(imth).eq."S".or.inflag(imth).eq."A")then
                  else if(inflag(imth).eq."f")then
c                   filter out the NEW Integrated Monthly Dataset (after QC)
                    orig(itarg, iyear, imth) = AMISS
                    inflag(imth) = 'M'
                  else  
                    orig(itarg, iyear, imth) = indata(imth)/rdiv
                    aval = aval + orig(itarg, iyear, imth)
                    nval = nval + 1
                  endif  
                  tflg(itarg,iyear,imth) = inflag(imth)
                end do
                if(nval .eq. 12) then
                  orig(itarg, iyear, 13) = aval / nval
                else
                  orig(itarg, iyear, 13) = amiss  
                endif  
              end if
          end do
  130     call perror(' Cannot open Normals data: ' // nfname)
c         make sure to close this unit!!!
  140     close(idunit)
  150   end do  
  
      else
        numsubs = maxstns
              
c       read Mmenne random series data
c       warning: assume 100 records per candidate-network
c       Annual series are positive
  160   if(irandom .gt. 0) then
          do k = 1, numyr
            read(rnunit,*,end=180,err=200)
     *        iy,(orig(ns,iy+begyr,13),ns=1,numsubs)
c           convert MM's missing to internal missing
            do ns = 1, numsubs
              if(orig(ns,iy+begyr,13) .le. ammis)
     *          orig(ns,iy+begyr,13)=amiss
            enddo  
          enddo
        else
c       monthly series are negative
          do k = 1, numyr
            do ns = 1, numsubs
              num(ns) = 0
              sum(ns) = 0.0
            enddo  
            do m = 1,12
              read(rnunit,*,end=180,err=200)
     *          iy,(orig(ns,k+begyr,m),ns=1,numsubs)
              if(((k-1)*12 + m) .ne. iy) then
                print *,'Monthly input and YYYY/MM out of sync'
                stop
              endif  
c             convert MM's missing to internal missing
              do ns = 1, numsubs
                if(orig(ns,k+begyr,m).le.ammis)then
                  orig(ns,k+begyr,m)=amiss
                else
                  num(ns) = num(ns) + 1
                  sum(ns) = sum(ns) + orig(ns,k+begyr,m)
                endif    
              enddo  
            enddo  
            do ns = 1, numsubs
              if(num(ns) .eq. 12) then
                orig(ns,k+begyr,13) = sum(ns)/ num(ns)
              else
                orig(ns,k+begyr,13) = amiss
              endif    
            enddo  
          enddo
        endif

c       read the stations and their nieghbors (1st diff corr sort)
  180   do istn = 1, numsubs
          intr(istn) = 1
          read(nnunit, 1000, end=190) ntstn(istn)
          read(nnunit, *, end = 190) (indx(i),i=1,nsnet)
          read(nnunit, *,end=190) adum
          k = indx(1)
          do i = 1,nsnet
            nindx(k,i) = indx(i)
          enddo  
        enddo  

c       test station number against begin and end stations
c        read(ntstn(1), fmt='(i6)') istn
c        if(istn .lt. begstn) goto 160
c        if(istn .gt. endstn) goto 190
        
      endif  

c      print *,'Last station (subnet) index:',numsubs
c      print *,' readnet orig(49,1950,1):',orig(49,1950,1)
      return
      
  190 print *,' Candidate network past end station'
      ieof = 1
      return

  200 print *,' End of candidate network records'
      ieof = 1
      return
    
      end

c     =======================================================================
      subroutine writsta(itarg,cstn,dval,aval,cval,aflag,otag,
     *  iounit)

c     common output routine for the Inhomogeneity techniques
c     cstn   - input - candidate (base) station
c     dval   - input - adjusted series data
c     tflg  - input - adjusted series flags
c     aval   - input - adjustments used for series
c     cval   - input - error of adj for series
c     aflag  - input - flags for which segment series is adjusted
c     otag  - input - two letter identifier of Inhomo Tech (LV,EP,MD)
c     iounit - input - output device number

c     pulled from the stnhist subroutine in ucpmonthly.v5
c     fix up to write out results....
c     current data array for all stations
c      real   montemp(maxstns,begyr:endyr,13)
c     original input flags for all stations
c      character*1 monflag(maxstns,begyr:endyr,13)      
c
c      character*6 ntstn(maxstns),shstn(nstns)
c      character*1 adjflag(begyr:endyr,13)
c      real outtemp(begyr:endyr,13),adjtemp(begyr:endyr,13)
c 
c            do iy = begyr, endyr
c              do im = 1, 13
c                adjflag(iy,im) = ' '
c                adjtemp(iy,im) = amiss
c                if(iy .ge. istpyr .and. 
c     *            montemp(itarg, iy, 13) .ne. amiss) then
c                  outtemp(iy,im) = montemp(itarg, iy, im)
c                  write(iounit,*) ' Cand: ', itarg, ntstn(itarg), iy, im,
c     *              outtemp(iy,im)
c                else
c                 outtemp(iy,im) = amiss  
c                endif  
c              enddo
c            enddo
c              
c          call writsta(shstn(1), outtemp, monflag, adjtemp, adjflag, 
c     *      'KW', idunit)

      include 'inhomog.comm.mthly.incl'
      include 'inhomog.parm.mthly.incl'
      include 'inhomog.restart.mthly.incl'
      parameter (NMTH = 13)
      character*6 cstn
      character*7 ostr
      character*2 otag
      character*1 aflag(begyr:endyr,NMTH)
      real dval(begyr:endyr,NMTH), aval(begyr:endyr,NMTH),
     *  cval(begyr:endyr,NMTH)
      real rmult /10./
      character*132 dfile, dpath, strpath
      integer idval(NMTH), iprnt(begyr:endyr)

      integer irstate(12,9) /
     *  01, 08, 09, 31, 38, 44,  0,  0,  0,  0,  0,  0,
     *  06, 07, 17, 18, 19, 27, 28, 30, 36, 37, 43,  0,
     *  11, 12, 15, 23, 33, 40, 46,  0,  0,  0,  0,  0,
     *  13, 20, 21, 47,  0,  0,  0,  0,  0,  0,  0,  0,
     *  03, 14, 16, 22, 34, 41,  0,  0,  0,  0,  0,  0,
     *  24, 25, 32, 39, 48,  0,  0,  0,  0,  0,  0,  0,
     *  02, 05, 29, 42,  0,  0,  0,  0,  0,  0,  0,  0,
     *  10, 35, 45,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     *  04, 26,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
     
c     if there is no region specified - write all stations
      if(iregion .eq. 0) goto 10

c     region given, make sure this station in the state list
      read(cstn,'(i6)') istn

      istate = istn/10000
      do ist = 1, 12
        if(istate .eq. irstate(ist, iregion)) goto 10
      enddo
c     not in this region's state list
      goto 200
      
c     data file first
   10 if(itarg .le. ncand) then
        dpath = outcand(1:lnblnk(outcand))
        strpath = 'CandOutDir: ' 
      else
        dpath = outcoop(1:lnblnk(outcoop))
        strpath = 'CoopOutDir: '
      endif
      ostr = otag // 's.' // otype
      dfile = dpath(1:lnblnk(dpath))// ostr(1:lnblnk(ostr)) //
     *    '/' // cstn // icelem // '.' // ostr(1:lnblnk(ostr))
      strpath = strpath(1:lnblnk(strpath)) // ostr(1:lnblnk(ostr)) //
     *    '/' // cstn // icelem // '.' // ostr(1:lnblnk(ostr))

c     open output file
      print *,' Writing: ', strpath(1:lnblnk(dfile))
      open(iounit, FILE=dfile, err= 300)
      
c     HOLD !!!!!!!!!!!!!!!!!!!!!!
c     write the data in HCN format
c      do iy = 1, endyr - begyr + 1
c        write(iounit,20,ERR=300) istn, iy + begyr - 1, inel,
c     *    (dval(im,iy),tflg(itarg,iy,im),im = 1,NMTH)
c   20   FORMAT(I6.6,1X,I4,1X,I1,'F',13(F6.2,1x,a1,2x))           
c      end do

c     write the data in Normals format
      do iy = begyr, endyr
        ndat = 0
        iprnt(iy) = 0
        do im = 1, NMTH
          if(dval(iy,im) .lt. AMISS + 1.0 ) then
            idval(im) = -9999
c           if a value has been made missing - print missing anyway
            if(tflg(itarg,iy,im) .ne. ' ') ndat = ndat + 1
          else
            idval(im) = nint(dval(iy,im)*rmult) 
            ndat = ndat + 1
          endif
        enddo   
        if(ndat .gt. 0) then
          write(iounit,1000,ERR=300) cstn, inel, iy,
     *      (idval(im),tflg(itarg,iy,im), im=1,NMTH)
 1000     FORMAT(a6,i1,I4,13(i6,a1))
          iprnt(iy) = 1
        endif  
      enddo
      close(iounit)
      
      goto 100

c     adjustment file next
      ostr = otag // 'a.' // otype
      dfile = dpath(1:lnblnk(dpath))// ostr(1:lnblnk(ostr)) //
     *    '/' // cstn // icelem // '.' // ostr(1:lnblnk(ostr))
c     open output file
      open(iounit, FILE=dfile, err= 300)
      
c     HOLD !!!!!!!!!!!!!!!!!!!!!!
c     write the data in HCN format
c      do iy = 1, endyr - begyr + 1
c        write(iounit,20,ERR=300) istn, iy + begyr - 1, inel,
c     *    (dval(im,iy),tflg(itarg,iy,im),im = 1,NMTH)
c   20   FORMAT(I6.6,1X,I4,1X,I1,'F',13(F6.2,1x,a1,2x))           
c      end do

c     write the data in Normals format
      do iy = begyr, endyr
        if(iprnt(iy) .eq. 1) then
          do im = 1, NMTH
            if(aval(iy,im) .lt. AMISS + 1.0 ) then
              idval(im) = -9999
            else
              idval(im) = nint(aval(iy,im)*rmult) 
            endif
          enddo   
          write(iounit,1000,ERR=300) cstn, inel, iy,
     *      (idval(im),aflag(iy,im), im=1,NMTH)
        endif
      enddo
      close(iounit)

c     err of adj file next
      ostr = otag // 'c.' // otype
      dfile = dpath(1:lnblnk(dpath))// ostr(1:lnblnk(ostr)) //
     *    '/' // cstn // icelem // '.' // ostr(1:lnblnk(ostr))
c     open output file
      open(iounit, FILE=dfile, err= 300)
      
c     write the data in Normals format
      do iy = begyr, endyr
        if(iprnt(iy) .eq. 1) then
          do im = 1, NMTH
            if(cval(iy,im) .lt. AMISS + 1.0 ) then
              idval(im) = -9999
            else
              idval(im) = nint(cval(iy,im)*rmult) 
            endif
          enddo   
          write(iounit,1000,ERR=300) cstn, inel, iy,
     *      (idval(im),aflag(iy,im), im=1,NMTH)
        endif
      enddo
      close(iounit)
  100 return
      
  200 print *,' Perimeter station not output: ', cstn
      return
      
  300 call perror(' ERROR: writing output: ' // dfile)
      stop
      
      end
  
c     =======================================================================
      subroutine write_restart(icpos)
      INCLUDE 'inhomog.parm.mthly.incl'
      include 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'

c     Current restart read-write compatible with resver = 'RSTRTv7 '
c                                                        15may06

c     icpos - indicates the position in the method/loop calling write_restart
c        1 = Beginning of the method/loop
c        2 = Just before the confirmfilt call
c        3 = Just before the Final stnhist call
      integer icpos
c     on a restart is equivalent to rentry
      
c     take out inhomog.comm.mthly.incl as soon as possible
c      integer nloop
      
      character*132 rsfile
      
c      character*64 unique
     
      write(rsfile,1000) unique(1:lnblnk(unique)), network, method, 
     *  nloop, icpos
 1000 format(a,'.',i4.4,'_',i1,'_',i2.2,'_',i1)
      open(40,file=rsfile,access='direct',recl=1,status='unknown')
c     resver is 8 bytes long
      write(40, rec=1) resver
c     these are the variables that define the array sizes (checked on read)
      write(40, rec=9) maxstns, nstns, begyr, endyr, ntech, ninh
c     these are the settings for the process variables and the arrays
      write(40,rec=33) network, method, nloop, icpos, nht, lht, ntr, 
     *    itimeres, nhits, sahist, temp, schgpt, nchgpt, ntest, nfound, 
     *    nspan, ndelete, nindx, tflg
      close(40)
      print *, ' JOURNAL WRITTEN: ', rsfile(1:lnblnk(rsfile))
      return
      end

c     =======================================================================
      subroutine read_restart(rsnet)
c      subroutine read_restart(jrnlfile, rsnet, nloop, unique)
      INCLUDE 'inhomog.parm.mthly.incl'
      INCLUDE 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'

c     these are read from the restart journal and checked against the
c       values set in the current code
      
c     need network (output) separated from rsnet (input)
      integer rsnet
      character*8 inver
      integer inmax, instns, inbeg, inend, intech, ininh

c     take out inhomog.comm.mthly.incl as soon as possible
c      integer nloop
      
      print *,' Unique identifier from Journal: ',
     *  jrnlfile(1:lnblnk(jrnlfile)-12) 
     
      open(40,file=jrnlfile,access='direct',recl=1,status='unknown')
c     Journal version is 8 bytes long
      read(40,rec=1) inver
      if(inver .ne. resver) then
        print *,' ------------- Restart Error 1 ---------- '
        print *,' Journal vs. Current Version Mismatch'
        print *,' Journal: ',inver
        print *,' Current: ', resver
        stop
      endif  
c     these are the variables that define the array sizes (checked on read)
      read(40, rec=9) inmax, instns, inbeg, inend, intech, ininh
      if(inmax .ne. maxstns .or. instns .ne. nstns .or. 
     *  inbeg .ne. begyr .or. inend .ne. endyr .or. 
     *  intech .ne. ntech .or. ininh .ne. ninh) then 
        print *,' ------------- Restart Error 2 ---------- '
        print *,' Journal vs. Current Parameter Mismatch'
        print *,'            maxstn, nstns, begyr, endyr'
        print *,' Journal: ', inmax, instns, inbeg, inend, intech, ininh
        print *,' Current: ', maxstns, nstns, begyr, endyr, ntech, ninh
        stop
      endif  
c     these are the settings for the process variables and the arrays
      read(40,rec=33,err=100) network, method, nloop, rentry, nht, lht, 
     *     ntr, itimeres, nhits, sahist, temp, schgpt, nchgpt, ntest,
     *     nfound, nspan, ndelete, nindx, tflg
      close(40)
      print *, ' JOURNAL READ: ', jrnlfile(1:lnblnk(jrnlfile))
      rsnet = network
      print *,' Journal version: ',inver
      print *,'   maxstn, nstns, begyr, endyr, ntech, ninh',
     *  inmax, instns, inbeg, inend, ntech, ninh
      print *,' >>> Restart Network: ', rsnet
      print *,' >>> Method: ', method
      print *,' >>> Nloop: ', nloop
      print *,' >>> Rentry: ', rentry
      print *,' >>> Timeres: ', itimeres
      
c      print *,' nhits>0 read: ',maxstns, nmo
c      do istn = 1, maxstns
c        do imo = 1, nmo
c          if(nhits(istn, imo) .ne. 0) print *,istn, imo, nhits(istn,imo)
c        enddo
c      enddo    
      
c      do k = 1, maxstns
c        do istn = 2, nstns                  
c          do imo = 1, nmo
c            do it = 1, intech
c              ichg = iachar(nspan(k,istn,imo,it))
c              if(ichg .gt. 0 .or. nfound(k,istn,imo,it) .ne. '0') 
c     *          print *,'NSPAN:',k,istn,imo,it,ichg,
c     *            ' ',nfound(k,istn,imo,it)
c            enddo
c          enddo    
c        enddo
c      enddo  
      
   30 return
   
  100 print *,' Error in reading process variables and arrays'
      stop
      end
          
        
