      program ushcn_fill_2004

c   version  date    description
c      4p   12may09  rearranged directory i/o structure for proceedure type
c               instead of state directories
c
c      4a   03mar09  added 3-flag option (-3) for data i/o
c
c      -    03mar08  added 'f' for failed QA check from GHCN monthly
c
c      4    28sep06  removed 's' flag for editting
c                    (see filnet_subs_2004.v3.f)
c
c      3    21nov05  added confidence interval reading & use
c            6jun06  add 's' and 'd' flags for removal
c
c     from fix.ms10 in the 1971-2000 normals
c     fix.ms10 is the same as the last normals_fill_2000.f
c       use with the full filnet_subs_2000.f

c   normals_fill_2000.v2.f reads in the filled-in hcn series,
c      one station at a time 
c      brings in the data with the hcn filnet
c      fills in all missing data between ibyear and ieyear

c     warning: changes to parameter list must be made to all 
c         parameter lists in subroutines as well
      include 'posthomog.parm.incl' 
      include 'posthomog.fill.incl' 

      real data(nmth,nyr), dconf(nmth,nyr), norm0(nmth), norms(nmth)
      real std0(nmth), stds(nmth), adiff(nmth)
      integer ngap(nmth), idbg /0/
      character*6 cstn
      character*132 argv, basedir/''/, netfile/''/, candfile/''/,
     *  reffile/''/
      character*4 celem(4)/'_max','_min','_avg','_pcp'/
      character*1 dflag(nmth, nyr), cflag(nmth, nyr)
      integer iargc, narg
      integer cmetafmt/0/, rmetafmt/0/

c     assume candidate list is one line format
      cmetafmt = 0
      inp3flag = 0
      netdir = ''
      candir = ''
      ctype = ''
      itype = ''
      ntype = ''
      otype = ''
      ikeep = 0

      print *,' fillin ushcn & coop stations - ushcn_fill_2004'
      
c     units 10, 11, 12, & 13 are reserved for main program! 
c      subroutines must not use these units except when passed as arguments!
c     i/o unit for candidate meta data (drives processing)
      icunit = 10
c     i/o unit for candidate - reference network definition file
      nnunit = 11
c     i/o unit for data input (used many times - open/close when using)
      idunit = 12
c     i/o unit for filled candidate data in hcn format
      iounit = 13

c     begin command line code      
c     get the number of command line arguments
      iarg = 1
      narg = iargc()
c     go thru command line arguments - keep the valid ones
c        and ignoring the undefined ones.
      do while(iarg .le. narg)
        call getarg(iarg, argv)
        iarg = iarg + 1
        if(argv .eq. '-c') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          candfile = argv
          print *, ' candidate stations input meta file :', candfile
        else if(argv .eq. '-n') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          netfile = argv
          print *, ' cand-ref network input file :', netfile
        else if(argv .eq. '-p') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          ctype = argv
          print *, ' processed stage for candidate input :', ctype
        else if(argv .eq. '-i') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          itype = argv
          print *, ' processed stage for cand conf int input :', itype
        else if(argv .eq. '-q') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          ntype = argv
          print *, ' processed stage for network input :', ntype
        else if(argv .eq. '-o') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          otype = argv
          print *, ' processed stage for fill data ouput :', otype
        else if(argv .eq. '-j') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          jtype = argv
          print *, ' processed stage for fill conf int ouput :', jtype
        else if(argv .eq. '-C') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          candir = argv
          ilen = lnblnk(candir)
          if(candir(ilen:ilen).ne.'/') candir(ilen+1:ilen+1)='/'
          print *, ' base directory for candidate data files :', candir
        else if(argv .eq. '-N') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          netdir = argv
          ilen = lnblnk(netdir)
          if(netdir(ilen:ilen).ne.'/') netdir(ilen+1:ilen+1)='/'
          print *, ' base directory for network data files :', netdir
        else if(argv .eq. '-e') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i3)') inel
          if(inel .gt. 4 .or. inel .lt. 1) then
            print *, ' element parameter out of range'
            go to 10
          endif 
          icelem = celem(inel) 
          print *, ' processing meteorological element:', inel, 
     *      ' : ', icelem          
        else if(argv .eq. '-d') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i1)') idbg
          print *, ' debug level:', idbg
        else if(argv .eq. '-k') then
          ikeep = 1
          print *, ' keep suspect data enabled'
        else if(argv .eq. '-3') then
          inp3flag = 1
          print *, ' Three Flag Input ENABLED'
        else if(argv .eq. '-u') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i4)') nc
          print *, ' candidate-neighbor merge - nc:',nc
        else
          print *, ' unknown argument :', argv, ': skipping'
        endif
      enddo
      
      if(candfile.ne.'' .and. netfile.ne.'' .and. candir .ne. '' .and.
     *   netdir.ne. '' .and. ctype .ne.'' .and. itype .ne.'' .and.
     *   jtype .ne.'' .and. ntype .ne. '' .and. otype.ne.'') goto 11
      
   10 print *,' fillin missing 1971-2000 normals routine'
      print *,' usage: normals_fill_2000 <-t> -e elem -c candfile ',
     *  '-n netfile -r reffile'
      print *,'   required parameters'
      print *,'    input files ---'
      print *,'     -u nc         number of candidate stations'
      print *,'     -c candfile   candidate stations file'      
      print *,'     -n netfile    cand-ref network stations file'
      print *,'     -p process    candidate input process level'
      print *,'     -i process    candidate confidence int proc level'
      print *,'     -q process    network input process level'
      print *,' --- new!!! '
      print *,'     -o process    fillin output process level'
      print *,'     -j process    fillin conf int output proc level'
      print *,'     -e elem       met elem(1=max,2=min,3=avg,4=pcp)'
      print *,'     -k            keep suspect data (default=fill)'
      print *,'    input directories ---'
      print *,'     -C candir    base candidate data directory'
      print *,'     -N netdir    base network data directory'
      print *,' note: raw data files assumed to be parsed into '
      print *,'        sudirectories by state (01,02...,50)'
      stop

   11 continue

      if(ikeep .eq. 0) print *,' suspect data removed & filled'

      if(inp3flag .eq. 0) print *,' 1 flag Normals input assumed'

c     assumptions:
c       changes from the normals to ushcn fill routine
c       candidate and reference metafiles in ushcn format
c
c       nc separates the ushcn from the coop stations in:
c       candidate == reference station file (all stations)
c       cand-ref network is all stations
c
      cmetafmt = 0
      rmetafmt = 0
      reffile = candfile

c     initialize reference station arrays for fillin
      call readrsta(reffile, idunit, rmetafmt, nc)

      open(nnunit, name=netfile)
      iskip = 0
      nwrite = 0
       
c     for each station
      do istn = 1, nstn
c       read station info into data
        call readsta(istn, cstn, data, dconf, dflag, cflag, isgood, 
     *    idunit, cmetafmt,idbg, nc)
        if(cstn .eq. "" .and. isgood .eq. 1) then
          print *," end of station list"
          go to 200
        end if  
        if(isgood .eq. 0) then
          print *," skipping "
          iskip = iskip + 1
          go to 100
        end if  
        
c       initialize the filnet data and network for current station
        call filinit(cstn, ista, nnunit, idunit, igood, idbg, inel)        
        
        if(igood .ne. 1) go to 100
            
        if(idbg .gt. 0) print *, "stn: ", cstn
c       ----- uncomment to calculate full normals
c        call calcnorm(data, norm0, 1, std0)
c        write (*,'(a6,"  norm0:", 12f8.2)'), 
c     *    cstn, norm0
c        write (*,'(a6,"   std0:", 12f8.2)'), 
c     *    cstn, std0
c
     
c       fill the gaps with filnet
        call filentry(data, dconf, dflag, cflag, ngap, 0, 0)
            
c        call calcnorm(data, norms, 0, stds)
c        write (*,'(a6,"  norms:", 12f8.2)'), 
c     *    cstn, norms

c       calculate stats for given station
c        do im = 1, nmth
c          adiff(im) = norm0(im) - norms(im)
c        end do
c        write (*,'(a6,"  diffs:", 12f8.2)'), 
c     *    cstn, adiff
c        write (*,'(a6,"  ngaps:", 12i8)'), 
c     *    cstn, ngap
     
c       10 year limit removed for nws48 - 31 jul 04 cw
c       if station has any month less than 10 values, do not print
c        do im = 1, nmth
c          if(ngap(im) .gt. 20) then
c            print *, cstn, ' less than 10 years - no output'
c            go to 100
c          endif  
c        enddo
c       ------ uncomment to calculate full normals

        nwrite = nwrite + 1  
        call writsta(istn,cstn,data,dconf,dflag,cflag,iounit,nc)

  100 end do ! end of stations loop
      
  200 print *,' number of serially complete stations ',
     *  nwrite 
      
c     close down input files
      close(icunit)
      close(nnunit)
      end      
  
c ***********************************************************
      subroutine readsta(istn,cstn,data,dconf,dflag,cflag,isgood,
     *  idunit,cmetafmt,idbg,nc )

      include 'posthomog.parm.incl'  
      include 'posthomog.fill.incl'

      character*6 cstn, instn
      character*132 dfile
      real ddata(nmth)
      integer idata(nmth)
      real data(nmth, nyr), dconf(nmth, nyr), aint
      character*1 dflag(nmth, nyr), cflag(nmth, nyr)
      character*3 iflag(nmth)
      character*1 aflag
      integer cmetafmt, itstn
      real adiv(2) /10.,100./
      integer imiss /-9999/
      
      imt2en = 0

      if(inel .ne. 4) then
        rdiv = adiv(1)
      else
        rdiv = adiv(2)
      endif  
            
c     initialize data array
      do im = 1, nmth
        do iy = 1, nyr
          data(im, iy) = amiss
          dconf(im,iy) = amiss
          dflag(im, iy) = ' '
          cflag(im,iy) = 'x'
        end do
      end do    
      
c     assume all data present
      isgood = 1

      write(cstn, '(i6.6)') rid(istn)
      if(idbg .gt. 0) 
     *  print *,' istn, nfname: ', istn,nfname(istn),inp3flag
     
c     open temporary i/o unit for data input
      open(idunit, name = nfname(istn))
      do while (1 .eq. 1)
c       get data values and day flags
        if(imt2en .eq. 0) then
          if(inp3flag .eq. 0) then
c           candidate station data in normals (1-flag) format 
            read(idunit,1000,end=100,err=100) instn,iyear,
     *        (idata(i),iflag(i),i=1,12)
 1000       format(a6,1x,i4,12(i6,a1))
            if(iyear .lt. ibyear .or. iyear .gt. ieyear) go to 80
            do i = 1, 12
c             decode observations and check QC flags
              if(idata(i) .eq. imiss) then
                data(i, iyear - ibyear + 1) = amiss
              else if(ikeep .eq. 0 .and.
     *           (iflag(i)(1:1) .eq. 'd')) then
                data(i, iyear - ibyear + 1) = amiss
              else if(iflag(i)(1:1) .eq. 'g' .or. 
     *          iflag(i)(1:1) .eq. 'S' .or. iflag(i)(1:1) .eq. 'X' .or.
     *          iflag(i)(1:1) .eq. 'f') then
                data(i, iyear - ibyear + 1) = amiss
              else
                data(i, iyear - ibyear + 1) = float(idata(i))/rdiv
              endif  
              dflag(i, iyear - ibyear + 1) = iflag(i)
            end do
          else
c           normals 3-flag format
            read(idunit,1500,end=100,err=100) instn,iyear,
     *        (idata(i),iflag(i),i=1,12)
 1500       format(a6,1x,i4,12(i6,a3))
            if(iyear .lt. ibyear .or. iyear .gt. ieyear) go to 80
            do i = 1, 12
              if(idata(i) .eq. imiss) then
                data(i, iyear - ibyear + 1) = amiss
              else if(iflag(i)(2:2) .ne. ' ') then
c               any flag means month value failed a QC test
                data(i, iyear - ibyear + 1) = amiss
                print *,' Editted: ',instn,iyear,i,iflag(i)
              else  
                data(i, iyear - ibyear + 1) = float(idata(i))/rdiv
              endif
              dflag(i, iyear - ibyear + 1) = iflag(i)
            end do    

          endif

        else
c         reference (hcn) style format (metric to english conv)
          read(idunit,2000,end=100,err=100)instn,iyear,(ddata(i),i=1,12)
 2000     format(a6,3x,i4,12f8.2)
          if(iyear .lt. ibyear .or. iyear .gt. ieyear) go to 80
          do i = 1,12
            if(ddata(i) .lt. (amiss + 1.0)) then
              data(i, iyear - ibyear + 1) = amiss
              if(idbg .eq. 2) print *, 'input missing:',instn,iyear,i
            else  
              if(inel .ne. 4) then
                data(i, iyear - ibyear + 1) = (ddata(i) * 1.8) + 32.0
              else
                data(i, iyear - ibyear + 1) = ddata(i) / 25.4
              endif
            endif  
          end do
        endif      

   80 end do
c     make sure to close this unit!!!
  100 close(idunit)
  
c     skip reading APH confidence intervals
      goto 195
      
c     open temporary i/o unit for data input
      open(idunit, name = ifname(istn))
      do while (1 .eq. 1)
c        originally ----
c        read(idunit, '(a6,1x,i1,1x,i4,12f8.2)', end = 100)
c     *    instn, isrc, iyear, indata
c       get data values and day flags
        if(imt2en .eq. 0) then
c         candidate station data in normals format          
          read(idunit,1000,end=190,err=190) instn,iyear,
     *      (idata(i),iflag(i),i=1,12)
          if(iyear .lt. ibyear .or. iyear .gt. ieyear) go to 180
          if(inel .ne. 4) then
            rdiv = adiv(1)
          else
            rdiv = adiv(2)
          endif    
          do i = 1, 12
            if(idata(i) .eq. imiss) then
              dconf(i, iyear - ibyear + 1) = amiss
            else  
              dconf(i, iyear - ibyear + 1) = float(idata(i))/rdiv
              cflag(i, iyear - ibyear + 1) = iflag(i)
            endif  
          end do
        else
c         reference (hcn) style format (metric to english conv)
          read(idunit,2000,end=190,err=190)instn,iyear,(ddata(i),i=1,12)
          if(iyear .lt. ibyear .or. iyear .gt. ieyear) go to 180
          do i = 1,12
            if(ddata(i) .lt. (amiss + 1.0)) then
              dconf(i, iyear - ibyear + 1) = amiss
              if(idbg .eq. 2) print *, 'input missing:',instn,iyear,i
            else  
              if(inel .ne. 4) then
                dconf(i, iyear - ibyear + 1) = (ddata(i) * 1.8) + 32.0
              else
                dconf(i, iyear - ibyear + 1) = ddata(i) / 25.4
              endif
            endif  
          end do
        endif      

  180 end do
c     make sure to close this unit!!!
  190 close(idunit)
      
c     initialize all of the conf intervals and flags starting at the 
c       present and going backwards
      aflag = 'a'
      aint = 0.0
      do iy = nyr, 1, -1
        do im = 12, 1, -1
          if(dconf(im,iy) .ne. amiss) then
            aflag = cflag(im,iy)
            aint = dconf(im,iy)
          else
            cflag(im,iy) = aflag
            dconf(im,iy) = aint  
          endif
        enddo
      enddo      

c     successful return
  195 return

  200 cstn = ""
      return
      
  300 print *,' error: reading station list'
      stop
      
  400 print *,' error: reading data file'
      stop        
      end
            
c ***********************************************************
      subroutine writsta(istn,cstn,data,dconf,dflag,cflag,iounit,nc)

      include 'posthomog.parm.incl'
      include 'posthomog.fill.incl'

      character*6 cstn
      character*1 dflag(nmth, nyr), cflag(nmth,nyr), annflag
      integer idval(nmth)
      real data(nmth, nyr), dconf(nmth,nyr)
      real rmult(4) /10.,10.,10.,100./
      character*132 dfile, cfile, outdir

c     assume output file has been opened somewhere else
      if(istn .le. nc) then
        outdir = candir
      else
        outdir = netdir
      endif    
      dfile = outdir(1:lnblnk(outdir)) // otype(1:lnblnk(otype)) //
     *  '/' // cstn // icelem(1:lnblnk(icelem)) // '.' // 
     *  otype(1:lnblnk(otype))
      cfile = outdir(1:lnblnk(outdir)) // jtype(1:lnblnk(jtype)) //
     *  '/' // cstn // icelem(1:lnblnk(icelem)) // '.' //
     *  jtype(1:lnblnk(jtype))

c     open data output file
      open(iounit, name=dfile)
      
c     write the data in normals format
      do iy = 1, ieyear - ibyear + 1
        annval = 0.0
        annnum = 0.0
        annflag = ' '
        ndat = 0
        do im = 1, nmth
          if(data(im,iy) .lt. amiss + 1.0 ) then
            idval(im) = -9999
          else
            ndat = ndat + 1
            idval(im) = nint(data(im,iy)*rmult(inel)) 
            annval = annval + data(im,iy)
            annnum = annnum + 1
            if(dflag(im,iy) .ne. ' ') annflag = 'i'
          endif
        enddo  
        if(ndat .eq. 12) then
          annual = nint(annval*rmult(inel)/annnum)
        else
          annual = -9999
        endif   
          
        if(ndat .gt. 0) then
          if(inp3flag .eq. 0) then
c           Normals 1-flag output
            write(iounit,1000,err=300) cstn, inel, iy + ibyear - 1,
     *        (idval(im),dflag(im,iy), im=1,nmth)
 1000       format(a6,i1,i4,13(i6,a1))
          else
c           GHCN 3-flag output
            write(iounit,2000,ERR=300) cstn, inel, iy + ibyear - 1,
     *        (idval(im),dflag(im,iy), im=1,nmth)
 2000       FORMAT(a6,i1,I4,13(i6,a3))
          endif
        endif  
      enddo

      close(iounit)
      
c     skip confidence interval print out
      goto 200

c     open conf int output file
      open(iounit, name=cfile)
      
c     write the data in normals format
      do iy = 1, ieyear - ibyear + 1
        ndat = 0
        do im = 1, nmth
          if(data(im,iy) .lt. amiss + 1.0 ) then
            idval(im) = -9999
          else
            ndat = ndat + 1
            idval(im) = nint(dconf(im,iy)*rmult(inel)) 
          endif
        enddo  
          
        if(ndat .gt. 0) then
          write(iounit,1000,err=300) cstn, inel, iy + ibyear - 1,
     *      (idval(im),cflag(im,iy), im=1,nmth)
        endif  
      enddo

  200 close(iounit)
      return
      
  300 print *,' error: writing output'
      stop
      
      end
            
c*************************************************************
      subroutine calcnorm(data, norms, stdopt, stds)

      include 'posthomog.parm.incl'  

      real data(nmth, nyr), norms(nmth), stds(nmth)
      real sum, num
      integer stdopt, inyr1
      
c     calculate the index of the beginning and end of the normals year
      inyr1 = 1971 - ibyear + 1
      inyr2 = 2000 - ibyear + 1
      
      do im = 1, nmth
        sum = 0.0
        sum2 = 0.0
        num = 0.0
        do iy = inyr1, inyr2
          if(data(im, iy) .gt. (amiss + 1.0)) then
            sum = sum + data(im, iy)
            if(stdopt .eq. 1)
     *        sum2 = sum2 + data(im, iy) * data(im, iy)
            num = num + 1
          endif  
        end do
        if(stdopt .eq. 1) then
          if(num .gt. 1) then
            stds(im) = sqrt((sum2 - (sum * sum / num)) / (num - 1))
          else
            stds(im) = amiss
          endif  
        endif  
        if(num .gt. 0) then
          norms(im) = sum / num
        else
          norms(im) = amiss
        endif
      end do
      return
      end
        
