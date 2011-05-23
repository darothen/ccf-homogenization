c***********************************************************************
c                          program normals_dist_2000.f
c                       from the ushcn program hcn_dist.f
c
c date of last change:  18 jan 2001
c 
c modifier:  claude williams
c
c version     modification                                     date    *
c -------     ------------                                  ---------- *
c  2004.v3    slight change in format and modified test for enough
c                stations                                   18 Dec 2008
c
c  2004.v2    skip the composite & candidate USHCN in neighbors
c                                                           15 Jul 2005
c  
c  2004       catch the 100 nearest neighbors               6 apr 2005
c
c  2000       version for the linux server                   1 feb 2005
c
c   1.1       sun workstation version                       10/18/1995 *
c                                                                      *
c   1.0       original                                                 *
c                                                                      *
c author:  claude williams                                             *
c                                                                      *
c modifiers:  pam hughes                                               *
c             david bowman                                             *
c                                                                      *
c usage:  hcn_dist                                                     *
c                                                                      *
c description:  this program determines the "numsrt" nearest neighbors *
c               to each normals candidate station.  this program is set up *
c               to handle up to "numstn" hcn stations, and only sorts  *
c               for the "numsrt" nearest neighbors.                    *
c               the major differences between this version and the 
c               hcn version are 
c             1) this routine uses separate lists for the candidates 
c                   and the network neighbors
c             2) reads from lat,lon lists instead of the ushcn station
c                   history file
c                                                                      *
c requirements:                                                        *
c                                                                      *
c   input file:  candidate lat,lon list
c                prospective neighbor ll list
c                                                                      *
c         **************************************************************
c notes:  * changes to any parameters must be made to the same         *
c         * parameters in any of the functions and/or subroutines.     *
c         **************************************************************
c                                                                      *
c         changed from calculating distances and sorting one candidate *
c         station at a time to doing all stations at one time.         *
c         - d. bowman:  01/21/94                                       *
c
c         changed back to doing one station at a time for flexibility
c         - c. williams 18 jan 2001
c                                                                      *
c references:  sun fortran reference guide, rev. a, 22 feb 1991        *
c                                                                      *
c              sun fortran user's guide, rev. a, 22 feb 1991           *
c                                                                      *
c              understanding fortran - michel boillet - 1978           *
c              isbn 0-8299-0355-0                                      *
c                                                                      *
c results:  produces file containing number of hcn stations, ids of    *
c           candidate station and nearest "numsrt" neighbors, pointers *
c           to nearest "numsrt" neighbors, and distances from          *
c           candidate to nearest "numsrt" neighbors.                   *
c                                                                      *
c   output file:  normals_distances
c                                                                      *
c parameters:                                                          *
c                                                                      *
c   numsrt     number of nearest neighbors to sort                     *
c                                                                      *
c   numstn     number of hcn stations arrays can hold                  *
c                                                                      *
c functions:  gcd                                                      *
c                                                                      *
c subroutines:  sortn                                                  *
c                                                                      *
c variables:                                                           *
c                                                                      *
c   dist       holds distances between candidate station and the
c              neighbor stations within "degsqr" degrees 
c   (numstn,numstn)                                                    *
c
c   degsqr     initially set to the inideg parameter (10 at present)
c              is the initial pass for the nieghbor stations
c
c   deginc     if not enough neighbor stations are withing degsqr then
c              it is incremented by this amount (2 at present)
c                                                                      *
c   fildir     directory where files are located                       *
c                                                                      *
c   i          loop counter                                            *
c                                                                      *
c   j          loop counter                                            *
c                                                                      *
c   lat        latitude                                                *
c   (numstn)                                                           *
c                                                                      *
c   latdeg     latitude degrees                                        *
c                                                                      *
c   latmin     latitude minutes                                        *
c                                                                      *
c   lon        longitude                                               *
c   (numstn)                                                           *
c                                                                      *
c   londeg     longitude degrees                                       *
c                                                                      *
c   lonmin     longitude minutes                                       *
c                                                                      *
c   nstn       number of hcn stations                                  *
c                                                                      *
c   ptr        index of stations sorted by closest distance of         *
c   (numstn,numstn)     "numsrt" nearest neighbors                     *
c                                                                      *
c   rcount     number of records read from station history file        *
c                                                                      *
c   stnid      station ids:       1 => station id                      *
c   (numstn,2) numstn stations    2 => position of last history record *
c                                                                      *
c   stnidr     station id from history record                          *
c                                                                      *
c***********************************************************************

      program hcndis
      include 'prehomog.parm.incl'  

      parameter ( maxnets = 17000, distini = 1000, 
     *  distinc = 500)

      parameter ( numsrt = ndist - 1 )
      
      real gcd,clat,clon
      real lat(maxnets),lon(maxnets)
      real dist(maxnets)
      integer stnid(maxnets), cstn, ptr(maxnets)
      character*6 comp(3)
      integer icomp(3)

c     command line variables
      character*132 argv, candfile, netfile, outfile
      integer cantype
c     input directory for default files (overridden with -i dirname)
c      originally - /'/sd4/dbowman/ushcn/files/'/
      integer iargc, narg

c     initialize variables
c     assume candidate format and neighbors are hcn
      cantype = 1
      reftype = 1
c     the distance from a station to itself is 0.0
      cdist = 0.0
      
      nc = 0
      
      candfile = ''
      netfile = ''
      outfile = ''

      print *,' generate nearest neighbor list - normals_dist_2004.v2'

c     begin command line code      
c     get the number of command line arguments
      narg = iargc()
      iarg = 1
c     go thru remaining command line arguments - keep the valid ones
c        and ignoring the undefined ones.
      do while(iarg .le. narg)
        call getarg(iarg, argv)
        iarg = iarg + 1
        if(argv .eq. '-c') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          candfile = argv
          print *, ' candidate station file :', candfile
        else if(argv .eq. '-o') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          outfile = argv
          print *, ' candidate w/network file :', outfile
        else if(argv .eq. '-n') then
          call getarg(iarg, argv)
          iarg = iarg + 1
          netfile = argv
          print *, ' network stations file :', netfile
        else if(argv .eq. '-h') then
          cantype = 0
          print *, ' enabled - non-hcn candidate file'
        else if(argv .eq. '-r') then
          reftype = 0
          print *, ' enabled - normals input from neighbor file'
        else if(argv .eq. '-u') then
          call getarg(iarg,argv)
          iarg = iarg + 1
          read(argv,fmt='(i4)') nc
          print *, ' Candidate-Neighbor merge - nc:',nc
        else
          print *, ' unknown argument :', argv, ': skipping'
        endif
      end do
      
          if(candfile.ne.'' .and. netfile.ne.'' .and. outfile.ne.'') 
     *  goto 11
      
      print *,' normals network neighborhood generation routine'
      print *,
     *  ' usage: normals_dist_2000 -c candfile -n netfile -o outfile'
      print *,'   required parameters'
      print *,'     -u nc         number of candidate stations'
      print *,'    input files ---'
      print *,'     -c candfile   candidate stations file'
      print *,'     -n netfile    network stations file'
      print *,'    output files ---'
      print *,'     -o outfile   candidate w/network file'
      print *,'   optional parameters'
      print *,'     -h           hcn type candidates'
      print *,'     -r           normals type neighbors'
      stop

   11 continue
      if(nc .eq. 0) print *, ' Candidate-Neighbor merge DISABLED'

c     read network station file
      open(10,file=netfile,status='old',err=240)

c     get list of station ids, lats & lons
      do nstn = 1, maxnets
        if(reftype .eq. 0) then
c         reading the neighbor station list in normals format
          read(10,*,end=20,err=240)idum,stnid(nstn),lat(nstn),lon(nstn)
        else
c         reading the neighbors in hcn format
          read(10,*,end=20,err=240) stnid(nstn), lat(nstn), lon(nstn)
        endif  
      end do

      print *,' too many network stations - increase maxnets'
      go to 360

c     save number of stations
   20 nstn = nstn - 1

c     close input file
      close(10,err=240)

c     open candidate station file
      open(10,file=candfile,status='old',err=250)

c     open output file
      open(16,file=outfile,err=260)

      istn = 0
      do while(1 .eq. 1)
        if(cantype .eq. 0) then
c         reading the non-hcn (coop) in hcn type format
          read(10,*,end=125,err=250) cstn, clat, clon
          do ic = 1, 3
            icomp(ic) = 0
          enddo  
        else
c         reading the candidates in hcn type format
          read(10,1000,end=125,err=250) cstn, clat, clon,  
     *      comp(1), comp(2), comp(3)
 1000     format(i6,f9.4,f10.4,5x,3(a6,1x))
c         updated for higher resolution meta lat/lon 
c 1000     format(i6,2f8.2,6x,3(a6,1x))
          do ic = 1,3
            if(comp(ic) .eq. '------') then
              icomp(ic) = 0
            else  
              read(comp(ic),'(i6)') icomp(ic)
            endif
          enddo    
        endif  
        istn = istn + 1

c        print *, cstn
        if(clat .gt. 90. .or. clat .lt. -90. .or. clon .gt. 360.
     *     .or. clon .lt. -180.) then
          print *,' cstn: ', cstn, ' lat/lon undefined: ', clat, clon
          go to 120
        endif  
          
c       compute distances between stations
        hidist = distini
        
   60   is = 0
        do 70 i = 1,nstn
c         skip the candidate station in the reference network
          if(stnid(i) .eq. cstn) go to 70
c         skip the composite station in the reference network
          do ic = 1,3
            if(stnid(i) .eq. icomp(ic)) goto 70
          enddo  
          adist = gcd(lat(i),lon(i),clat,clon)
          if(adist .lt. hidist) then
            is = is + 1
            dist(is) = adist
            ptr(is) = i
c            print *,stnid(i), adist, is
          endif  
   70   continue
   
c       see if there were enough stations, if not, crank up cut off
        if(is .lt. numsrt .and. is .lt. nstn-1) then
          print *, cstn, ' has too few neighbors within: ', hidist,
     *      ' try: ', hidist + distinc
          hidist = hidist + distinc
          go to 60
        endif  

c       sort distance and id arrays
        call sortn(dist,is,ptr)
        
c       output station ids, pointers, and distances
        write(16,90,err=260) cstn,(stnid(ptr(j)),j=1,numsrt)
   90   format(<numsrt+1>(i6.6,1x))

        write(16,100,err=260) istn,(ptr(j)+nc,j=1,numsrt)
  100   format(<numsrt+1>(i6,1x))

        write(16,110,err=260) cdist,(dist(j),j=1,numsrt)
  110   format(<numsrt+1>(f6.1,1x))
  120   continue
      end do

c     close input file
  125 close(10,err=250)

      if(nc .eq. 0) goto 225
c     ----- Now go through the Neighbors, adding to the cand-neigh list -----
      do istn = 1, nstn
        cstn = stnid(istn)
        clat = lat(istn)
        clon = lon(istn)

c        print *, cstn
        if(clat .gt. 90. .or. clat .lt. -90. .or. clon .gt. 360.
     *     .or. clon .lt. -180.) then
          print *,' cstn: ', cstn, ' lat/lon undefined: ', clat, clon
          go to 220
        endif  
          
c       compute distances between stations
        hidist = distini
        
  160   is = 0
        do 170 i = 1,nstn
c         skip the candidate station in the reference network
          if(stnid(i) .eq. cstn) go to 170
          adist = gcd(lat(i),lon(i),clat,clon)
          if(adist .lt. hidist) then
            is = is + 1
            dist(is) = adist
            ptr(is) = i
c            print *,stnid(i), adist, is
          endif  
  170   continue
   
c       see if there were enough stations, if not, crank up cut off
        if(is .lt. numsrt) then
          print *, cstn, ' has too few neighbors within: ', hidist,
     *      ' try: ', hidist + distinc
          hidist = hidist + distinc
          go to 160
        endif  

        print *, 'Cand: ', cstn, clat, clon, ' neigh: ', is

c       sort distance and id arrays
        call sortn(dist,is,ptr)
        
c       output station ids, pointers, and distances

        write(16,90,err=260) cstn,(stnid(ptr(j)),j=1,numsrt)

        write(16,100,err=260) istn+nc,(ptr(j)+nc,j=1,numsrt)

        write(16,110,err=260) cdist,(dist(j),j=1,numsrt)
  220   continue
      end do

c     close input file
  225 close(10,err=250)

c     close output file
      close(16,err=260)

c     notify user of successful completion
      write(6,230)
  230 format('hcn_dist has completed successfully!')

c     skip error messages and end program
      go to 380

c     error messages and aborts
  240 call perror ( "cannot open/read/close: " // netfile )
      go to 360

  250 call perror ( "cannot open/read/close: " // candfile )
      go to 360

  260 call perror ( "cannot open/read/close: " // outfile )

  360 write(6,370)
  370 format('program is aborting.')

  380 stop
      end

c***********************************************************************
c end of program hcndis.                                               *
c***********************************************************************


c***********************************************************************
c                             function gcd                             *
c                                                                      *
c date of last change:  21 january 1994                                *
c                                                                      *
c modifier:  david bowman                                              *
c                                                                      *
c version     modification                                     date    *
c -------     ------------                                  ---------- *
c   1.1       sun workstation version                       01/21/1994 *
c                                                                      *
c   1.0       original                                                 *
c                                                                      *
c author:  claude williams                                             *
c                                                                      *
c modifiers:  pam hughes                                               *
c             david bowman                                             *
c                                                                      *
c usage:  gcd(double precision alat1, double precision alon1,          *
c             double precision alat2, double precision alon2)          *
c                                                                      *
c description:  this program calculates the distance in kilometers     *
c               between two stations, given the latitude and longitude *
c               of each.                                               *
c                                                                      *
c notes:  none.                                                        *
c                                                                      *
c results:  returns distance in kilometers between two stations.       *
c                                                                      *
c variables:                                                           *
c                                                                      *
c   alat1      latitude of first station                               *
c                                                                      *
c   alat2      latitude of second station                              *
c                                                                      *
c   alon1      longitude of first station                              *
c                                                                      *
c   alon2      longitude of second station                             *
c                                                                      *
c   cd2r       degrees to radians constant                             *
c                                                                      *
c   radius     radius of earth in kilometers                           *
c                                                                      *
c   rdlat      intermediate distance calculation                       *
c                                                                      *
c   rdlon      intermediate distance calculation                       *
c                                                                      *
c   rolat      intermediate distance calculation                       *
c                                                                      *
c   temp       intermediate distance calculation                       *
c                                                                      *
c***********************************************************************

      function gcd(alat1,alon1,alat2,alon2)

c     initialize constants

      cd2r = 3.14159265 / 180.
      radius = 6371.

c     compute distance between two stations

      rolat = cd2r * (alat1 + alat2) / 2.
      rdlat = cd2r * (alat1 - alat2)
c     watchout for different hemispheres
      aldiff = alon1 - alon2
      if(aldiff .lt. -180.) aldiff = aldiff + 360.
      if(aldiff .gt. 180.) aldiff = aldiff - 360.
          
      rdlon = cd2r * aldiff
      temp = acos(cos(rdlat) * cos(rdlon * cos(rolat)))
      gcd = abs(temp) * radius

      return

      end

c***********************************************************************
c end of function gcd.                                                 *
c***********************************************************************


c***********************************************************************
c                           subroutine sortn                           *
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
c usage:  sortn(double precision dist(numstn,numstn), integer nstn,    *
c               integer nsort, integer ptr(numstn,numstn))             *
c                                                                      *
c description:  this subroutine determines and sorts the nsort nearest *
c               neighbors for each candidate stations.                  *
c                                                                      *
c notes:  changed from sorting one candidate station at a time to      *
c         sorting all candidates at one time.  subroutine now sorts    *
c         distance array as well as pointer array, and sorts only for  *
c         the nsort nearest neighbors. - d. bowman:  01/21/94          *
c        
c         changed back to one sort at a time, and sorts the whole 
c         array coming in.... 
c         changed sorting algorithm to heapsort (see
c             numerical recipes (p.231)   - c. williams   19jan2001
c                                                                      *
c results:  sorts the distance and pointer arrays for the nsort        *
c           nearest neighbors for each candidate station.              *
c                                                                      *
c parameters:                                                          *
c                                                                      *
c   numstn     number of hcn stations for arrays to hold               *
c                                                                      *
c variables:                                                           *
c                                                                      *
c   dist       holds distances between all hcn stations                *
c   (numstn,numstn)                                                    *
c                                                                      *
c   distmp     holds distance value when sorting                       *
c                                                                      *
c   i          loop counter                                            *
c                                                                      *
c   j          loop counter                                            *
c                                                                      *
c   k          loop counter                                            *
c                                                                      *
c   min        index of minimum value                                  *
c                                                                      *
c   nsort      number of stations to sort                              *
c                                                                      *
c   nstn       number of hcn stations                                  *
c                                                                      *
c   ptr        index of stations sorted by closest distance of nsort   *
c   (numstn,numstn)     nearest neighbors                              *
c                                                                      *
c   ptrtmp     holds pointer value when sorting                        *
c                                                                      *
c***********************************************************************

      subroutine sortn(dist,nstn,ptr)

      real dist(nstn)
      integer ptr(nstn)
      
      l = nstn/2 + 1
      ir = nstn
      
   10 continue
      if(l .gt. 1) then      
        l = l - 1
        rdist = dist(l)
        rptr = ptr(l)
      else
        rdist = dist(ir)
        rptr = ptr(ir)  
        dist(ir) = dist(1)
        ptr(ir) = ptr(1)
        ir = ir - 1
        if(ir .eq. 1) then
          dist(1) = rdist
          ptr(1) = rptr
          return
        endif
      endif
      i = l
      j = l+l
   20 if(j .le. ir) then
        if(j .lt. ir) then
          if(dist(j) .lt. dist(j+1)) j = j + 1
        endif
        if(rdist .lt. dist(j)) then
          dist(i) = dist(j)
          ptr(i) = ptr(j)
          i = j
          j = j + j
        else
          j = ir + 1
        endif
        go to 20
      endif
      
      dist(i) = rdist
      ptr(i) = rptr
      go to 10
      end    
 
c***********************************************************************
c end of subroutine sortn.                                             *
c***********************************************************************
