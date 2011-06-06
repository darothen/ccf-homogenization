c ---------------------------------------------------------------
c                         Subroutine readrsta
c
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
c
c   4p        rearranged directory i/o structure for proceedure type
c               instead of state directories                  18may09
c
c   doe2a     added 3-flag i/o option (-3)                    03mar09
c
c   doe2      added flag "g" for date to remove from fillin
c               due to GAP causing UNSTBL series              03apr07
c
c   doe       reverted to any number of years
c
c    3        limited fillin to 10 years from last good value 28sep06
c
c    2        added confidence interval reading & use        21nov2005
c             confidence flag redefined for fill_2004
c
c    1        Original                                       24Jan2001
c
c Author:   Claude Williams
c
c USAGE:
c
C DESCRIPTION: Read in the network reference stations meta data
c                used for adjustment/fillins.

      subroutine readrsta(netfile, iunit, irfmt,nc)

c     parameters from calling routine
      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.fill.incl'

      character*132 netfile
      character*6 cstn
      
c     read network station file list
      OPEN(iunit, FILE = netfile,STATUS = 'old',ERR = 150)

C     GET LIST OF reference STATION IDS, and datafile names
      DO 30 Inet = 1,NUMREF
        read(iunit,20,err=170,end=35) rid(inet)
   20   format(i6)
        write(cstn, '(i6.6)') rid(inet)
c       fill out the pathname with the appropriate directory early
c       nc indicates the number of candidates in the metafile from
c       the command line
        if(inet .le. nc) then
          nfname(inet) = candir(1:lnblnk(candir)) // 
     *      ctype(1:lnblnk(ctype)) // '/' // cstn // 
     *      icelem(1:lnblnk(icelem)) // '.' // ctype(1:lnblnk(ctype))
          ifname(inet) = candir(1:lnblnk(candir)) //
     *      itype(1:lnblnk(itype)) // '/' // cstn // 
     *      icelem(1:lnblnk(icelem)) // '.' // itype(1:lnblnk(itype))
        else
          nfname(inet) = netdir(1:lnblnk(netdir)) //
     *      ntype(1:lnblnk(ntype)) // '/' // cstn //
     *      icelem(1:lnblnk(icelem)) // '.' // ntype(1:lnblnk(ntype))
          ifname(inet) = netdir(1:lnblnk(netdir)) // 
     *      itype(1:lnblnk(itype)) // '/' // cstn //
     *      icelem(1:lnblnk(icelem)) // '.' // itype(1:lnblnk(itype))
        endif

   30 CONTINUE

   35 nstn = inet - 1

C     CLOSE network station FILE
      CLOSE(iunit)

      return
      
  150 print *,' Error in opening netfile: ', netfile
      stop 999
      
  170 print *,' Error in reading netfile: ', netfile
      stop 999
      
      end

c --------------------------------------------------------------
c                         Subroutine filinit
c
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
c  1.0        Original                                     18 Jan 2001
c
c Author:   Claude Williams
c
c USAGE:
C
C DESCRIPTION: This routine initializes the commons and other parameters
c                used throughout the Filnet subroutines
c      Mainly, the DEL array for the Network of stations around the
c        Candidate will be setup here so that they will not have to
c        be reread in the filentry routine
c
      subroutine filinit(cstn, ista, nnunit, idunit, igood, idbg, nel)

c     cstn is the desired station id
c     ista is the station's index

c     parameters from calling routine
      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.fill.incl'
c     parameters inherent to filnet_subs (inumyr = nyr)
      INCLUDE 'posthomog.comm.incl'  
      
c    --------------- COMMON /BLK3/

c   FILENTRY subroutine argument from calling program
C   DEBUG      Print Out Variable Values:      0 => No Debug Printouts *
C                                            > 0 => Debug Printouts    *

c   Following common variables are generated in TTFILL or RSFILL for
c     the optimum network (in common mainly for debug output in main)
C   IOPIK      Number of Neighboring stations consituting the optimum network
C   IOPSRT(inmstn)  Order of ISRT for optimum network
C   IOPTA      Number of Years in After Part of Window of Station with *
C              Lowest Confidence Interval                              *
C   IOPTB      Number of Years in Before Part of Window of Station     *
C              with Lowest Confidence Interval                         *
C   JY1OP      Optimized Value of JY1                                  *
C   JY2OP      Optimized Value of JY2                                  *
C   NAY(inmstn)  Number of After Years in Window for Each Network Station
C   NBY(inmstn)  Number of Before Years in Window for Each Network Station       *
C   OPBAN      Optimized Value of BANF Associated with Station with    *
C              Lowest Confidence Interval                              *
C   OPCCF      Optimized Value of OFF Associated with Station with     *
C              Lowest Confidence Interval                              *
C   OPCI       Optimized (Lowest) Value of CONF Associated with        *
C              Station with Lowest Confidence Interval                 *
C   OPCONF(inmstn)  CONF of Station Network at Time of Optimization         *
C   OPNET      Sum of Ratios of a Candidate's Neighbor's Monthly Data  *
C              for the Current Year to the Confidence Intervals at     *
C              Time of Optimization                                    *
C   OPNUM      Sum of Reciprocals of CONF at Time of Optimization      *

c   generated in FILENTRY wrt candidate data
C   LY         Current Year Being Processed

      DIMENSION IOPSRT(inmstn),NAY(inmstn),NBY(inmstn)
      DIMENSION OPCONF(inmstn)
      INTEGER DEBUG,OPBAN
      COMMON/BLK3/DEBUG,IOPIK,IOPSRT,IOPTA,IOPTB,JY1OP,JY2OP,LY,NAY,NBY,
     *            OPBAN,OPCCF,OPCI,OPCONF,OPNET,OPNUM,
     *            relnet,stdnet,pcpnet
      
c   ----------------- COMMON /BLK4/

c   generated in STNBIN
C   XJY(20)        Candidate's Neighbors' Precipitation Values             *
C   XKY(20)        Candidate's Precipitation Values                        *
      DIMENSION XJY(20),XKY(20)
      COMMON/BLK4/XJY,XKY
      
c   -----------------  COMMON /BLK5/
C   ICONF      Confidence Factor of the Network - Used only with Precipitation
c     as of the May 1996 FILNET this value has been hard coded to 1. 
      COMMON/BLK5/ICONF

      
c   -----------------  COMMON /BLK6/
c   generated in SRTDVR 
C   XSDIF(numrat)      Array of Precipitation Ratios 
      DIMENSION XSDIF(numrat)
      COMMON/BLK6/XSDIF

c   ----------------- COMMON /BLK7/
c   generated in FILCF and FILRS
C   CONF(inmstn)  Confidence Intervals for a Candidate's Nearest Neighbors
      DIMENSION CONF(inmstn)
      COMMON/BLK7/CONF

c   ----------------- common /netdata/
c     delm will contain the network data
      real DELM(inmstn+1,inumyr,13)
      character*6 candsta
      real crel(13), cstd(13), cpcp(13)
      common/netdata/ delm, crel, cstd, cpcp, candsta
      
      integer indata(nmth), imiss /-9999/
      integer ptr(inmstn), igood, idbg
      character*6 cstn, cnstn/'      '/, instn
      character*3 inflag(nmth)

      candsta = cstn

      do iy = 1, inumyr
        do imth = 1, 13
          do istn = 1, inmstn+1
            delm(istn,iy,imth) = AMISS
          enddo  
        end do
      end do
          
c     Set the element to the incoming parameter and input resolution
      inel = nel
      if(inel .ne. 4) then
        rdiv = 10.
      else
        rdiv = 100.
      endif    

c     ICONF is hard wired to 1 in the May 1996 FILNET
      iconf = 1
      
c     DEBUG printout value gets set here
      debug = idbg
      
C     GET STATION NETWORK pointers from the candidate/network file
C     Since we will be reading it each call, the cand/net file is opened 
c       in main
      do while(cstn .ne. cnstn)
        read(nnunit, 105,end=200) cnstn
  105   format(a6)      
        read(nnunit,110, end=200) (PTR(K),K=1,INMSTN)
  110   FORMAT(<INMSTN>(I6.6,1X),/)
      enddo  

c     for each reference station read in its data
  103 DO K = 2,inmstn
c       Check to see if neighbor not available
        if(ptr(k) .eq. 0) goto 150
          
c       Open temporary i/o unit for data input
        open(idunit, name = nfname(ptr(K)))
        do while (1 .eq. 1)
          if(inp3flag .eq. 0) then
c            original HCN format
c            read(idunit, '(a6,1x,i1,1x,i4,12f8.2)', end = 100)
c     *        instn, isrc, iyear, indata
            read(idunit, 1000, end = 100) instn, iyear, 
     *        (indata(imth),inflag(imth),imth=1,12)
 1000       format(a6,1x,i4,12(i6,a1))
            if(iyear .ge. IBYEAR .and. iyear .le. IEYEAR) then
              do imth = 1, NMTH
                if(inflag(imth).eq.'S' .or. inflag(imth).eq.'X') then
                  delm(k, iyear - IBYEAR + 1, imth) = AMISS                
                else if(indata(imth) .eq. IMISS) then
                  delm(k, iyear - IBYEAR + 1, imth) = AMISS
                else if(inflag(imth) .eq. 'd' .or. 
     *            inflag(imth) .eq. 'g') then
                  delm(k, iyear - IBYEAR + 1, imth) = AMISS                
                else  
                  delm(k, iyear - IBYEAR + 1, imth) = indata(imth)/rdiv
                endif  
              end do
            end if
          else
c           normals 3-flag format
            read(idunit,1500,end=100,err=100) instn,iyear,
     *        (indata(imth),inflag(imth),imth=1,12)
 1500       format(a6,1x,i4,12(i6,a3))
            if(iyear .ge. IBYEAR .and. iyear .le. IEYEAR) then
              do imth = 1, 12
                if(indata(imth) .eq. imiss) then
                  delm(k, iyear - IBYEAR + 1, imth) = AMISS
                else if(inflag(imth)(2:2) .ne. ' ') then
c                 any flag means month value failed a QC test
                  delm(k, iyear - IBYEAR + 1, imth) = AMISS                
                else  
                  delm(k, iyear - IBYEAR + 1, imth) = indata(imth)/rdiv
                endif
              end do 
            endif   
          endif
        end do
c       make sure to close this unit!!!
  100   close(idunit)
        igood = 1
      end do  
      return
      
  150 if(k .eq. 2) then
        print *,' Station: ', cstn, ' has no stations in network'
        igood = 0
      else
        print *,' Station: ',cstn,' has only ',k-1,' neighbors'
        igood = 1
      endif
      return  
      
  200 print *,' Station: ', cstn, ' is not equal ', cnstn,
     *  ' has no network or not sorted'
      igood = 0
      end

C***********************************************************************
c                         Subroutine filentry
c
c 
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
c  1.0        Original                                     18 Jan 2001
c
c Author:   Claude Williams
c
c USAGE:
C                                                                      *
C DESCRIPTION: This routine is the wrapper for the Filnet subroutines
c                used with the 1971-2000 Normals
c              Input is the Candidate station and its Network neighbors
c              For the testing phase the Candidate is an HCN station
c                 with data artificial removed (see rnum.HCNgaps.f)
c
c     Need to define LY, IM, and INEL
c
      subroutine filentry(cdata,dconf,dflag,cflag,ngap,inmth,inyr)

c     parameters from calling routine
      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.fill.incl'  
c     parameters internal to filnet_subs
      INCLUDE 'posthomog.comm.incl'  
 
c    --------------- COMMON /BLK3/
      DIMENSION IOPSRT(inmstn),NAY(inmstn),NBY(inmstn)
      DIMENSION OPCONF(inmstn)
      INTEGER DEBUG,OPBAN
      COMMON/BLK3/DEBUG,IOPIK,IOPSRT,IOPTA,IOPTB,JY1OP,JY2OP,LY,NAY,NBY,
     *            OPBAN,OPCCF,OPCI,OPCONF,OPNET,OPNUM,
     *            relnet,stdnet,pcpnet

c   ----------------- common /netdata/
c     delm will contain the network data
      real DELM(inmstn+1,inumyr,13)
      character*6 candsta
      real crel(13), cstd(13), cpcp(13)
      common/netdata/ delm, crel, cstd, cpcp, candsta
      
c     cdata is the candidate station's data 
      real cdata(12, inumyr), dconf(12, inumyr)
      CHARACTER*40 STALST
      
      integer ngap(NMTH)
      character*1 dflag(nmth, inumyr), cflag(nmth,nyr)
      
      do imth = 1, nmth
        ngap(imth) = 0
      end do  
      
c     move the incoming candidate data into the delm array
      do ly = 1, nyr
        do imth = 1, nmth
          delm(1, ly, imth) = cdata(imth, ly)
        end do
      end do  
      
c     assuming that filinit has taken care of the rest of the network
c       find the missing data, sending each instance to the 
c       appropriate fillin routine (TTFILL for temp, RSFILL for precip)
c     if inyr and/or inmth are non-zero only test fillin for those yr/mth
      if(inyr .eq. 0) then
        ly1 = 1
        ly2 = nyr
      else
        ly1 = inyr
        ly2 = inyr
      endif
      if(inmth .eq. 0) then
        im1 = 1
        im2 = nmth
      else
        im1 = inmth
        im2 = inmth
      endif    
      do ly = ly1, ly2
        do im = im1, im2 
          if(delm(1, ly, im) .le. AMISS + 1.0) then
            if(debug .ge. 2) print *,' Filling year:', ly, ' mth:', im
            if(ly+ibyear-1 .gt. 1970) ngap(im) = ngap(im) + 1
            DINT = 99.99
C           DETERMINE "OPTIMUM HOMOGENEOUS NETWORK" TO DETERMINE 
C           MISSING DATUM INITIALIZE VARIABLES
            JY1OP = 0
            JY2OP = 0
            OPBAN = 0
            OPCCF = 0.
            OPCI = 99.99

c           set the network before and after window years at maximum
            do is = 2, inmstn
              nby(is) = ly
              nay(is) = nyr - ly
            end do  

            IF(INEL .LT. 4) THEN
              CALL TTFILL(DELm(1, 1, im))
            ELSE
              CALL RSFILL(DELm(1, 1, im))
            END IF
            
C           INITIALIZE LIST OF NEIGHBORS
            STALST = '1                                       '

c           pass through network information
            relnet = crel(im)
            stdnet = cstd(im)
            pcpnet = cpcp(im)

C           MARK STATIONS USED IN NETWORK
            DO 110 IT = 1,IOPIK
              STALST(IOPSRT(IT):IOPSRT(IT)) = '*'
  110       CONTINUE

C           SET CONFIDENCE INTERVAL AND DATA VALUE TO MISSING IF OPTIMIZED
C           CONFIDENCE INTERVAL IS 99.99
  120       iym = ly * 100 + im
            lastym = 11203
            IF(OPCI .EQ. 99.99 .and. ly .le. lastym) THEN
              dconf(im, ly) = amiss
              cdata(im, ly) = AMISS
              dflag(im, ly) = 'X'
              cflag(im, ly) = ' '
              if(debug .ge. 1) 
     *           print *,' WARNING: FILNET FAILED year: ',ly,' mth:',im

            ELSE IF(IOPIK .GE. 2) THEN
C             CALCULATE DATA VALUE AND CONFIDENCE INTERVAL
C             MUST HAVE AT LEAST TWO NEIGHBORING STATIONS
              IF(INEL .LE. 3) THEN
                cdata(im, ly) = OPCCF + OPNET / OPNUM
                dflag(im, ly) = 'M'
                dint = dconf(im,ly)
                IF(dint .EQ. 99.99) THEN
                  dconf(im,ly) = OPCI
                ELSE
                  dconf(im,ly) = SQRT(dint*dint + OPCI*OPCI)
                END IF

C             PRECIPITATION MUST HAVE AT LEAST FIVE NEIGHBORING STATIONS OR
C             NUMBER OF STATIONS PLUS WINDOW TIME PERIOD MUST BE AT LEAST 10;
C             E.G., TWO NEIGHBORS MUST HAVE AT LEAST EIGHT YEARS; THREE
C             NEIGHBORS MUST HAVE AT LEAST SEVEN YEARS; ETC.
              ELSE IF(IOPIK .GE. 5 .OR.
     *          (JY2OP - JY1OP + IOPIK) .GE. 10) THEN
                cdata(im, ly) = EXP(OPCCF) * OPNET / OPNUM
                dflag(im, ly) = 'M'
                cflag(im, ly) = 'M'

C               DETERMINE IF PRECIPITATION VALUE SHOULD BE 0.
                IF(cdata(im, ly) .LT. 0.014997) THEN
                  PCPCNT = 0
                  DO 130 III = 1,IOPIK
c                   Fix reversal of ioptb and iopta - cw 10aug01
                    DO 130 JJJ = LY-IOPTB+1,LY+IOPTA
                      IF(DELM(IOPSRT(III),JJJ,im) .EQ. 0.0)
     *                  PCPCNT=PCPCNT+1
  130             CONTINUE
                  IF((2*PCPCNT) .GT. (IOPIK*(IOPTB+IOPTA)))
     *              cdata(im, ly) = 0. 
                END IF

C               CALCULATE MEAN OF OPTIMIZED NETWORK
C               INITIALIZE VARIABLES
                AVGOPN = 0.
                DATMAX = 0.
                DATMIN = 99.99

C               CALCULATE MEAN; FIND MAXIMUM AND MINIMUM VALUES
                DO 140 IT = 1,IOPIK
                  AVGOPN = AVGOPN + DELm(IOPSRT(IT),LY,im)
                  IF(DELm(IOPSRT(IT),LY,im) .GT. DATMAX)
     *              DATMAX = DELm(IOPSRT(IT),LY,im)
                  IF(DELm(IOPSRT(IT),LY,im) .LT. DATMIN)
     *              DATMIN = DELm(IOPSRT(IT),LY,im)
  140           CONTINUE
  
C               CALCULATE MEAN
                AVGOPN = AVGOPN / FLOAT(IOPIK)

C               CALCULATE RANGE
                RANGE = DATMAX - DATMIN

C               ENSURE THAT THE ADJUSTED VALUE IS NO MORE THAN 1.5 TIMES THE MEAN
C               + 1/2 THE RANGE, OR NO LESS THAN 2/3 THE MEAN - 1/2 THE RANGE
                CALC = 1.5 * AVGOPN + 0.5 * RANGE
                IF(cdata(im, ly) .GT. CALC) THEN
                  cdata(im, ly) = CALC
                ELSE
                  CALC = 2./3. * AVGOPN - 0.5 * RANGE
                  IF(CALC .LT. 0.) CALC = 0.  
                  IF(cdata(im, ly) .LT. CALC) cdata(im, ly) = CALC
                END IF
              END IF

C             DETERMINE CONFIDENCE LEVEL
              dint = dconf(im,ly)
              IF(DINT .EQ. 99.99) THEN
                dconf(im,ly) = OPCI
              ELSE
                dconf(im,ly) = SQRT(DINT * DINT + OPCI * OPCI)
              END IF
            END IF

C           DEBUG PRINTOUT

            IF(DEBUG .Ge. 1) WRITE(6,150) candsta, IbYEAR+LY-1,im,
     *         cdata(im, ly),
     *         IBYEAR-JY1OP+1,IBYEAR+JY2OP-1,OPCCF,
     *         OPNET/OPNUM,DINT,OPCI,
     *         dconf(im,ly),STALST
  150       FORMAT(a6,I5,i3,1X,F6.2,1X,2(I4,1X),2(F6.2,1X),
     *         3(F6.2,1X),A40)

          end if
        end do
      end do
      end  

C***********************************************************************
C                          Subroutine TTFILL                           *
C                                                                      *
C DATE OF LAST CHANGE:  10 August 1994                                 *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         08/10/94 *
C                                                                      *
C   1.0       Original Version                                01/02/86 *
C                                                                      *
C AUTHORS:  Claude Williams                                            *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  TTFILL(real DEL(inmstn+1,inumyr))                            *
C                                                                      *
C DESCRIPTION:  This subroutine determines the optimum nearest         *
C               neighbor interval for the correction of data with      *
C               respect to station moves, using confidence intervals   *
C               associated with the T-Test.                            *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Determines the optimum nearest neighbor interval for the   *
C           correction of data with respect to station moves, using    *
C           confidence intervals associated with the T-Test.           *
C                                                                      *
C PARAMETERS:                                                          *
C                                                                      *
C   inmstn     Network of Stations Including Candidate                 *
C                                                                      *
C   inumyr     Last Cell of Arrays that Hold Year Data; Currently Set  *
C              to Handle Through the Year 2025                         *
C                                                                      *
C   minsta     Minimum Number of Stations for Window Not to be Broken  *
C                                                                      *
C   minyr      Minimum Number of Years in Window for Homogeneity       *
C                                                                      *
C SUBROUTINES:  CALCTA                                                 *
C               FILCF                                                  *
C               MRGBN                                                  *
C               SRTCNF                                                 *
C                                                                      *
C VARIABLES:
c
c  for common variables see filinit
C                                                                      *
C   B2I        Sum of Temperature Differences Squared Between          *
C              Candidate and Neighbors Over Years in Window            *
C                                                                      *
C   BAN        Total Number of Non-Missing Seasonal Values in Window   *
C              Less One (Degrees of Freedom)                           *
C                                                                      *
C   BANF       Array of Indices for Determining Value of T-Alpha from  *
C   (inmstn)   the Confidence Array RVT67 for a Candidate's Nearest    *
C              Neighbors                                               *
C                                                                      *
C   BI         Sum of Temperature Differences Between Candidate and    *
C              Neighbors Over Years in Window                          *
C                                                                      *
C   BN         Count of Non-Missing Candidate's and Neighbors' Monthly *
C              Values that are Not Flagged as "X"/"U" in Window        *
C                                                                      *
C   DEL        Monthly Data                          -99.99 => Missing *
C   (inmstn+1,inumyr)  inmstn+1 Stations, inumyr Years                 *
C                                                                      *
C   IK         Loop Counter                                            *
C                                                                      *
C   IK1        Holds Station Index Value                               *
C                                                                      *
C   INDTA      Index for Determining Value of T-Alpha from the         *
C              Confidence Interval Array RVT67                         *
C                                                                      *
C   IS         Loop Counter                                            *
C                                                                      *
C   ISRT       Array of Index Numbers for CONF, Which Sorts CONF from  *
C   (inmstn)   Lowest to Highest                                       *
C                                                                      *
C   IT         Loop Counter                                            *
C                                                                      *
C   ITRIG      Determines if Window Limits Changed:                    *
C                                     0 => No Changes to Window Limits *
C                                     1 => Changes to Window Limits    *
C                                                                      *
C   JY1        Beginning of Window                                     *
C                                                                      *
C   JY2        End of Window                                           *
C                                                                      *
C   LEDIT      Index of Neighboring Station that is the Least          *
C              Correlated                                              *
C                                                                      *
C   MINBAN     Minimum Value of BANF for Neighboring Stations          *
C                                                                      *
C   OFF        Array of Confidence Correction Factors for a            *
C   (inmstn)   Candidate's Nearest Neighbors                           *
C                                                                      *
C   TA         Critical Value of T-Alpha from Confidence Interval      *
C              Array Used in Calculating the Confidence Interval       *
C                                                                      *
C   WCI        Confidence Interval                                     *
C                                                                      *
C***********************************************************************

      SUBROUTINE TTFILL(DEL)

      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.comm.incl'  

      DIMENSION DEL(inmstn+1,inumyr)

      COMMON/BLK3/DEBUG,IOPIK,IOPSRT,IOPTA,IOPTB,JY1OP,JY2OP,LY,NAY,NBY,
     *            OPBAN,OPCCF,OPCI,OPCONF,OPNET,OPNUM,
     *            relnet,stdnet,pcpnet
      COMMON/BLK7/CONF

      DIMENSION CONF(inmstn),IOPSRT(inmstn),NAY(inmstn),NBY(inmstn)
      DIMENSION OPCONF(inmstn)
      INTEGER DEBUG,OPBAN

      DIMENSION ISRT(inmstn),OFF(inmstn)
      DOUBLE PRECISION B2I
      INTEGER BANF(inmstn)

C     INITIALIZE VARIABLE

      IOPSTA = 0

C     INITIALIZE ALL STATIONS CONFIDENCE INTERVALS AND CONFIDENCE
C     CORRECTION FACTORS
   10 CONF(1) = 99.99
      OFF(1) = 0.

      DO 20 IS = 2,inmstn
        CONF(IS) = 99.99
        OFF(IS) = 0.

C       CHECK FOR DATA NOT MISSING
        IF(DEL(IS,LY) .gt. -99) THEN

C         CALCULATE CONFIDENCE INTERVALS IF ENOUGH YEARS
          IF((NAY(IS)+NBY(IS)) .GE. minyr)
     *       CALL FILCF(BANF,DEL,IS,LY,NAY(IS),NBY(IS),OFF,WCI)
        END IF
   20 CONTINUE

C     SORT CONFIDENCE INTERVALS FROM LOWEST TO HIGHEST
   30 CALL SRTCNF(ISRT)

C     INITIALIZE OPTIMIZATION VARIABLES TO VALUES FROM STATION WITH
C     LOWEST CONFIDENCE INTERVAL
      IF(OPCI .EQ. 99.99) THEN
        IOPIK = 1
        IOPSRT(1) = ISRT(1)
        OPBAN = BANF(ISRT(1))   
        OPCCF = OFF(ISRT(1))
        OPCI = CONF(ISRT(1))
        OPNET = DEL(ISRT(1),LY)
        OPNUM = 1.
      END IF

      IOPTA = NAY(ISRT(1))
      IOPTB = NBY(ISRT(1))
      LOPTA = IOPTA
      LOPTB = IOPTB

C     DEBUG PRINTOUT
      IF(DEBUG .GE. 2) WRITE(6,40) BANF,CONF,ISRT,NAY,NBY
   40   FORMAT('BANF= ',<inmstn>i6,/,'CONF= ',<inmstn>F6.3/,
     *    'ISRT= ',<inmstn>I6,/,'NAY = ',<inmstn>I6,/,
     *    'NBY = ',<inmstn>I6)

C     CONTINUE IF CANDIDATE STATION DOES NOT HAVE THE LOWEST CONFIDENCE
C     INTERVAL
      IF(ISRT(1) .NE. 1) THEN

C       HELP STABILIZE STATION NETWORK
        DO 130 IK = 2,inmstn-1

C         RETURN IF CANDIDATE STATION
          IF(ISRT(IK) .EQ. 1) GO TO 140

          IK1 = IK

C         RESET NUMBER OF BEFORE/AFTER YEARS IF LESS THAN OPTIMUM STATION'S
C         YEARS
          IF(NBY(ISRT(IK)) .LT. IOPTB) THEN
            IK1 = 1
            IOPTB = NBY(ISRT(IK))
          END IF

          IF(NAY(ISRT(IK)) .LT. IOPTA) THEN
            IK1 = 1
            IOPTA = NAY(ISRT(IK))
          END IF

C         RESET BEGINNING AND END OF WINDOW
          JY1 = LY - IOPTB + 1
          JY2 = LY + IOPTA

C         INITIALIZE STATION VARIABLES AND RE-SORT CONFIDENCE INTERVALS IF
C         TOO FEW YEARS
          IF((IOPTA+IOPTB) .LT. minyr) THEN

C           DEBUG PRINTOUT
            IF(DEBUG .GE. 2) WRITE(6,50) IK,LY,JY1,JY2
   50                        FORMAT('IK = ',I6,' LY = ',I6,' JY1 = ',I6,
     *                              ' JY2 = ',I6)
            CONF(ISRT(IK)) = 99.99
            NAY(ISRT(IK)) = 0
            NBY(ISRT(IK)) = 0
            OFF(ISRT(IK)) = 0.
            GO TO 30
          END IF

C         RECALCULATE CONFIDENCE INTERVAL AND CONFIDENCE CORRECTION FACTOR
C         FOR CURRENT STATION (IK1 = IK), OR FOR ALL STATIONS FROM THE FIRST
C         THROUGH THE CURRENT STATION (IK1 = 1)
          DO 70 IT = IK1,IK
            IS = ISRT(IT)
            CALL FILCF(BANF,DEL,IS,LY,IOPTA,IOPTB,OFF,WCI)
            IF(DEBUG .ge. 3) WRITE(6,60) IS,IT,WCI,CONF(IS),OFF(IS)
   60                        FORMAT('IS = ',I6,' IT = ',I6,' WCI = ',
     *                              F10.5,' CONF(IS) = ',F10.5,
     *                              ' OFF(IS) = ',F10.5)
   70     CONTINUE
   
c          print *,'jy2 ', jy2

          CALL MRGBN(B2I,BI,BN,DEL,IK,ISRT,JY1,JY2,LY)

C         DEBUG PRINTOUT
          IF(DEBUG .ge. 3) WRITE(6,80) B2I,BI,BN,JY1,JY2
   80                      FORMAT('B2I = ',F10.5,' BI = ',F10.5,
     *                            ' BN = ',F10.5,' JY1 = ',I6,' JY2 = ',
     *                            I6)

C         CALCULATE BAN
          BAN = BN - 1.

C         FOLLOWING CODE IS TO HELP STABILIZE NETWORK STATION SELECTION
          IF(nint(BN) .LT. minyr) THEN

            ITRIG = 0

C           SET THE WINDOWS OF ALL STATIONS WITH WINDOWS WIDER THAN THE
C           CURRENT OPTIMIZED VALUE TO THE CURRENT OPTIMIZED VALUE
            DO 90 IS = 2,inmstn
              IF(NBY(IS) .GT. IOPTB) THEN
                ITRIG = 1
                NBY(IS) = IOPTB
              END IF

              IF(NAY(IS) .GT. IOPTA) THEN
                ITRIG = 1
                NAY(IS) = IOPTA
              END IF
   90       CONTINUE

C           NO CHANGES MADE TO WINDOW LIMITS
            IF(ITRIG .EQ. 0) THEN

C             INITIALLY BIAS TOWARD LONG PERIOD (AND FEWER STATIONS)
C             RETURN IF AT LEAST MINIMUM STATIONS
              IF(IK .GE. minsta) GO TO 140

C             FIND LOWEST CORRELATED STATION
              LEDIT = 1
              MINBAN = BANF(ISRT(1))

              DO 100 IT = 2,IK
                IF(BANF(ISRT(IT)) .LE. MINBAN) THEN
                  LEDIT = IT
                  MINBAN = BANF(ISRT(IT))
                END IF
  100         CONTINUE

C             REMOVE LOWEST CORRELATED STATION FROM THE ALGORITHM
              NAY(ISRT(LEDIT)) = 0
              NBY(ISRT(LEDIT)) = 0
            END IF

C           START SUBROUTINE OVER FROM THE BEGINNING
            GO TO 10
          END IF

C         SAVE OPTIMAL WINDOW VALUES
          LOPTA = IOPTA
          LOPTB = IOPTB

C         CALCULATE CONFIDENCE INTERVAL
          CALL CALCTA(BAN,INDTA,TA)
          
          if((BN * B2I - BI * BI) .le. 0) then
            wci = 0.01
          else  
            WCI = (TA / BN) * SQRT((BN * B2I - BI * BI) / (BN - 1.))
          endif  
          
C         STORE OPTIMIZED VALUES
          IF(WCI .LE. OPCI .OR. IOPSTA .LE. minsta) THEN

C           DEBUG PRINTOUT
            IF(DEBUG .GE. 2) 
     *        WRITE(6,110) WCI,BI/BN,IOPTB,IOPTA,
     *          (ISRT(IT),CONF(ISRT(IT)),OFF(ISRT(IT)),IT=1,IK)
  110         FORMAT('WCI: ',F7.4,' OFF: ',F7.3,
     *                              ' IOPTB: ',I2,' IOPTA: ',I2,
     *                              ' (STA,CONF,OFF):',
     *                              3(/,7(I3,2(',',F7.3),': ')))
            IOPIK = IK
            IOPSTA = IK
            JY1OP = JY1
            JY2OP = JY2
            OPBAN = INDTA
            OPCCF = BI / BN
            OPCI = WCI

C           INITIALIZE VALUES USED TO COMPUTE OPTIMIZED WEIGHTED AVERAGE
            OPNET = 0.
            OPNUM = 0.

C           COMPUTE OPTIMIZED WEIGHTED AVERAGE
            DO 120 IT = 1,IK
              IS = ISRT(IT)
              IOPSRT(IT) = IS
              OPCONF(IS) = CONF(IS)
              OPNET = OPNET + DEL(IS,LY) / CONF(IS)
              OPNUM = OPNUM + 1. / CONF(IS)
  120       CONTINUE
            if(debug.ge.2) print *,' Num sta:', iopik, ' Estimate:', 
     *        OPCCF + OPNET / OPNUM
          END IF
c         if not enough stations - STOP!
          if(ik .gt. minsta .and. ik .gt. iopik + 2) go to 140
  130   CONTINUE
      END IF

  140 RETURN

      END

C***********************************************************************
C End of Subroutine TTFILL.                                            *
C***********************************************************************

C***********************************************************************
C                          Subroutine RSFILL                           *
C                                                                      *
C DATE OF LAST CHANGE:  08 September 1994                              *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         09/08/94 *
C                                                                      *
C   1.0       Original Version                                01/02/86 *
C                                                                      *
C AUTHORS:  Claude Williams                                            *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  RSFILL(real DEL(inmstn+1,inumyr))                            *
C                                                                      *
C DESCRIPTION:  This subroutine determines the optimum nearest         *
C               neighbor interval for the correction of data with      *
C               respect to station moves, using confidence intervals   *
C               associated with the rank sum test.                     *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Determines the optimum nearest neighbor interval for the   *
C           correction of data with respect to station moves, using    *
C           confidence intervals associated with the rank sum test.    *
C                                                                      *
C PARAMETERS:                                                          *
C                                                                      *
C   inmstn     Network of Stations Including Candidate                 *
C                                                                      *
C   inumyr     Last Cell of Arrays that Hold Year Data; Currently Set  *
C              to Handle Through the Year 2025                         *
C                                                                      *
C   minsta     Minimum Number of Stations for Window Not to be Broken  *
C                                                                      *
C   minyr      Minimum Number of Years in Window for Homogeneity       *
C                                                                      *
C SUBROUTINES:  CALCZA                                                 *
C               FILRS                                                  *
C               MRGBN                                                  *
C               SRTCNF                                                 *
C               SRTDVR                                                 *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
c  for common variables see filinit
c
C   B2I        Not Used by This Subroutine                             *
C                                                                      *
C   BAN        BN * BN                                                 *
C                                                                      *
C   BANF       Values of BAN for Each Network Station                  *
C   (inmstn)                                                           *
C                                                                      *
C   BI         Not Used by This Subroutine                             *
C                                                                      *
C   BN         Count of Non-Missing Candidate's and Neighbors' Monthly *
C              Values that are Not Flagged as "X"/"U" in Window        *
C                                                                      *
C   DEL        Monthly Data                          -99.99 => Missing *
C   (inmstn+1,inumyr)  inmstn+1 Stations, inumyr Years                 *
C                                                                      *
C   IK         Loop Counter                                            *
C                                                                      *
C   IK1        Holds Station Index Value                               *
C                                                                      *
C   IS         Loop Counter                                            *
C                                                                      *
C   ISRT       Array of Index Numbers for CONF, Which Sorts CONF from  *
C   (inmstn)   Lowest to Highest                                       *
C                                                                      *
C   IT         Loop Counter                                            *
C                                                                      *
C   ITRIG      Determines if Window Limits Changed:                    *
C                                     0 => No Changes to Window Limits *
C                                     1 => Changes to Window Limits    *
C                                                                      *
C   JY1        Beginning of Window                                     *
C                                                                      *
C   JY2        End of Window                                           *
C                                                                      *
C   LEDIT      Index of Neighboring Station that is the Least          *
C              Correlated                                              *
C                                                                      *
C   MINBAN     Minimum Value of BANF for Neighboring Stations          *
C                                                                      *
C   OFF        Array of Confidence Correction Factors for a            *
C   (inmstn)   Candidate's Nearest Neighbors                           *
C                                                                      *
C   WCCF       Confidence Correction Factor                            *
C                                                                      *
C   WCI        Confidence Interval                                     *
C                                                                      *
C***********************************************************************

      SUBROUTINE RSFILL(DEL)

      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.comm.incl'  

      DIMENSION DEL(inmstn+1,inumyr)

      COMMON/BLK3/DEBUG,IOPIK,IOPSRT,IOPTA,IOPTB,JY1OP,JY2OP,LY,NAY,NBY,
     *            OPBAN,OPCCF,OPCI,OPCONF,OPNET,OPNUM,
     *            relnet,stdnet,pcpnet
      COMMON/BLK7/CONF

      DIMENSION CONF(inmstn),IOPSRT(inmstn),NAY(inmstn),NBY(inmstn)
      DIMENSION OPCONF(inmstn)
      INTEGER DEBUG,OPBAN

      DOUBLE PRECISION B2I
      DIMENSION ISRT(inmstn),OFF(inmstn)
      INTEGER BANF(inmstn)

C     INITIALIZE ALL STATIONS CONFIDENCE INTERVALS AND CONFIDENCE
C     CORRECTION FACTORS

   10 CONF(1) = 99.99
      OFF(1) = 0.
      WCI = 99.99

      DO 20 IS = 2,inmstn
        CONF(IS) = 99.99
        OFF(IS) = 0.

C     CALCULATE CONFIDENCE INTERVALS IF DATA NOT MISSING AND ENOUGH
C     YEARS

        IF(DEL(IS,LY) .gt. -99. .AND. (NAY(IS) + NBY(IS)) .GT. minyr)
     *     CALL FILRS(BANF,DEL,IS,LY,NAY(IS),NBY(IS),OFF,WCI)
   20 CONTINUE

C     SORT CONFIDENCE INTERVALS FROM LOWEST TO HIGHEST

   30 CALL SRTCNF(ISRT)

C     CONTINUE IF CANDIDATE STATION DOES NOT HAVE THE LOWEST CONFIDENCE
C     INTERVAL

      IF(ISRT(1) .NE. 1) THEN

C     INITIALIZE OPTIMIZATION VARIABLES TO VALUES FROM STATION WITH
C     LOWEST CONFIDENCE INTERVAL

        IOPTA = NAY(ISRT(1))
        IOPTB = NBY(ISRT(1))
        LOPTA = IOPTA
        LOPTB = IOPTB
        OPNET = DEL(ISRT(1),LY)
        OPNUM = 1.

C     DEBUG PRINTOUT

        IF(DEBUG .GE. 2) WRITE(6,40) BANF,CONF,ISRT,NAY,NBY
   40                    FORMAT('BANF = ',<inmstn>(I3,1X),/,
     *                          'CONF = ',<inmstn>(F10.6,1X),/,
     *                          'ISRT = ',<inmstn>(I6,1X),/,
     *                          'NAY = ',<inmstn>(I6,1X),/,
     *                          'NBY = ',<inmstn>(I6,1X))

C     HELP STABLIZE STATION NETWORK

        DO 130 IK = 2,inmstn-1

C     RETURN IF CANDIDATE STATION

          IF(ISRT(IK) .EQ. 1) GO TO 140

C     RESET NUMBER OF BEFORE/AFTER YEARS IF LESS THAN OPTIMUM STATION'S
C     YEARS

          IF(NAY(ISRT(IK)) .LT. IOPTA) IOPTA = NAY(ISRT(IK))
          IF(NBY(ISRT(IK)) .LT. IOPTB) IOPTB = NBY(ISRT(IK))

C     RESET BEGINNING AND END OF WINDOW

          JY1 = LY - IOPTB + 1
          JY2 = LY + IOPTA

C     INITIALIZE STATION VARIABLES AND RE-SORT CONFIDENCE INTERVALS IF
C     TOO FEW YEARS

          IF((JY2-JY1) .LT. minyr) THEN

C     DEBUG PRINTOUT

            IF(DEBUG .GE. 2) WRITE(6,50) IK,LY,JY1,JY2
   50                        FORMAT('IK = ',I6,' LY = ',I6,' JY1 = ',I6,
     *                              ' JY2 = ',I6)
            CONF(ISRT(IK)) = 99.99
            NAY(ISRT(IK)) = 0
            NBY(ISRT(IK)) = 0
            OFF(ISRT(IK)) = 0.
            GO TO 30
          END IF

C     DETERMINE IF RECALCULATE CONFIDENCE INTERVAL AND CORRECTION FACTOR
C     FOR CURRENT STATION OR FOR ALL PRIOR STATIONS

          IF(IOPTA .EQ. LOPTA .AND. IOPTB .EQ. LOPTB) THEN
            IK1 = IK
          ELSE
            IK1 = 1
            LOPTA = IOPTA
            LOPTB = IOPTB
          END IF

C     RECALCULATE CONFIDENCE INTERVAL AND CONFIDENCE CORRECTION FACTOR
C     FOR CURRENT STATION (IK1 = IK), OR FOR ALL STATIONS FROM THE FIRST
C     THROUGH THE CURRENT STATION (IK1 = 1)

          DO 70 IT = IK1,IK
            IS = ISRT(IT)
            CALL FILRS(BANF,DEL,IS,LY,IOPTA,IOPTB,OFF,WCI)

C     DEBUG PRINTOUT

            IF(DEBUG .ge. 2) WRITE(6,60) IS,IT,WCI,CONF(IS),OFF(IS)
   60                        FORMAT('IS = ',I6,' IT = ',I6,' WCI = ',
     *                              F10.5,' CONF(IS) = ',F10.5,
     *                              ' OFF(IS) = ',F10.5)
   70     CONTINUE

          CALL MRGBN(B2I,BI,BN,DEL,IK,ISRT,JY1,JY2,LY)

C     DEBUG PRINTOUT

          IF(DEBUG .ge. 3) WRITE(6,80) BN,JY1,JY2
   80                      FORMAT('BN = ',F4.0,' JY1 =',I4,' JY2 =',I4)

C     FOLLOWING CODE IS TO HELP STABILIZE NETWORK STATION SELECTION

          IF(nint(BN) .LT. minyr) THEN

            ITRIG = 0

C     SET THE WINDOWS OF ALL STATIONS WITH WINDOWS WIDER THAN THE
C     CURRENT OPTIMIZED VALUE TO THE CURRENT OPTIMIZED VALUE

            DO 90 IS = 2,inmstn
              IF(NAY(IS) .GT. IOPTA) THEN
                ITRIG = 1
                NAY(IS) = IOPTA
              END IF

              IF(NBY(IS) .GT. IOPTB) THEN
                ITRIG = 1
                NBY(IS) = IOPTB
              END IF
   90       CONTINUE

C     NO CHANGES MADE TO WINDOW LIMITS

            IF(ITRIG .EQ. 0) THEN

C     INITIALLY BIAS TOWARD LONG PERIOD (AND FEWER STATIONS)
C     RETURN IF AT LEAST MINIMUM STATIONS

              IF(IK .GE. minsta) GO TO 140

C     FIND LOWEST CORRELATED STATION

              LEDIT = 1
              MINBAN = BANF(ISRT(1))

              DO 100 IT = 2,IK
                IF(BANF(ISRT(IT)) .LE. MINBAN) THEN
                  LEDIT = IT
                  MINBAN = BANF(ISRT(IT))
                END IF
  100         CONTINUE

C     REMOVE LOWEST CORRELATED STATION FROM THE ALGORITHM

              NAY(ISRT(LEDIT)) = 0
              NBY(ISRT(LEDIT)) = 0
            END IF

C     START SUBROUTINE OVER FROM THE BEGINNING

            GO TO 10
          END IF

C     SAVE OPTIMAL WINDOW VALUES

          LOPTA = IOPTA
          LOPTB = IOPTB

C     CALCULATE CONFIDENCE INTERVAL

          CALL SRTDVR(BAN,BN)
          CALL CALCZA(BAN,BN,WCCF,WCI)

C     STORE OPTIMIZED VALUES
          IF(WCI .LE. OPCI .OR. IOPSTA .LE. minsta) THEN

C     DEBUG PRINTOUT

            IF(DEBUG .GE. 2) WRITE(6,110) WCI,WCCF,IOPTB,IOPTA,
     *                                    (ISRT(IT),CONF(ISRT(IT)),
     *                                     OFF(ISRT(IT)),IT=1,IK)
  110                        FORMAT('WCI: ',F6.2,' OFF: ',F7.3,
     *                              ' IOPTB: ',I2,' IOPTA: ',I2,
     *                              ' (STA,CONF,OFF): ',
     *                              3(/,7(I3,2(',',F7.3),': ')))
            IOPIK = IK
            IOPSTA = IK
            JY1OP = JY1
            JY2OP = JY2
            OPBAN = NINT(BAN)
            OPCCF = WCCF
            OPCI = WCI

C     INITIALIZE VALUES USED TO COMPUTE OPTIMIZED WEIGHTED AVERAGE

            OPNET = 0.
            OPNUM = 0.

C     COMPUTE OPTIMIZED WEIGHTED AVERAGE

            DO 120 IT = 1,IK
              IS = ISRT(IT)
              IOPSRT(IT) = IS
              OPCONF(IS) = CONF(IS)
              OPNET = OPNET + DEL(IS,LY) / CONF(IS)
              OPNUM = OPNUM + 1. / CONF(IS)
  120       CONTINUE
          END IF
c         if not enough stations - STOP!
          if(ik .gt. minsta .and. ik .gt. iopik + 2) go to 140
  130   CONTINUE
      END IF

  140 RETURN

      END

C***********************************************************************
C End of Subroutine RSFILL.                                            *
C***********************************************************************


C***********************************************************************
C                          Subroutine SRTCNF                           *
C                                                                      *
C DATE OF LAST CHANGE:  22 April 1994                                  *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         04/22/94 *
C                                                                      *
C   1.0       Original Version                                01/02/86 *
C                                                                      *
C AUTHORS:  Claude Williams                                            *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  SRTCNF(integer ISRT(inmstn))                                 *
C                                                                      *
C DESCRIPTION:  Sorts the confidence array CONF from lowest to highest *
C               values.                                                *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Sorts the array ISRT, which points to the confidence array *
C           CONF; sorted from lowest to highest.                       *
C                                                                      *
C PARAMETERS:                                                          *
C                                                                      *
C   inmstn     Network of Stations Including Candidate                 *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
c  for common variables see filinit
c
C   I          Loop Counter                                            *
C                                                                      *
C   J          Loop Counter                                            *
C                                                                      *
C   K          Loop Counter                                            *
C                                                                      *
C   ISRT       Array of Index Numbers for CONF, Which Sorts CONF from  *
C   (inmstn)   Lowest to Highest                                       *
C                                                                      *
C***********************************************************************

      SUBROUTINE SRTCNF(ISRT)

      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.comm.incl'  

      DIMENSION ISRT(inmstn)

      COMMON/BLK7/CONF

      DIMENSION CONF(inmstn)

C     SORT LOWEST TO HIGHEST

      ISRT(1) = 1

      DO 30 K = 2,inmstn
        DO 20 I = 1,K-1
          IF(CONF(K) .LT. CONF(ISRT(I))) THEN
            DO 10 J = K-1,I,-1
              ISRT(J+1) = ISRT(J)
   10       CONTINUE

            ISRT(I) = K
            GO TO 30
          END IF
   20   CONTINUE

        ISRT(K) = K
   30 CONTINUE

      RETURN

      END

C***********************************************************************
C End of Subroutine SRTCNF.                                            *
C***********************************************************************


C***********************************************************************
C                          Subroutine SRTDVR                           *
C                                                                      *
C DATE OF LAST CHANGE:  16 August 1994                                 *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         09/16/94 *
C                                                                      *
C   1.0       Original Version                                01/02/86 *
C                                                                      *
C AUTHORS:  Claude Williams                                            *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  SRTDVR(real BAN, real BN)                                    *
C                                                                      *
C DESCRIPTION:  This subroutine calculates the ratios of the           *
C               candidate's precipitation values to those of its       *
C               neighbors for the window time of period and sorts them *
C               from lowest to highest.                                *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Returns the ratios of the candidate's precipitation values *
C           to those of its neighbors for the window time of period,   *
C           sorted from lowest to highest.                             *
C                                                                      *
C PARAMETERS:                                                          *
C                                                                      *
C   numrat     Number of Precipitation Ratios of Before Window Values  *
C              to After Window Values that Array XSDIF Will Hold       *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
c  for common variables see filinit
c
C   BAN        BN * BN                                                 *
C                                                                      *
C   BN         Count of Non-Missing Candidate's and Neighbors' Monthly *
C              Values that are Not Flagged as "X"/"U" in Window        *
C                                                                      *
C   ISD        Counter for Number of Ratios Calculated                 *
C                                                                      *
C   IW1        Loop Counter                                            *
C                                                                      *
C   IW2        Loop Counter                                            *
C                                                                      *
C   MIN        Index of Minimum Precipitation Ratio                    *
C                                                                      *
C   XSDTMP     Holds Precipitation Ratio Being Swapped Out             *
C                                                                      *
C***********************************************************************

      SUBROUTINE SRTDVR(BAN,BN)

      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.comm.incl'  

      COMMON/BLK4/XJY,XKY
      COMMON/BLK6/XSDIF

      DIMENSION XJY(20),XKY(20)
c      DIMENSION XSDIF(numrat), tsdif(numrat)
      DIMENSION XSDIF(numrat)

C     INITIALIZE VARIABLES

      ibn = bn
      iBAN = iBN * iBN
      BAN = iBAN
      ISD = 0

      DO 10 IW1 = 1,iBN
        DO 10 IW2 = 1,iBN
          ISD = ISD + 1
          XSDIF(ISD) = XKY(IW2) / XJY(IW1)
c          tsdif(isd) = xsdif(isd)
   10 CONTINUE

C     SORT RATIOS FROM LOWEST TO HIGHEST
      call hsort(xsdif, iban)

      RETURN

      END

C***********************************************************************
C End of Subroutine SRTDVR.                                            *
C***********************************************************************


C***********************************************************************
C                           Subroutine hsort                           *
C                                                                      *
c   Heapsort (see Numerical Recipes (p.231)- C. Williams 19Jan2001
C                                                                      *
C***********************************************************************

      SUBROUTINE HSORT(ARR,NVAL)

      dimension ARR(NVAL)
      
      l = NVAL/2 + 1
      ir = NVAL
      
   10 continue
      if(l .gt. 1) then      
        l = l - 1
        rARR = ARR(l)
      else
        rARR = ARR(ir)
        ARR(ir) = ARR(1)
        ir = ir - 1
        if(ir .eq. 1) then
          ARR(1) = rARR
          return
        endif
      endif
      i = l
      j = l+l
   20 if(j .le. ir) then
        if(j .lt. ir) then
          if(ARR(j) .lt. ARR(j+1)) j = j + 1
        endif
        if(rARR .lt. ARR(j)) then
          ARR(i) = ARR(j)
          i = j
          j = j + j
        else
          j = ir + 1
        endif
        go to 20
      endif
      
      ARR(i) = rARR
      go to 10
      end    
 
C***********************************************************************
C End of Subroutine SORTN.                                             *
C***********************************************************************


C***********************************************************************
C                           Subroutine MRGBN                           *
C                                                                      *
C DATE OF LAST CHANGE:  09 August 1994                                 *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         08/09/94 *
C                                                                      *
C   1.0       Original Version                                01/02/86 *
C                                                                      *
C AUTHORS:  Claude Williams                                            *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  MRGBN(double precision C2I, real CI, real CN,                *
C               real DEL(inmstn+1,inumyr), integer IK,                 *
C               integer ISRT(inmstn), integer MY1, integer MY2,        *
C               integer LY)                                            *
C                                                                      *
C DESCRIPTION:  This subroutine calculates the sum of differences and  *
C               differences squared between the candidate's data and   *
C               the weighted data of the candidate's neighbors for     *
C               temperature, and stores the candidate's and            *
C               candidate's neighbors' data for precipitation.         *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Returns sum of differences and differences squared between *
C           the candidate's data and the weighted data of the          *
C           candidate's neighbors for temperature, and stores the      *
C           candidate's and candidate's neighbors' data for            *
C           precipitation.                                             *
C                                                                      *
C PARAMETERS:                                                          *
C                                                                      *
C   inmstn     Network of Stations Including Candidate                 *
C                                                                      *
C   inumyr     Last Cell of Arrays that Hold Year Data; Currently Set  *
C              to Handle Through the Year 2025                         *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
c  for common variables see filinit
c
C   C2I        Sum of Temperature Differences Squared Between          *
C              Candidate and Neighbors Over Years in Window            *
C                                                                      *
C   CI         Sum of Temperature Differences Between Candidate and    *
C              Neighbors Over Years in Window                          *
C                                                                      *
C   CN         Count of Non-Missing Candidate's and Neighbors' Monthly *
C              Values that are Not Flagged as "X"/"U" in Window        *
C                                                                      *
C   DEL        Monthly Data                          -99.99 => Missing *
C   (inmstn+1,inumyr)  inmstn+1 Stations, inumyr Years                 *
C                                                                      *
C   IK         Index of Candidate's Neighbors Up to Which Calculations *
C              are Made                                                *
C                                                                      *
C   ISRT       Array of Index Numbers for CONF, Which Sorts CONF from  *
C   (inmstn)   Lowest to Highest                                       *
C                                                                      *
C   IST        Index of Station with Next Lowest Confidence Interval   *
C                                                                      *
C   IT         Loop Counter                                            *
C                                                                      *
C   IW         Loop Counter                                            *
C                                                                      *
C   LY         Current Year Being Processed                            *
C                                                                      *
C   MY1        Beginning of Window                                     *
C                                                                      *
C   MY2        End of Window                                           *
C                                                                      *
C   SD         Sum of Ratios of a Neighbor's Monthly Value to the      *
C              Confidence Interval                                     *
C                                                                      *
C   SW         Sum of Reciprocals of Confidence Intervals              *
C                                                                      *
C   X          Intermediate Calculation                                *
C                                                                      *
C***********************************************************************

      SUBROUTINE MRGBN(C2I,CI,CN,DEL,IK,ISRT,MY1,MY2,LY)

      INCLUDE 'posthomog.parm.incl'
      INCLUDE 'posthomog.fill.incl'  
      INCLUDE 'posthomog.comm.incl'

      DIMENSION ISRT(inmstn)
      DIMENSION DEL(inmstn+1,inumyr)
      DOUBLE PRECISION C2I

      COMMON/BLK4/XJY,XKY
      COMMON/BLK7/CONF

      DIMENSION CONF(inmstn),XJY(20),XKY(20)

C     INITIALIZE VARIABLES
      C2I = 0.
      CI = 0.
      iCN = 0

C     YEARS BEFORE CURRENT YEAR
      DO 20 IW = LY-1,MY1,-1

C       INITIALIZE SUM VARIABLES
        SD = 0.
        SW = 0.

C       CHECK FOR CANDIDATE'S MONTHLY DATA VALUE NOT MISSING NOR FLAGGED
C       AS "X"/"U"
        IF(DEL(1,IW) .gt. -99.) THEN

C         COMPUTE INTERMEDIATE CALCULATIONS FOR WEIGHTED VALUES
          DO 10 IT = 1,IK
            IST = ISRT(IT)
            IF(DEL(IST,IW) .le. -99.) GO TO 20
            SD = SD + DEL(IST,IW) / CONF(IST)
            SW = SW + 1. / CONF(IST)
   10     CONTINUE

C         DO CALCULATIONS SINCE STATIONS FOUND FROM ABOVE LOOP
          IF(SW .GT. 0.) THEN

C           COUNT NUMBER OF YEARS STATIONS FOUND
            iCN = iCN + 1

C           RETURN IF COUNT MORE THAN 20
            IF(iCN .GT. 20) THEN
              iCN = 20
              GO TO 50
            END IF

C           CALCULATE SUM OF DIFFERENCES AND DIFFERENCES SQUARED FOR
C           TEMPERATURE
            IF(INEL .LT. 4) THEN
              X = DEL(1,IW) - (SD / SW)
              CI = CI + X
              C2I = C2I + X * X

C           STORE CANDIDATE'S AND CANDIDATE'S NEIGHBORS' DATA FOR
C           PRECIPITATION
            ELSE
              XJY(iCN) = SD / SW
              XKY(iCN) = DEL(1,IW)
              if(xky(icn) .eq. 0.0) xky(icn) = 0.0001
              if(xjy(icn) .eq. 0.0) xjy(icn) = 0.0001
            END IF
          END IF
        END IF
   20 CONTINUE

C     YEARS AFTER CURRENT YEAR
      DO 40 IW = LY+1,MY2

C       INITIALIZE SUM VARIABLES
        SD = 0.
        SW = 0.

C       CHECK FOR CANDIDATE'S MONTHLY DATA VALUE NOT MISSING NOR FLAGGED
C       AS "X"/"U"
        IF(DEL(1,IW) .gt. -99.) THEN

C         COMPUTE INTERMEDIATE CALCULATIONS FOR WEIGHTED VALUES
          DO 30 IT = 1,IK
            IST = ISRT(IT)
            IF(DEL(IST,IW) .le. -99.) GO TO 40
            SD = SD + DEL(IST,IW) / CONF(IST)
            SW = SW + 1. / CONF(IST)
   30     CONTINUE

C         DO CALCULATIONS SINCE STATIONS FOUND FROM ABOVE LOOP
          IF(SW .GT. 0.) THEN

C           COUNT NUMBER OF YEARS STATIONS FOUND
            iCN = iCN + 1

C           RETURN IF COUNT MORE THAN 20
            IF(iCN .GT. 20) THEN
              iCN = 20
              GO TO 50
            END IF

C           CALCULATE SUM OF DIFFERENCES AND DIFFERENCES SQUARED FOR
C           TEMPERATURE
            IF(INEL .LT. 4) THEN
              X = DEL(1,IW) - (SD / SW)
              CI = CI + X
              C2I = C2I + X * X

C           STORE CANDIDATE'S AND CANDIDATE'S NEIGHBORS' DATA FOR
C           PRECIPITATION
            ELSE
              XJY(iCN) = SD / SW
              XKY(iCN) = DEL(1,IW)
              if(xky(icn) .eq. 0.0) xky(icn) = 0.0001
              if(xjy(icn) .eq. 0.0) xjy(icn) = 0.0001
            END IF
          END IF
        END IF
   40 CONTINUE

   50 cn = icn
      RETURN

      END

C***********************************************************************
C End of Subroutine MRGBN.                                             *
C***********************************************************************

C***********************************************************************
C                           Subroutine FILRS                           *
C                                                                      *
C DATE OF LAST CHANGE:  21 April 1994                                  *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         04/21/94 *
C                                                                      *
C   1.0       Original Version                                01/02/86 *
C                                                                      *
C AUTHORS:  Claude Williams                                            *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  FILRS(integer BANF(inmstn), real DEL(inmstn+1,inumyr),       *
C               integer IS, integer LY, integer INAY, integer INBY,    *
C               real OFF(inmstn), real WCI)                            *
C                                                                      *
C DESCRIPTION:  This subroutine calculates the number of nonmissing    *
C               monthly values in the window, the confidence           *
C               correction factors, and the confidence intervals.      *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Calculates number of nonmissing monthly values in the      *
C           window, confidence correction factors, and confidence      *
C           intervals.                                                 *
C                                                                      *
C PARAMETERS:                                                          *
C                                                                      *
C   inmstn     Network of Stations Including Candidate                 *
C                                                                      *
C   inumyr     Last Cell of Arrays that Hold Year Data; Currently Set  *
C              to Handle Through the Year 2025                         *
C                                                                      *
C   minyr      Minimum Number of Years in Window for Homogeneity       *
C                                                                      *
C SUBROUTINES:  CALCZA                                                 *
C               SRTDVR                                                 *
C               STNBN                                                  *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
c  for common variables see filinit
c
C   B2I        Not Used by this Subroutine                             *
C                                                                      *
C   BAN        BN * BN                                                 *
C                                                                      *
C   BANF       Values of BAN for Each Network Station                  *
C   (inmstn)                                                           *
C                                                                      *
C   BI         Not Used by this Subroutine                             *
C                                                                      *
C   BN         Count of Non-Missing Candidate's and Neighbors' Monthly *
C              Values that are Not Flagged as "X"/"U" in Window        *
C                                                                      *
C   DEL        Monthly Data                          -99.99 => Missing *
C   (inmstn+1,inumyr)  inmstn+1 Stations, inumyr Years                 *
C                                                                      *
C   INAY       Number of Years in After Part of Window for Network     *
C              Station IS                                              *
C                                                                      *
C   INBY       Number of Years in Before Part of Window for Network    *
C              Station IS                                              *
C                                                                      *
C   IS         Index of Candidate's Neighboring Station                *
C                                                                      *
C   JY1        Beginning of Window                                     *
C                                                                      *
C   JY2        End of Window                                           *
C                                                                      *
C   LY         Current Year Being Processed                            *
C                                                                      *
C   OFF        Array of Confidence Correction Factors for a            *
C   (inmstn)   Candidate's Nearest Neighbors                           *
C                                                                      *
C   WCCF       Confidence Correction Factor                            *
C                                                                      *
C   WCI        Confidence Interval                                     *
C                                                                      *
C***********************************************************************

      SUBROUTINE FILRS(BANF,DEL,IS,LY,INAY,INBY,OFF,WCI)

      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.comm.incl'  

      DIMENSION OFF(inmstn)
      DIMENSION DEL(inmstn+1,inumyr)
      INTEGER BANF(inmstn)

      COMMON/BLK7/CONF

      DIMENSION CONF(inmstn)
      DOUBLE PRECISION B2I

C     CALCULATE BEGINNING AND END OF WINDOW
      JY1 = LY - INBY + 1
      JY2 = LY + INAY

      CALL STNBN(B2I,BI,BN,DEL,IS,JY1,JY2,LY)

C     SET VARIABLES IF NOT ENOUGH YEARS
      IF(NINT(BN) .LT. minyr) THEN
        CONF(IS) = 99.99
        OFF(IS) = 0.
        WCI = 99.99

C     ENOUGH YEARS:  CALCULATE AND STORE BN * BN, CONFIDENCE CORRECTION
C     FACTORS, AND CONFIDENCE INTERVALS
      ELSE
        CALL SRTDVR(BAN,BN)
        CALL CALCZA(BAN,BN,WCCF,WCI)

        BANF(IS) = NINT(BAN)
        if(WCI .eq. 0.0) then
          conf(is) = 0.01
        else  
          CONF(IS) = WCI
        endif  
        OFF(IS) = WCCF
      END IF

      RETURN

      END

C***********************************************************************
C End of Subroutine FILRS.                                             *
C***********************************************************************

C***********************************************************************
C                           Subroutine FILCF                           *
C                                                                      *
C DATE OF LAST CHANGE:  09 August 1994                                 *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         08/09/94 *
C                                                                      *
C   1.0       Original Version                                01/02/86 *
C                                                                      *
C AUTHORS:  Claude Williams                                            *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  FILCF(integer BANF(inmstn), real DEL(inmstn+1,inumyr),       *
C               integer IS, integer LY, integer INAY, integer INBY,    *
C               real OFF(inmstn), real WCI)                            *
C                                                                      *
C DESCRIPTION:  This subroutine calculates the number of nonmissing    *
C               monthly values in the window, the confidence           *
C               correction factors, and the confidence intervals.      *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Calculates number of nonmissing monthly values in the      *
C           window, confidence correction factors, and confidence      *
C           intervals.                                                 *
C                                                                      *
C PARAMETERS:                                                          *
C                                                                      *
C   inmstn     Network of Stations Including Candidate                 *
C                                                                      *
C   inumyr     Last Cell of Arrays that Hold Year Data; Currently Set  *
C              to Handle Through the Year 2025                         *
C                                                                      *
C   minyr      Minimum Number of Years in Window for Homogeneity       *
C                                                                      *
C SUBROUTINES:  CALCTA                                                 *
C               STNBN                                                  *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
c  for common variables see filinit
c
C   B2I        Sum of Temperature Differences Squared Between          *
C              Candidate and Neighbors Over Years in Window            *
C                                                                      *
C   BAN        Total Number of Non-Missing Seasonal Values in Window   *
C              Less One (Degrees of Freedom)                           *
C                                                                      *
C   BANF       Array of Indices for Determining Value of T-Alpha from  *
C   (inmstn)   the Confidence Array RVT67 for a Candidate's Nearest    *
C              Neighbors                                               *
C                                                                      *
C   BI         Sum of Temperature Differences Between Candidate and    *
C              Neighbors Over Years in Window                          *
C                                                                      *
C   BN         Count of Non-Missing Candidate's and Neighbors' Monthly *
C              Values that are Not Flagged as "X"/"U" in Window        *
C                                                                      *
C   DEL        Monthly Data                          -99.99 => Missing *
C   (inmstn+1,inumyr)  inmstn+1 Stations, inumyr Years                 *
C                                                                      *
C   INAY       Number of Years in After Part of Window for Network     *
C              Station IS                                              *
C                                                                      *
C   INBY       Number of Years in Before Part of Window for Network    *
C              Station IS                                              *
C                                                                      *
C   INDTA      Index for Determining Value of T-Alpha from the         *
C              Confidence Interval Array RVT67                         *
C                                                                      *
C   IS         Index of Candidate's Neighboring Station                *
C                                                                      *
C   JY1        Beginning of Window                                     *
C                                                                      *
C   JY2        End of Window                                           *
C                                                                      *
C   LY         Current Year Being Processed                            *
C                                                                      *
C   OFF        Array of Confidence Correction Factors for a            *
C   (inmstn)   Candidate's Nearest Neighbors                           *
C                                                                      *
C   TA         Critical Value of T-Alpha from Confidence Interval      *
C              Array Used in Calculating the Confidence Interval       *
C                                                                      *
C   WCI        Confidence Interval                                     *
C                                                                      *
C***********************************************************************

      SUBROUTINE FILCF(BANF,DEL,IS,LY,INAY,INBY,OFF,WCI)

      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.comm.incl'  

      DIMENSION OFF(inmstn)
      DIMENSION DEL(inmstn+1,inumyr)
      INTEGER BANF(inmstn)

      COMMON/BLK7/CONF

      DIMENSION CONF(inmstn)

      DOUBLE PRECISION B2I

C     CALCULATE BEGINNING AND END OF WINDOW

      JY1 = LY - INBY + 1
      JY2 = LY + INAY

      CALL STNBN(B2I,BI,BN,DEL,IS,JY1,JY2,LY)

C     CALCULATE "BAN"

      BAN = BN - 1.

C     SET VARIABLES IF NOT ENOUGH YEARS

      IF(NINT(BN) .LT. minyr) THEN

        BANF(IS) = NINT(BAN)
        CONF(IS) = 99.99
        OFF(IS) = 0.
        WCI = 99.99

C     ENOUGH YEARS:  CALCULATE AND STORE CONFIDENCE CORRECTION FACTORS
C     AND CONFIDENCE INTERVALS

      ELSE
        CALL CALCTA(BAN,INDTA,TA)
        BANF(IS) = INDTA
        OFF(IS) = BI / BN
        if(B2I * BN - BI * BI .le. 0.0) then
          wci = 0.01
        else  
          WCI = (TA / BN) * SQRT((B2I * BN - BI * BI) / (BN - 1.))
        endif
        conf(is) = wci  
      END IF

      RETURN

      END

C***********************************************************************
C End of Subroutine FILCF.                                             *
C***********************************************************************

C***********************************************************************
C                          Subroutine CALCTA                           *
C                                                                      *
C DATE OF LAST CHANGE:  08 August 1994                                 *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         08/08/94 *
C                                                                      *
C   1.0       Original                                                 *
C                                                                      *
C AUTHOR:  Claude Williams                                             *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  CALCTA(real BAN, integer INDTA, real TA)                     *
C                                                                      *
C DESCRIPTION:  Determines the critical value of T-Alpha used in       *
C               calculating the confidence interval.                   *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Returns the critical value of T-Alpha used in calculating  *
C           the confidence interval.                                   *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
C   BAN        Total Number of Non-Missing Seasonal Values in Before   *
C              Window Less One (Degrees of Freedom)                    *
C                                                                      *
C   I          Integer Value of BAN                                    *
C                                                                      *
C   INDTA      Index for Determining Value of T-Alpha from the         *
C              Confidence Interval Array RVT67                         *
C                                                                      *
C   RVT67      Array of T-Alpha Values Used in Calculating the         *
C   (33)       Confidence Interval                                     *
C                                                                      *
C   TA         Critical Value of T-Alpha from Confidence Interval      *
C              Array Used in Calculating the Confidence Interval       *
C                                                                      *
C***********************************************************************

      SUBROUTINE CALCTA(BAN,INDTA,TA)

      real      RVT67(33)/1.7923,1.3294,1.2030,1.1466,1.1149,1.0943,
     *                    1.0801,1.0698,1.0622,1.0557,1.0511,1.0465,
     *                    1.0426,1.0398,1.0383,1.0352,1.0331,1.0313,
     *                    1.0303,1.0285,1.0275,1.0257,1.0248,1.0238,
     *                    1.0228,1.0228,1.0218,1.0210,1.0200,1.0200,
     *                    1.0150,1.0110,1.0060/

C     USE NUMBER OF NON-MISSING SEASONAL VALUES TO DETERMINE CONFIDENCE
C     INTERVAL ARRAY INDEX

      I = NINT(BAN)

      IF(I .LE. 30) THEN
        INDTA = I
      ELSE IF(I .LE. 35) THEN
        INDTA = 30
      ELSE IF(I .LE. 50) THEN
        INDTA = 31
      ELSE IF(I .LE. 90) THEN
        INDTA = 32
      ELSE
        INDTA = 33
      END IF

C     DETERMINE CONFIDENCE INTERVAL VALUE

      TA = RVT67(INDTA)

      RETURN

      END

C***********************************************************************
C End of Subroutine CALCTA.                                            *
C***********************************************************************


C***********************************************************************
C                          Subroutine CALCZA                           *
C                                                                      *
C DATE OF LAST CHANGE:  19 September 1994                              *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         09/19/94 *
C                                                                      *
C   1.0       Original                                                 *
C                                                                      *
C AUTHOR:     Claude Williams                                          *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  CALCZA(real BAN, real BN, real WCCF, real WCI)               *
C                                                                      *
C DESCRIPTION:  This subroutine calculates the confidence correction   *
C               factor and confidence interval for the precipitation   *
C               values.                                                *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Returns confidence correction factor and confidence        *
C           interval for precipitation values.                         *
C                                                                      *
C PARAMETER:                                                           *
C                                                                      *
C   numrat     Number of Precipitation Ratios of Before Window Values  *
C              to After Window Values that Array XSDIF Will Hold       *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
C   BAN        BN * BN                                                 *
C                                                                      *
C   BN         Count of Non-Missing Candidate's and Neighbors' Monthly *
C              Values that are Not Flagged as "X"/"U" in Window        *
C                                                                      *
C   C          Calculation for Finding the Upper and Lower Boundaries  *
C              for Calculating the Confidence Correction Factor and    *
C              the Confidence Interval                                 *
C                                                                      *
C   ICONF      Confidence Factor of the Network                        *
C                                                                      *
C   IRH        Upper Boundary Used in Calculating the Confidence       *
C              Correction Factor and the Confidence Interval           *
C                                                                      *
C   IRL        Lower Boundary Used in Calculating the Confidence       *
C              Correction Factor and the Confidence Interval           *
C                                                                      *
C   TA         Value Used in Calculating C; Dependent on ICONF         *
C                                                                      *
C   WCCF       Confidence Correction Factor                            *
C                                                                      *
C   WCH        Natural Log of XSDIF(IRH)                               *
C                                                                      *
C   WCI        Confidence Interval                                     *
C                                                                      *
C   WCL        Natural Log of XSDIF(IRL)                               *
C                                                                      *
C   XSDIF      Array of Precipitation Ratios                           *
C   (numrat)                                                           *
C                                                                      *
C***********************************************************************

      SUBROUTINE CALCZA(BAN,BN,WCCF,WCI)

      INCLUDE 'posthomog.parm.incl'  
      INCLUDE 'posthomog.comm.incl'  

      COMMON/BLK5/ICONF
      COMMON/BLK6/XSDIF

      DIMENSION XSDIF(numrat)

C     VALUE OF "TA" BASED ON CONFIDENCE LEVEL

      IF(ICONF .EQ. 2) THEN
        TA = 1.95
      ELSE IF(ICONF .EQ. 3) THEN
        TA = 2.57
      ELSE IF(ICONF .EQ. 4) THEN
        TA = 3.30
      ELSE IF(ICONF .EQ. 5) THEN
        TA = 4.05
      ELSE
        TA = 1.00
      END IF

C     CALCULATE UPPER AND LOWER VALUES USED IN DETERMINING THE
C     CONFIDENCE CORRECTION FACTOR AND CONFIDENCE INTERVAL

      C = (BAN / 2.) + TA * SQRT(BAN * (BN + (BAN / BN) + 1.) / 12.)
      IRH = NINT(C)
      IRL = NINT(BAN - C + 1.)

      IF(XSDIF(IRH) .EQ. 0. .OR. XSDIF(IRL) .EQ. 0.) THEN
        WRITE(6,10) C,IRH,IRL,XSDIF(IRH),XSDIF(IRL)
   10   FORMAT('XSDIF = ZERO!',/,'C = ',F10.5,' IRH = ',I6,' IRL = ',I6,
     *         ' XSDIF(IRH) = ',F10.5,' XSDIF(IRL) = ',F10.5)

C     TAKE NATURAL LOGS OF RATIOS AND CALCULATE CONFIDENCE CORRECTION
C     FACTOR AND CONFIDENCE INTERVAL

      ELSE
   20   WCH = ALOG(XSDIF(IRH))
        WCL = ALOG(XSDIF(IRL))
        WCCF = (WCH + WCL) / 2.
        WCI = (WCH - WCL) / 2.

C     MAKE INTERVAL WIDER IF (WCH = WCL) => WCI = 0

        IF(ABS(WCI) .LT. 0.00003) THEN
          IRH = IRH + 1
          IRL = IRL - 1

C     IF ALL VALUES ARE THE SAME, INITIALIZE CONFIDENCE INTERVAL

          IF(IRH .GT. BAN .OR. IRL .EQ. 0) THEN
            WCI = 0.001

C     RECALCULATE CONFIDENCE CORRECTION FACTOR AND CONFIDENCE INTERVAL

          ELSE
            GO TO 20
          END IF
        END IF 
      END IF

      RETURN

      END

C***********************************************************************
C End of Subroutine CALCZA.                                            *
C***********************************************************************

C***********************************************************************
C                           Subroutine STNBN                           *
C                                                                      *
C DATE OF LAST CHANGE:  09 August 1994                                 *
C                                                                      *
C MODIFIER:  David Bowman                                              *
C                                                                      *
C VERSION     MODIFICATION                                      DATE   *
C -------     ------------                                    -------- *
C   1.1       SUN Workstation Version                         08/09/94 *
C                                                                      *
C   1.0       Original Version                                01/02/86 *
C                                                                      *
C AUTHORS:  Claude Williams                                            *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  STNBN(double precision C2I, real CI, real CN,                *
C               real DEL(inmstn+1,inumyr), integer IS, integer MY1,    *
C               integer MY2, integer LY)                               *
C                                                                      *
C DESCRIPTION:  Calculates (for temperature) the sum of the            *
C               differences and the sum of the differences squared     *
C               between the candidate and its neighbors.  Stores (for  *
C               precipitation) the candidate's and neighbors'          *
C               precipitation values.                                  *
C                                                                      *
C NOTES:  None.                                                        *
C                                                                      *
C RESULTS:  Returns (for temperature) the sum of the differences and   *
C           the sum of the differences squared between the candidate   *
C           and its neighbors; and returns (for precipitation) the     *
C           candidate's and neighbors' precipitation values.           *
C                                                                      *
C PARAMETERS:                                                          *
C                                                                      *
C   inmstn     Network of Stations Including Candidate                 *
C                                                                      *
C   inumyr     Last Cell of Arrays that Hold Year Data; Currently Set  *
C              to Handle Through the Year 2025                         *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
c  for common variables see filinit
c
C   C2I        Sum of Temperature Differences Squared Between          *
C              Candidate and Neighbors Over Years in Window            *
C                                                                      *
C   CI         Sum of Temperature Differences Between Candidate and    *
C              Neighbors Over Years in Window                          *
C                                                                      *
C   CN         Count of Non-Missing Candidate's and Neighbors' Monthly *
C              Values that are Not Flagged as "X"/"U" in Window        *
C                                                                      *
C   DEL        Monthly Data                          -99.99 => Missing *
C   (inmstn+1,inumyr)  inmstn+1 Stations, inumyr Years                 *
C                                                                      *
C   IS         Station Index Number                                    *
C                                                                      *
C   IW         Loop Counter                                            *
C                                                                      *
C   LY         Current Year Being Processed                            *
C                                                                      *
C   MY1        Beginning of Window                                     *
C                                                                      *
C   MY2        End of Window                                           *
C                                                                      *
C   X          Temperature Difference Between Candidate and Neighbor   *
C                                                                      *
C***********************************************************************

      SUBROUTINE STNBN(C2I,CI,CN,DEL,IS,MY1,MY2,LY)

      INCLUDE 'posthomog.parm.incl' 
      INCLUDE 'posthomog.fill.incl'  
      INCLUDE 'posthomog.comm.incl'

      DIMENSION DEL(inmstn+1,inumyr)
      DOUBLE PRECISION C2I

      COMMON/BLK4/XJY,XKY

      DIMENSION XJY(20),XKY(20)

C     INITIALIZE VARIABLES
      C2I = 0.
      CI = 0.
      iCN = 0

C     CALCULATE TEMPERATURE DIFFERENCES OR STORE PRECIPITATION VALUES
C     FOR YEARS IN WINDOW
      DO 10 IW = LY-1,MY1,-1

C       IF NEITHER CANDIDATE'S NOR NEIGHBOR'S MONTHLY VALUE IS MISSING NOR
C       AS "X"/"U", COUNT AND CALCULATE DIFFERENCES OR STORE PRECIPITATION
C       RATIOS
        IF(DEL(1,IW) .gt. -99. .AND. DEL(IS,IW) .gt. -99.) THEN

C         COUNT NUMBER OF DATA VALUES
          iCN = iCN + 1

C         SET COUNT TO 20 AND RETURN IF COUNT IS MORE THAN 20
          IF(iCN .GT. 20) THEN
            iCN = 20
            GO TO 30
          END IF

C         CALCULATE AND SUM DIFFERENCES AND DIFFERENCES SQUARED FOR
C         TEMPERATURE
          IF(INEL .LT. 4) THEN
            X = DEL(1,IW) - DEL(IS,IW)
            CI = CI + X
            C2I = C2I + X * X

C         STORE CANDIDATE'S AND NEIGHBORS' PRECIPITATION VALUES
          ELSE
            XJY(iCN) = DEL(IS,IW)
            XKY(iCN) = DEL(1,IW)
              if(xky(icn) .eq. 0.0) xky(icn) = 0.0001
              if(xjy(icn) .eq. 0.0) xjy(icn) = 0.0001
          END IF
        END IF
   10 CONTINUE

C     CALCULATE TEMPERATURE DIFFERENCES OR STORE PRECIPITATION VALUES
C     FOR YEARS IN WINDOW
      DO 20 IW = LY+1,MY2

C       IF NEITHER CANDIDATE'S NOR NEIGHBOR'S MONTHLY VALUE IS MISSING NOR
C       AS "X"/"U", COUNT AND CALCULATE DIFFERENCES OR STORE PRECIPITATION
C       RATIOS
        IF(DEL(1,IW) .gt. -99. .AND. DEL(IS,IW) .gt. -99.) THEN

C         COUNT NUMBER OF DATA VALUES
          iCN = iCN + 1

C         SET COUNT TO 20 AND RETURN IF COUNT IS MORE THAN 20
          IF(iCN .GT. 20) THEN
            iCN = 20
            GO TO 30
          END IF

C         CALCULATE AND SUM DIFFERENCES AND DIFFERENCES SQUARED FOR
C         TEMPERATURE
          IF(INEL .LT. 4) THEN
            X = DEL(1,IW) - DEL(IS,IW)
            CI = CI + X
            C2I = C2I + X * X

C         STORE CANDIDATE'S AND NEIGHBORS' PRECIPITATION VALUES
          ELSE
            XJY(iCN) = DEL(IS,IW)
            XKY(iCN) = DEL(1,IW)
              if(xky(icn) .eq. 0.0) xky(icn) = 0.0001
              if(xjy(icn) .eq. 0.0) xjy(icn) = 0.0001
          END IF
        END IF
   20 CONTINUE

   30 cn = icn
      RETURN

C***********************************************************************
C End of Subroutine STNBN.                                             *
C***********************************************************************
      END
