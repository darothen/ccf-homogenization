      subroutine read_hist(mmunit,mounit,ntstn,nht,mmdates)
      
c     Read the history data for the current network from either 
c       the legacy USHCN format (irandom == 0)
c       or Matt's meta.txt file format (irandom != 0). 
c     Store these data into the UCP's Station History "Hit" array (sahist)
c
c   NOTE: CDMP history is the COOP data files has been moved from source = 3
c              to source = 1 and resorted.                      08aug05
c
c
c   Subroutine Parameters
c     mmunit - SHAP formatted input file unit
c     mounit - Output of compressed meta record's (MM's metadata format)
c     ntstn - Candidate stations sub-network
c     nht - Number of suspect changepoints
c
      INCLUDE 'inhomog.parm.mthly.incl'  
      include 'inhomog.comm.mthly.incl'

c     setup the location and amount arrays
      real amt(ninh),slope(ninh)
      integer loc(ninh)

c     network station list (base station first)
      character*6 ntstn(maxstns), cin

c     move info to convey from MM's metadata to the alignment routine
      integer move(ninh), mday(ninh)
      
c     mmdates contains the dates of going to/from MMTS
      integer mmdates(maxstns,2)

      lin = 0
      
c     Read Matt's meta.txt
c       since, in the realworld, there are no amounts - they are ignored...
      if(irandom .ne. 0) then
      
c       go thru each station
        do istn = 1,maxstns
        
          do ichg = 1, ninh
            move(ichg) = 0
            amt(ichg) = amiss
            slope(ichg) = 0
            mday(ichg) = 31
          enddo  
          mday(1) = 1

c         setup the first move to begin at the beginning
          move(1) = 1

c         read the station history
          read(mmunit,*,end=999) in,is,nc,(move(ic),amt(ic),ic=2,nc+1)
          
          write(cin,'(i6)') in
          
          move(nc+2) = 1200

c         align SHF moves with missing and minlen data segments
          call alignmoves(cin, istn, mounit, nc+2, move, amt, 
     *      mday)
        enddo

      else
c       get the legacy USHCN history data
        call inshp4(ntstn, inread, amt, mmdates, mounit)  
        if(inread .eq. 0) stop
      endif
      return
      
c     read error in Matt's meta.txt
  999 call perror(' Unable to read Matt meta.txt: ' // mattmeta)
      stop
      end

C***********************************************************************
C                          Subroutine INSHP4                           *
C                                                                      *
C                                                                      *
C VERSION     MODIFICATION                                     DATE    *
C -------     ------------                                  ---------- *
c
c INSHP4  incorporated into the UCPMONTHLY program 
c             NOTE: removes data arrays, HISTORY ONLY!!!!!
c                                                           19 jan 2005
c
C   1.1       SUN Workstation Version                       10/24/1995 *
C                                                                      *
C   1.0       Original Version                                         *
C                                                                      *
C AUTHOR:  Claude Williams                                             *
C                                                                      *
C MODIFIERS:  Pam Hughes                                               *
C             David Bowman                                             *
C                                                                      *
C USAGE:  INSHP4(character*6 ntstns, integer INEL, 
c                   integer sahist(ninh,nstns), integer inread)
C                                                                      *
C DESCRIPTION:  This subroutine is passed an array containing stations *
C               on which to get data, an array of the station          *
C               identification numbers of the candidate station and    *
C               its 20 most correlated closest neighbors, begin year,  *
C               end year, and element.  The subroutine reads in the    *
C                history information for   *
C               the appropriate stations, and determines station moves *
C               and move dates (day of month).                         *
C                                                                      *
C REQUIREMENTS:                                                        *
C                                                                      *
C   Input Files:  hcn_mmts_max_data                                    *
C                 hcn_mmts_mean_data                                   *
C                 hcn_mmts_min_data                                    *
C                 hcn_mmts_pcp_data                                    *
C                 hcn_station_history                                  *
C                                                                      *
C NOTES:  Fixed bugs in decoding of observation time, in determining   *
C         move using distance and direction from previous location,    *
C         and in initializing of IHDAY array. - D. Bowman: 09/24/93    *
C                                                                      *
C         Added an initialization loop for candidate station's flags   *
C         after realization that the absence of this reinitialization  *
C         during recursive SHAP runs jumbles the O/P flags, and        *
C         retains candidate station's flags on adjusted record after   *
C         grab of the data value. - P. Hughes:  09/16/93               *
C                                                                      *
C         Changed back to 20 neighbors (from 40). - P. Hughes:         *
C         04/09/93                                                     *
C                                                                      *
C         C. Williams originally let loop through 22, but NTR was only *
C         declared as 21.  Finagled NTR in calling program and         *
C         subroutine to go with larger indexing. - P. Hughes:          *
C         12/17/92                                                     *
C                                                                      *
C RESULTS:  Stores the history move *
C           dates (day of month) for the array of stations.     *
C                                                                      *
C PARAMETERS:                                                          *
C                                                                      *
C   nstns     Network of Stations Including Candidate                 *
C                                                                      *
C   numyr     Last Cell of Arrays that Hold Year Data; Currently Set  *
C              to Handle Through the Year 2025                         *
C                                                                      *
C FUNCTIONS:                                                           *
C                                                                      *
C   DIRTOI     Converts Character Direction to Numeric Value           *
C                                                                      *
C VARIABLES:                                                           *
C                                                                      *
C   BEGD       Begin Day from Station History Record                   *
C                                                                      *
C   BEGM       Begin Month from Station History Record                 *
C                                                                      *
C   BEGY       Begin Year from Station History Record                  *
C                                                                      *
C   BM         Begin Month                                             *
C                                                                      *
C   DREC       Number of Data Record to Read                           *
C                                                                      *
C   DYEAR      Year Read from Data Record                              *
C                                                                      *
C   ELEV       Elevation from Station History Record                   *
C                                                                      *
C   EM         End Month                                               *
C                                                                      *
C   ENDD       End Day from Station History Record                     *
C                                                                      *
C   ENDM       End Month from Station History Record                   *
C                                                                      *
C   ENDY       End Year from Station History Record                    *
C                                                                      *
C   FIRSTY     First Year of Data if Later Year Than begyr               *
C                                                                      *
C   HISTCD     History Code              0 => No Move                  *
C                                        1 => Real Move                *
C                                        2 => Observer Move Only       *
C                                                                      *
C   HREC       Number of History Record to Read                        *
C                                                                      *
C   I          Loop Counter                                            *
C                                                                      *
C   endyr        Latest Year to Run Program                              *
C                                                                      *
C   IDIRDF     Absolute Value of Difference in Directions to Determine *
C              if Post Office Direction Changed by 90 Degrees or More  *
C                                                                      *
C   IHDAY      Move Number Array; Holds Begin Day of Move When Begin   *
C   (30)       Day is Not = "1st"                                      *
C                                                                      *
C   begyr        Earliest Year to Begin Program                          *
C                                                                      *
C   INEL       Element                        1 => Maximum Temperature *
C                                             2 => Minimum Temperature *
C                                             3 => Mean Temperature    *
C                                             4 => Precipitation       *
C                                                                      *
C   INSHT      Instrument Height from Station History Record:          *
C                          (1:2) Holds Precipitation Instrument Height *
C                          (3:4) Holds Temperature Instrument Height   *
C                                                                      *
C   INSTR      Instruments from Station History Record                 *
C                                                                      *
C   LAT        Latitude from Station History Record                    *
C                                                                      *
C   LDIR       Direction from Previous Location from Station History   *
C              Record                                                  *
C                                                                      *
C   LDIS       Distance from Previous Location from Station History    *
C              Record                                                  *
C                                                                      *
C   LELEV      Previous Record's ELEV                                  *
C                                                                      *
C   LINSHT     Previous Record's INSHT                                 *
C                                                                      *
C   LINSTR     Previous Record's INSTR                                 *
C                                                                      *
C   LLAT       Previous Record's LAT                                   *
C                                                                      *
C   LLON       Previous Record's LON                                   *
C                                                                      *
C   LOBSER     Previous Record's OBSER                                 *
C                                                                      *
C   LON        Longitude from Station History Record                   *
C                                                                      *
C   LPODIB     Previous Record's PODIB                                 *
C                                                                      *
C   LPODIS     Previous Record's PODIS                                 *
C                                                                      *
C   LTMPOB     Previous Record's TMPOB                                 *
C                                                                      *
C   eyear     Minimum of ENDY and endyr                                 *
C                                                                      *
C   MTHIND     Index for Month Component of Arrays                     *
C                                                                      *
C   byear     Maximum of BEGY and begyr                                 *
C                                                                      *
C   MOVNUM     Number of Station Moves                                 *
C                                                                      *
C   NH         Net History                                             *
C   (nstns,numyr,12)  nstns Stations, numyr Years, 12 Months       *
C                                    0 => No Move                      *
C                                    1 => Real Move and Begin Day = 01 *
C                                    2 => Real Move and Begin Day > 01 *
C                                    9 => Missing                      *
C                                                                      *
C   NOTMOV     Determines if Station Moved:                            *
C              (1:1) = 0 => Change in Observation Time (Temperature)   *
C              (1:1) = 1 => No Change in Observation Time or Current/  *
C                           Previous Observation Time Unknown          *
C                           (Temperature)                              *
C              (1:1) = 1 => (Precipitation)                            *
C              (2:2) = 0 => Change in Instrument Height                *
C              (2:2) = 1 => No Change in Instrument Height, or         *
C                           Current/Previous Instrument Height Unknown *
C              (3:3) = 0 => Change in Instrument "CRS" or "NSS"        *
C                           (Temperature)                              *
C              (3:3) = 1 => No Change in Instrument "CRS" and "NSS"    *
C                           (Temperature)                              *
C              (3:3) = 0 => Change in Instrument "FP", "NSRG", "RRNG", *
C                           or "TB" (Precipitation)                    *
C              (3:3) = 1 => No Change in Instrument "FP" and "NSRG"    *
C                           and "RRNG" and "TB" (Precipitation)        *
C              (4:4) = 0 => Change in Distance and/or Direction from   *
C                           Previous Location                          *
C              (4:4) = 1 => No Change in Distance and Direction from   *
C                           Previous Location, or Unknown Distance     *
C                           from Previous Location                     *
C                                                                      *
C   OBSER      Observer from Station History Record                    *
C                                                                      *
C   OBTIM      Observation Time from Station History Record            *
C                                                                      *
C   iPASS       Flag to Check for Real Move:                            *
C                0 => First Record of Station History; Do Not Check    *
C                     for Move                                         *
C                1 => Other than First Record of Station History;      *
C                     Check for Move                                   *
C                                                                      *
C   PODIB      Direction (and Block Column) from Post Office from      *
C              Station History Record                                  *
C                                                                      *
C   PODIS      Distance from Post Office from Station History Record   *
C                                                                      *
C   POSCHG     Station Move                             0 => No Move   *
C                                                     > 0 => Real Move *
C                                                                      *
C   RECDAT     Monthly Data Value from Data Record                     *
C   (12)                                                               *
C                                                                      *
C   istn     Index for Station Component of Arrays                   *
C                                                                      *
C   TMPOB      Temperature Observation Time Decoded from OBTIM         *
C                                                                      *
C   YEAR       Index for Year to Look for Data/History                 *
C                                                                      *
C***********************************************************************

      SUBROUTINE INSHP4(ntstn, inread, amt, mmdates, mounit)

c      IMPLICIT INTEGER (A-Z)

      include 'inhomog.parm.mthly.incl'
      include 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'
      
c     instrument string parameters
      PARAMETER (ninstr = 37, maxnin = 11)

c     network station list (base station first)
      character*6 ntstn(maxstns)

      character*30 movstr
      DIMENSION IHDAY(30)
      
      integer mmdates(maxstns,2)

      CHARACTER*3 LDIR
      CHARACTER*4 INSHT,LPODIB,LTMPOB,NOTMOV,OBTIM,PODIB,TMPOB
      character*4 linsht
      character*2 insth, linsth/'  '/
      CHARACTER*6 LAT,LLAT, cstnstr
      CHARACTER*7 LLON,LON
      character*11 ldisdir
      CHARACTER*46 LOBSER,OBSER
      character*132 histfile, elemfile
      REAL RECDAT(12)
      integer src, move(ninh), mday(ninh)
      real amt(ninh)
      integer begm, begy, begd, endy, endm, endd, byear, eyear
      integer histcd, elev
c     because of the measurement changes in Lat/Lon there is a dramatic
c        change in lat/lon resolution as of 1998. In the new code, this 
c        will have to be dealt with better, probably with testing the
c        measurement type, then the amount of change (=45 arcseconds)
      real lleps /0.0125/
      real latdeg,latmin,latsec,londeg,lonmin,lonsec
      real alat,allat,alon,allon,aallat,aallon

      character*5 cstrng(ninstr) / 'AI   ','CRS  ','DT   ','EVA  ',
     *     'FP   ','HYTHG','MXMN ','NRIG ','NSRG ','NSS  ','RRIG ',
     *     'RRNG ','SDE  ','SG   ','SRG  ','SS   ','TG   ','DGT  ',
     *     'TB   ','EVO  ','MMTS ','TELSY','HYGRO','HY6  ','HY8  ',
     *     'SFP  ','SRRNG','SSG  ','SSRG ','STB  ','AMOS ','AUTOB',
     *     'PSY  ','ASOS ','PLAST','STO  ','HYGR '/
      character*5 cstnin(maxnin)
      integer INSTR(ninstr),LINSTR(ninstr)
c      character*5 prcp1/'     '/, tmpr1/'     '/, 
c     *  lprcp1/'     '/, ltmpr1/'     '/

      inread = 1
      ihdbug = 0
      
c     tolerance limit for lat/lon changes wrt GPS ONLY entry (45 arcseconds)
      tolgps = 45./60./60.
      
      maxmo = nmo
      if(itimeres .eq. 0) maxmo = numyr

C     INITIALIZE VARIABLES
      FIRSTY = 0
      DO 10 I = 1,30
        IHDAY(I) = 0
   10 CONTINUE

C     NOTE - HISTORY ARRAY INITIALIZED ELSEWHERE

C     GET STATION NETWORK HISTORY
      nohist = 0
      DO 150 istn = 1,numsubs
        cstnstr = ntstn(istn)
c       define paths to USHCN & CDMP-3220 datasets (and current type - icoop)
c       NOTE: ncand==0 assumes USHCN format for "COOP" stations
        if(istn .le. ncand .or. ncand .eq. 0) then
          histfile = incand(1:lnblnk(incand)) // 'his/' // 
     *      cstnstr // '.his'
          icoop = 0
        else
          histfile = incoop(1:lnblnk(incoop)) // 'his/' // 
     *      cstnstr // '.his'
          icoop = 1
        endif    

C       OPEN STATION HISTORY FILE
        OPEN(9,FILE=histfile,STATUS='old',ERR=50)
        goto 55
          
   50   nohist = nohist + 1
        print *, 'No history: ', histfile(1:lnblnk(histfile))
        goto 150
c       print *,' No history set (no file) for: ',istn,' ',cstnstr
c       print *, histfile
c       goto 150

C       INITIALIZE VARIABLES
   55   MOVNUM = 0
        iPASS = 0
c       last dates for the HCN source (=0)
        lys0 = 0
        lms0 = 0
        lds0 = 0

        do ichg = 1, ninh
          move(ichg) = 0
          amt(ichg) = amiss
          mday(ichg) = 0
        enddo  
 
C       STATION HISTORY FILE LOOP
C       READ STATION HISTORY RECORDS
        do while (1 .eq. 1)
c         initialize instrument string
          do iin = 1, ninstr
            linstr(iin) = instr(iin)
            instr(iin) = 0
          end do
              
   85     read(9,90,end = 120,err=220) src,begy,begm,begd,endy,endm,
     *       endd, latdeg, latmin, latsec,londeg, lonmin, lonsec,
     *       ldisdir,elev,insht,obtim,cstnin
   90     format(i1,7x,2(1x,i4,2i2),1x,f3.0,2(1x,f2.0),1x,f4.0,
     *           2(1x,f2.0),1x,a11,1x,i5,2x,a4,1x,a4,5x,13(a5,1x))
          histcd = 0

c         The daily time of observation history records are not used here...
          if(src .eq. 1) go to 85
c         CDMP data should not be used with HCN stations!
          if(src .eq. 3 .and. icoop .eq. 0) goto 85
            
c         MSHR history records are used after the last begin date of the
c           USHCN & CDMP history records
          if(src .eq. 2 .and. (begy .lt. lys0 .or. 
     *      (begy.eq.lys0 .and. begm.lt.lms0) .or.
     *      (begy.eq.lys0 .and. begm.eq.lms0 .and. begd.le.lds0))) 
     *      goto 85
     
c         shift the instrument strings for the CDMP_3220 history
c          (Coop data (coop = 1: src=0)
          if(src .eq. 0 .and. icoop .eq. 1) then
c           there are only 10 max instruments
            do inin = 1, 10
              if(lnblnk(cstnin(inin)) .eq. 1) goto 92
              do while (cstnin(inin)(1:1).ne.' ')
                cstnin(inin)(1:4) = cstnin(inin)(2:5)
                cstnin(inin)(5:5) = ' '
              enddo
            enddo
          endif      
     
c         Get lat/lon ready
   92     alat = latdeg + ((latmin + latsec/60.0)/60.0)
          alon = londeg + ((lonmin + lonsec/60.0)/60.0)

c         decode last dir/dir from the USHCN
          if(src .eq. 0 .and. icoop .eq. 0) then
c           see if there is a decimal in distance
            if(ldisdir(1:1) .eq. '.' .or. ldisdir(2:2) .eq. '.' .or.
     *        ldisdir(3:3) .eq. '.') then
              read(ldisdir(1:3), '(f3.1)') dtemp
              ldis = int(dtemp * 10.)
            else
              read(ldisdir(1:3),'(i3)') ldis
            endif  
            ldir = ldisdir(5:7)
          endif  

c         find instruments
          do inin= 1, maxnin
            do iin = 1, ninstr
              if(cstnin(inin) .eq. cstrng(iin)) then
                instr(iin) = 1
                go to 91
              endif
            end do    
   91       continue
          end do  

c         find instrument height
C         POSITIONS 3 AND 4 HOLD TEMPERATURE INSTRUMENT HEIGHT
C         POSITIONS 1 AND 2 HOLD PRECIPITATION INSTRUMENT HEIGHT
          IF(INEL .LT. 4) THEN
            insth = insht(3:4)
          ELSE
            insth = insht(1:2)
          END IF

C         GET HISTORY FOR CORRECT TIME FRAME
          IF(ENDY .GE. begyr .AND. BEGY .LE. endyr) THEN

C           DECODE TEMPERATURE OBSERVATION TIME:  OBSERVATION TIME = "xxHR" OR
C           "TRID" => COPY TO TEMPERATURE TIME; OBSERVATION TIME = "9xx9" =>
C           COPY "xx" TO TEMPERATURE TIME; OTHERWISE, COPY TO TEMPERATURE
C           TIME.
            IF(INEL .LT. 4) THEN
              IF(OBTIM(3:4) .EQ. 'HR' .OR. OBTIM(3:4) .EQ. 'ID') THEN
                TMPOB = OBTIM
              ELSE IF(OBTIM(1:1) .EQ. '9' .AND.
     *                OBTIM(2:2) .NE. '9') THEN
                TMPOB = OBTIM(2:3)
              ELSE
                TMPOB = OBTIM(3:4)
              END IF
            END IF

C           FILL IN MISSING DATES
            IF(BEGM .EQ. 99) BEGM = 6
            IF(BEGD .EQ. 99) BEGD = 15
            IF(ENDM .EQ. 99 .AND. ENDD .EQ. 99 .AND.
     *         ENDY .EQ. 9999) THEN
              ENDM = 12
              ENDD = 31
              ENDY = endyr
            ELSE
              IF(ENDM .EQ. 99) ENDM = 6
              IF(ENDD .EQ. 99) ENDD = 15
            END IF

C           MERGE INSTRUMENTS

C*      *****************************BACKGROUND*******************************
C
C             MODIFIED THE INSTRUMENT TYPES "EQUATES" FOR "HYGROS" AND
C             "SHIELDS". - P. Hughes 07/15/91
C
C             SET ALL "EQUATES" TO INSTR(2) ("CRS") - "REFRESHES" SAVED LIST 
C             IN THE EVENT AN INSTRUMENT CHANGE IN "EQUATES" OCCURS.
C             - P. Hughes 12/16/91
C
C             "MMTS" HAVE BEEN ADJUSTED A PRIORI
C             "HYGROS" WILL BE EVALUATED BY "SHAP"
C             CHANGES FROM "DT" TO "CRS", ETC. ARE EVALUATED
C             CHANGES FROM "CRS", ETC. TO "HYGRO", ETC. ARE EVALUATED
C             - P. Hughes 04/09/93
C
C             DECISION MADE TO TREAT "DT"s WITH/WITHOUT SHELTERS AS SAME
C             - P. Hughes 04/09/93
c
c             ----------------------- The 2000 update ------------------------
c             Am changing back some of the instrumentation before to 90's.....
c             1) Setting the MMTS apart again (#21)
c             Changes for the Millenneum are:
c             1) Added ASOS as a instrument (#34)
c             2) For Precip also added 4"Plastic (#35) and STO in Hiwaii (#36)
c             3) Source 2 (the MSHR) generic Hygrothemometer (assumed HY8) (#37)
c
c             Master Station History Records (src=2) have two new fields:
c             1) the primary precip instrument and
c             2) the primary temperature instrument
c             USHCN SHF (src=1) primary instruments are determined the from the mix
c               of instruments at any given time.
c                                                            09 May 2001 cw

c           MXMN == HYTHG == SS == TG == DGT
c           E2 - modification remove MMTS from "equate"
c                  i.e. add back MMTS only record as changepoint
            IF(INSTR(7) .EQ. 1 .OR.
     *         INSTR(6) .EQ. 1 .OR. INSTR(16) .EQ. 1 .OR.
     *         INSTR(17) .EQ. 1 .OR. INSTR(18) .EQ. 1) then
              INSTR(2) = 1
            ELSE
              INSTR(2) = 0
            END IF
              
c           HYGR == HY8
            If(instr(37) .eq. 1) instr(25) = 1

C           ADDED "EQUATES" FOR THE SHIELDED PRECIPITATION INSTRUMENTS
c           SFP == FP : SRRNG == RRNG : SSRG == SRG
            IF(INSTR(26) .EQ. 1) INSTR(5) = 1
            IF(INSTR(27) .EQ. 1) INSTR(12) = 1
            IF(INSTR(29) .EQ. 1) INSTR(15) = 1

C           DETERMINE BEGIN AND END YEARS FOR STATION MOVES LOOP
            if(begy .lt. begyr) then
              begm = 1
              begd = 1
              byear = begyr
            else
              byear = begy
            end if

            if(endy .gt. endyr) then
              endm = 12
              endd = 31
              eyear = endyr
            else
              eyear = endy
            end if

C           LOOP FOR DETERMINING STATION MOVES
            movstr = ''
            byear = byear

C           TEST FOR REAL MOVE
            IF(iPASS .EQ. 1) THEN

C             INITIALIZE VARIABLE
              NOTMOV = '0000'

C             CHECK OBSERVATION TIME
C             NO OBSERVATION TIME CHANGE IF:  PRECIPITATION; UNKNOWN
C             OBSERVATION TIME; UNKNOWN PREVIOUS OBSERVATION TIME; OR
C             OBSERVATION TIME IS SAME AS PREVIOUS OBSERVATION TIME.
              IF(INEL.EQ.4 .OR. TMPOB.EQ.'99' .OR. tmpob.eq.'  '
     *            .or. LTMPOB.EQ.'99' .OR. TMPOB.EQ.LTMPOB) then
                NOTMOV(1:1) = '1'
              else
                movstr = movstr(1:lnblnk(movstr)) // ' OBT'
              endif  

C             CHECK INSTRUMENT HEIGHT
C             NO INSTRUMENT HEIGHT CHANGE IF:  UNKNOWN INSTRUMENT HEIGHT;
C             UNKNOWN PREVIOUS INSTRUMENT HEIGHT; OR INSTRUMENT HEIGHT IS SAME
C             AS PREVIOUS INSTRUMENT HEIGHT.

c              print *,' insth,linsth: ', insth, ' ', linsth
              IF(insth.EQ.'99' .OR. insth.eq.'  ' .or.
     *           linsth.EQ.'99' .OR. linsth.eq.'  ' .or.
     *           insth .EQ. linsth) THEN
                NOTMOV(2:2) = '1'
              ELSE
                HISTCD = 1
                movstr = movstr(1:lnblnk(movstr))  // ' IHT'
              END IF

C             TEST INSTRUMENTS
C             NO INSTRUMENT CHANGE IF:  TEMPERATURE, AND NO CHANGE IN "CRS" AND
C             NO CHANGE IN "NSS"; OR PRECIPITATION, AND NO CHANGE IN "SRG" OR NO
C             CHANGE IN "FP" AND NO CHANGE IN "NSRG" AND NO CHANGE IN "RRNG" AND
C             NO CHANGE IN "TB".
              IF(INEL .LT. 4) THEN
c               if current && last ASOS inst then NO MOVE
                if(instr(34).eq.1 .and. 
     *            linstr(34).eq.1) then
                  NOTMOV(3:3) = '1'
c               then CRS
                else if(instr(34).eq.linstr(34) .and.
     *            INSTR(2).eq.1.and.LINSTR(2).eq.1) then
c                 if ASOS not changed and current && last CRS inst then NO MOVE
                  NOTMOV(3:3) = '1'
c               last - everybody else
                else if(instr(34).eq.linstr(34) .and.
     *            INSTR(2).eq.LINSTR(2) .and.
     *            instr(21) .eq. linstr(21) .and.
     *            instr(22) .eq. linstr(22) .and.
     *            INSTR(10) .EQ. LINSTR(10)) THEN
                  NOTMOV(3:3) = '1'
                ELSE
                  HISTCD = 1
                  if(instr(35).eq.1.and.
     *              linstr(35).eq.0) then
                    movstr = movstr(1:lnblnk(movstr)) // ' ASOS'
                  else
                    movstr = movstr(1:lnblnk(movstr)) // ' INST'
                  endif
                END IF
              ELSE
                IF(INSTR(15).eq.1 .and. LINSTR(15) .eq. 1) then
                  NOTMOV(3:3) = '1'
                else if(INSTR(15) .eq. LINSTR(15) .and.
     *            INSTR(5).eq.1 .and. LINSTR(5).eq.1) then
                  NOTMOV(3:3) = '1'
                else if(INSTR(15) .eq. LINSTR(15) .and.
     *            INSTR(5) .eq. LINSTR(5) .and.
     *            INSTR(9) .EQ. LINSTR(9) .AND.
     *            INSTR(12) .EQ. LINSTR(12) .AND.
     *            INSTR(19) .EQ. LINSTR(19) .and.
     *            instr(24) .eq. linstr(24) .and.
     *            instr(25) .eq. linstr(25)) THEN
                  NOTMOV(3:3) = '1'
                ELSE
                  movstr = movstr(1:lnblnk(movstr)) // ' INST'
                  HISTCD = 1
                END IF
              END IF

C             TEST LAST DISTANCE AND DIRECTION
C             TEMPERATURE MOVED);
c             with the newer keyed sources (MSHR & CDMP) - 
c               nonblank ldisdir == MOVE
              if(src.eq.2 .or. (src.eq.0 .and. icoop.eq.1)) then
                if(ldisdir .ne. "           ") then
                  movstr = movstr(1:lnblnk(movstr)) // ' LDIS'
                  HISTCD = 1
                endif 
              else IF(src.eq.0 .and. icoop.eq.0) then
                if(ldis.eq.999) then
                  movstr = movstr(1:lnblnk(movstr)) // ' LDIS'
                  HISTCD = 1
                else if((LDIS .EQ. 0 .AND. LDIR .EQ. '000') .OR.
     *             (INEL .LT. 4 .AND. (LDIS / 100) .EQ. 8) .OR.
     *             (INEL .EQ. 4 .AND. (LDIS / 100) .EQ. 9)) THEN
c                 with the USHCN if ldis != 999 then
C                 NO DIS/DIR CHANGE IF:  NO CHANGE IN DISTANCE AND NO
C                 CHANGE IN DIRECTION; TEMP inst AND DIST only if "8XX" (I.E.,
C                 PRECIP MOVE); PREC inst AND DIST only if "9XX"
                  NOTMOV(4:4) = '1'
                ELSE
                  movstr = movstr(1:lnblnk(movstr)) // ' LDIS'
                  HISTCD = 1
                END IF
              END IF

C             INITIALIZE STATION MOVE VARIABLE
              POSCHG = 0

c -----------------Post Office Position removed 08 May 2001 -----------------

C             TEST ELEVATION, LATITUDE, AND LONGITUDE
C             MOVE IF:  CHANGE IN ELEVATION, LATITUDE, OR LONGITUDE
              IF(ELEV.NE.LELEV) then
                movstr = movstr(1:lnblnk(movstr)) // ' ELEV'
                POSCHG = POSCHG+1
              endif  
                    
              aallat = abs(aLAT-aLLAT)
              aallon = abs(aLON-aLLON)
              if(aallat.gt.lleps .OR. aallon.gt.lleps) then
                movstr = movstr(1:lnblnk(movstr)) // ' LALO'
                POSCHG = POSCHG+1
              endif  

              IF(POSCHG .EQ. 0) THEN
                NOTMOV(4:4) = '1'
              ELSE
C               REAL MOVE OCCURRED
                movstr = movstr(1:lnblnk(movstr)) // ' MOVE'
                HISTCD = 1
              END IF

              if(HISTCD .eq. 1) go to 110

C             NO MOVE OCCURRED
              HISTCD = 0

c ---------------Observer change removed 08 May 2001 -----------------
            END IF

c           last shot to affect history record - is this a GPS only entry
  110       if(byear .gt. 1995 .and. movstr .eq. ' LALO MOVE' .and.
     *        llatsec .eq. 0 .and. llonsec .eq. 0 .and. 
     *        aallat .le. tolgps .and. aallon .le. tolgps) then
              histcd = 0
              movstr = ' GPS'
            endif
            
c           need to save off the beginning of the history records
            if(movnum .eq. 0) then
              movnum = movnum + 1
              call iym2imo(byear, begm, movym)
c             if the stnhist is before the beginning, set to the beginning
              if(movym .le. 0) then
                move(movnum) = 1
                mday(movnum) = 1
              else  
                move(movnum) = movym
                mday(movnum) = begd
              endif  
              if(ihdbug .eq. 1) 
     *          write(6,1001)cstnstr,istn,movnum,byear,begm,begd,movstr
            end if              

            call iym2imo(byear, begm, movym)
            if(histcd .eq. 0) then
              if(ihdbug .eq. 1) 
     *          write(6,1000)cstnstr,istn,histcd,byear,begm,movstr
 1000         format(a6,i4,' SHF STAY ',i5,1x,i4,1x,i2,1x,a)                    
            else
c             accumulate move dates from the shf records 
c              - use in next step with data
c             if there are more than one before the beginning, skip 
              if(movym .gt. 1 .and. movym .le. maxmo .and.
     *          movym .ne. move(movnum)) then
                movnum = movnum + 1
                move(movnum) = movym
                if(itimeres .eq. 0) then
c                 scale annual as if monthly
                  mday(movnum) = (begm * 30 + begd) / 12
                else  
                  mday(movnum) = begd
                endif  
                if(ihdbug .eq. 1) write(6,1001)cstnstr,
     *            istn,movnum,byear,begm,begd,movstr
 1001           format(a6,i4,' SHF MOVE ',i5,i5,2i3,1x,a)
              endif
            end if
c           test for MMTS instr and whether to update mmdates
            if(instr(21) .eq. 1) then
c             MMTS in use - set begin date if not already set
              if(mmdates(istn,1) .eq. 9999) then
c               MMDATES MUST BE A MONTHLY INDEX
                ltimeres = itimeres
                itimeres = 1
                call iym2imo(byear, begm, movym)
                itimeres = ltimeres
                mmdates(istn,1)=movym
                print *,' Turn on MMTS',istn,byear,begm,movym
              endif
c            else
c             MMTS not in use - if begin date was set and end date is
c              not set then set end date
c              if(mmdates(istn,1) .ne. 9999 .and. 
c     *          mmdates(istn,2) .eq. 9999) then
c                mmdates(istn,2)=movym
c                print *,' Turn off MMTS',istn,movym
c              endif  
            endif  

c            print *,' pass,histcd: ', ipass, ' ', histcd
c            if(ipass .eq. 0 .or. histcd .eq. 1) then
C             PAST INITIAL STATION HISTORY RECORD:  INITIALIZE PASS VARIABLE SO
C             TEST FOR REAL MOVES
              iPASS = 1

C             SAVE HISTORY DATA
              LELEV = ELEV
              LINSHT = INSHT
              LLAT = LAT
              LLON = LON
              LOBSER = OBSER
              LPODIB = PODIB
              LPODIS = PODIS
              LTMPOB = TMPOB
              linsth = insth
              allat = alat
              allon = alon
              llatsec = latsec
              llonsec = lonsec
              lbegy = begy
              lbegm = begm
              lbegd = begd
              lendy = endy
              lendm = endm
              lendd = endd
c             this is needed for both USHCN & CDMP
              if(src.eq.0) then
                lys0 = begy
                lms0 = begm
                lds0 = begd
              endif  
c            endif  

          END IF ! end of process record within begyr-endyr
        end do ! end of read station history data loop

c       align SHF moves with missing and minlen data segments
  120 if(movnum .gt. 0) then
        call alignmoves(cstnstr, istn, mounit, movnum, move, amt, mday)
      else
        print *,'SHF: ',istn,' Station has no moves:', cstnstr
      endif  
        
  150 CONTINUE ! end of station loop
      print *,nohist, ' stations have no history'

C     CLOSE STATION HISTORY FILE
      CLOSE(9,ERR=240)

      RETURN

C     PRINT ERROR MESSAGE AND ABORT PROGRAM

  160 call perror('CANNOT OPEN DATA FILE:' //  elemfile)
      GO TO 260

  180 call perror('CANNOT READ DATA FILE:' // elemfile )
      GO TO 260

  200 call perror('CANNOT OPEN HCN FILE:' // histfile )
      GO TO 260

  220 call perror('CANNOT READ HCN FILE:' // histfile )
      GO TO 260

  240 call perror('CANNOT CLOSE HCN FILE:' // histfile )
  
  260 WRITE(6,270)
c     let calling program know that the data cannot be retrieved
      inread = 0
  270 FORMAT('SKIPPING STATION')

      return
      END

C***********************************************************************
C End of Subroutine INSHP3.                                            *
C***********************************************************************
 

      subroutine alignmoves(cstnstr, istn, mounit, movnum, move, amt, 
     *  mday)
c       SHF records have been interpreted for moves, align with the data
c         in the UCP routines, the changepoint location is the LAST MONTH
c         of the SEGMENT with data BEFORE (EARLIER) the change.
c       Rules:
c         A) see if the current month can be used, depending upon the day of the move
c           1) if the day > 25 the current month is used with the earlier segment
c           2) if the 5 < day < 25 the current month is unusable (make missing)
c                and use the previous month
c           3) if the day < 5 then the current month is used in the later segment
c                and use the previous month
c         B) Adjust the move backward to the first month with data
c         C) Remove any segments without sufficient data for adjustment calculation
c               LATER: MAY BE ABLE TO USE THE #HITS TO DECIDE ANOTHER METHOD

c     subroutine arguments
c     cstnstr - station identifier
c     istn - station index
c     mounit - output file number

c     Input (raw metadata values)
c     Output (compressed series)
c       movnum - total number of moves
c       move - "ym" of each move
c       amt - amplitude of each move (unknown for SHF)
c       mday - day of month for each move

      include 'inhomog.parm.mthly.incl'
      include 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'
      
      character*6 cstnstr
      integer move(ninh), mday(ninh)
      real amt(ninh)

      integer isegym(ninh)
      real amtym(ninh)
      
      maxmo = nmo
      if(itimeres .eq. 0) maxmo = numyr

      do ichg = 1, ninh
        isegym(ichg) = 0
        amtym(ichg) = 0.0
      enddo  
      
c     go through the station history data, adjusting for the day of
c       change - REMEMBER the move is the last month of a segment
c       if day <= 5 then set move back 1 month
c       if day > 5 && <= 25 then set move back 1 month and delete month data
c       if day > 25 then leave move were it is and keep data 
      do imove = 1, movnum
        call imo2iym(iy,im,move(imove))
        if(move(imove) .gt. 0 .and. move(imove) .lt. maxmo) then
          if(mday(imove) .le. 5) then
            if(imove .gt. 1) then
c           ... if the move was within the first five days, this
c              month will go with the Next segment
              move(imove) = move(imove) - 1
c              print *,' day<5: ', move(imove)
            endif  
          else if(mday(imove) .lt. 25) then
c           ... if the move was between the 5th and 25th - this
c                 month is no good to anyone, delete it
            temp(istn,iy,im) = amiss
            tflg(istn,iy,im) = 'r'
            move(imove) = move(imove) - 1
c            print *,' day>5<25 ', move(imove)
c         ... else if the move is after the 25th - then 1) leave metadata
c               alone and 2) use this month in the THIS segment
          endif
        endif  
        if(idebug .ge. 2)
     *    print *,' Hist move: ', imove, istn, iy, im
      enddo  
      
c     test the segments between moves for enough data and adjust the meta
c       to the last non-missing value of a segment           
      jmove = 0
      sumamt = 0.0
      firstval = amiss
      ifym = 0
      iseg = 0
      imove = 1

c     Go through the data
      do iy = begyr+1, endyr
        do im = lomth, himth
          call iym2imo(iy,im,iym)
          dval = temp(istn, iy, im)

          if(dval .ne. amiss) then
            iseg = iseg + 1
            lgood = iym
          endif

c         set the first data position in series as the first move
          if(firstval .lt. amiss + 1) then
            firstval = dval
            ifym = iym
c            print *, iym, imove, move(imove)
            if(iym .eq. move(imove)) imove = imove + 1
          else if(iym .eq. move(imove) .or. iym .eq. maxmo) then
c           past the first meta - if the segment long enough
            if(iseg .ge. minlenshf) then
c             set the meta to the last good value location
              jmove = jmove + 1
              isegym(jmove) = lgood
              amtym(jmove) = amt(imove)
              lsegym = lgood
            else if(jmove .gt. 0) then
c             segment too short - delete data in segment & move shf info
              call imo2iym(ky,km,iym)
              call imo2iym(ly,lm,lsegym+1)
              if(idebug .ge. 2)
     *          print *,cstnstr,istn,' Alignmoves - Delete segment: ',
     *            ly, lm, lsegym+1,' to ',ky,km,iym
              do jym = lsegym + 1, iym
                ntest(istn, jym) = 0
                nchgpt(istn, jym) = 0
                call imo2iym(jy, jm, jym)
                temp(istn, jy, jm) = amiss
                if(imove .le. movnum) tflg(istn, jy, jm) = 'a'
              enddo
c             if annual time res - remove monthly values as well
              if(itimeres .eq. 0) then
                do jy = ly,ky
                  do jm = 1, 12
                    temp(istn, jy, jm) = amiss
                    if(imove .le. movnum) tflg(istn, jy, jm) = 'a'
                  enddo
                enddo
              endif      
              
c             accum move adj (if available)
              if(amt(imove) .ne. amiss) 
     *          amtym(jmove) = amtym(jmove) + amt(imove)
            else
c             first segment(s) not long enough - reset firstval, ifym
              call imo2iym(ky,km,iym)
              call imo2iym(ly,lm,ifym)
              if(idebug .ge. 2) 
     *          print *,cstnstr,istn,' Alignmoves - Del 1st segment: ',
     *            ly, lm, ifym,' to ',ky,km,iym
              do jym = ifym, iym
                ntest(istn, jym) = 0
                nchgpt(istn, jym) = 0
                call imo2iym(jy, jm, jym)
                temp(istn, jy, jm) = amiss
                if(imove .le. movnum) tflg(istn, jy, jm) = 'a'
              enddo
c             if annual time res - remove monthly values as well
              if(itimeres .eq. 0) then
                do jy = ly,ky
                  do jm = 1, 12
                    temp(istn, jy, jm) = amiss
                    if(imove .le. movnum) tflg(istn, jy, jm) = 'a'
                  enddo
                enddo
              endif      
              firstval = amiss
              ifym = 0
            endif  ! end of segment check at move
              
c           reset for next segment
            iseg = 0
            imove = imove + 1
c            if(imove .eq. movnum) goto 130
          endif  
        enddo ! end of month loop
      enddo ! end of year loop    

c     wrap up this station by filling the sahist array
  130 if(jmove .eq. 0) then
        if(idebug .gt. 1)
     *    print *,cstnstr, istn,' SHF/DATA too corrupt for UCP'
        goto 150
      endif
        
c     get the SHF input ready for the alignmove output
      do inh = 1, ninh
        move(inh) = 0
        amt(inh) = 0.0
        mday(inh) = 0
      enddo  
      movnum = jmove
        
      if(idebug .ge. 2) then
        print *,cstnstr,istn,' SHF/DATA move summary: '
        call imo2iym(iy,im,ifym)
        print *,' First data value: ', iy,im
      endif  
      do imove = 1, movnum
        sahist(istn,imove) = isegym(imove)
        move(imove) = isegym(imove)
        amt(imove) = amtym(imove)
        mday(imove) = 31
        call imo2iym(iy,im,isegym(imove))
        if(idebug .ge. 2) then
          if(imove .ne. jmove) then  
            print *,' End seg: ',imove - 1,' ym: ',iy,im,
     *        move(imove), amt(imove)
          else  
            print *,'   End segment ym: ', iy,im,move(imove)
          endif
        endif  
        nht = nht + 1
      enddo
        
      if(mounit .gt. 0)
     *  write(mounit,1000) cstnstr, istn, jmove, 
     *    (move(i),amt(i),i=1,jmove)
 1000 format(a6,i6,i4,30(i7,f9.3))
      
  150 return
      end
