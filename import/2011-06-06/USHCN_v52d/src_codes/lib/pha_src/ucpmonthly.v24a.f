      program ucpmonthly

c    Undocumented Changepoint Driver began as a duplicate of octpart.v5d
c      on 10 May 2004. This is the major work into Phase II of the Inhomogeneity
c      project for the USHCN.
c
c     Version    Date     Description
c
c    --------------------------- Official Version USHCNv2 ----------------------
c
c    (52d) v24a 03aug08  Minor change in debug output formats
c
c          v24  29jul08  Minor change to ADJ output
c
c    (52a) v23  21mar08  Added output for HOFN graphs via commline -H hyear
c                          with output to -O outdir
c
c    (52a) v22a 06apr07  Created 'g' output flag for UNSTBL segment to remove from
c                          fillin (shortest segment)
c                        Fixed bug that excluded stations < 24 years!
c
c    (v52) v22  05jan07  Made distinct output idents for CONSHF, UNSTBL, etc.
c
c          v20e 20dec06  No Meta loop /    Sig Loop
c    (51d 51e 51f)
c
c    (v51) v20d 20dec06  Use ONLY TPR0 adjustments in ESTAMT at NPASS
c                                  fixed sigtrnd := 0.0 in ESTAMT
c    (51a 51b 51c)          Meta loop / No Sig loop
c
c    (50s) v20c 18dec06  Remove "not sig" step again
c
c    (50r) v20b 18dec06  If a local trend is removed then remove med adj (sig or not)
c
c 51a(50q) v20a 18dec06   Added paired series positive correlation requirement
c
c v51(50p) v20  18dec06  Attempt to use a skewed interquartile sig test 1.5/1.0 
c                          Changed usage of ichgout as temporary for each pass
c                          Fixed metadata KEPT when not tested
c                          Reset lower chgpt to be evaluated
c                          Changed Not sig minsta from 2 to 3
c                          Restored "not sig" step
c
c 50f     v19f  13dec06   Estamt - reverted to 2 steps i.e. removed "not sig" 
c
c 50e(o)  v19e  12dec06   Estamt - removed "good neigh" pass for step 3
c                                 removed 10 yr restriction on trend est
c
c 50d     v19d  12dec06   reduced interquartile trim scale from 1.5 to 1.0 
c
c 50c     v19c  06dec06   Remove local trend & interim adj
c
c v50(n)  v19b  05dec06   nindx migrated to Restart journal
c                        mmts adjustments to monthly in an Annual run fixed
c
c 50a     (minslpseg=18)added step for not-significant chgpt removal before
c 50b     (minslpseg=36)      adjustment step
c
c      v19a    28nov06   correct for removal of UNSTBL changepoints
c
c      v19     12nov06   corrections from the evaluation process
c                         1) move up timeres setting for himth/lomth 
c
c      v18a3   02nov06   additions for evaluation - 
c                        1) added nhits to WM99 output in estamt
c                        2) added mmts dates for correction to read_hist
c
c      v18a2   02nov06   reset estamt npass=2
c
c      v18a1   26oct06   Set estamt npass=4 and uncommented attendant lines
c
c      v18a    20oct06   Expanded merge year range for configfilt conflation
c
c      v18     18oct06   Consolidated removal of chgpt logic in Estamt
c
c      v17j    17oct06   Modified chgpt adj zscore to min mly anom std
c                          dev(temp1,temp2) from std dev of diff series
c
c      v17i4   06oct06   Removed shock-and-recovery modification from 
c                          non-significant slopefail median decisions (ZERSLP) 
c
c      v17i3   05oct06   Moved shock-and-recovery mod to median pooled estimate
c                          for slopefail contributions
c
c      v17i2   04oct06   Remove LMQ SINE test......
c
c                        Current solution to the "shock-and-recovery"
c                          senarios that we first began to address in MAY06!
c                        06oct06 - corrected equation to:
c                          zzadj(j) = adj(j) - 0.5 * (trend(j) + trend(j+1))
c
c      v17i1   04oct06   TEST - Try the LMQ SINE test again
c
c      v17i    04oct06   reverted back to 17g - updated all of the small
c                          errors found in the interim 
c                          Retained sigtest in final adjustment (tchgpt>tcrit)
c
c      v17g2   03oct06   Added significance testing of Targ-Pair in final
c                          adjustment using tchgpt > tcrit from MINBIC(iopt=4)
c
c      v17g1   03oct06   Rewrote final adjustment for one pass that broke the
c                          series back into non-slope and slope-fail threads,
c                          lost total offset for entire series (dstep - dtrend)
c
c      v17h    02oct06   TEST - Added back LMQ SINE test
c
c      v17g    29sep06   Revert to 17e3 - and remove extreme slope test (i.e.
c                          allow all slopes) and fixed cleanup for non-slopefail
c
c      v17e4   28sep06   TEST - Removed significance testing for adjustments
c                     (take all estimates, no test for median other than zero)
c
c      v17e3   27sep06   Fixed final adj/slope loop cleanup for non-slopefail
c                          and removed a bunch of trimming for slopes
c                          and removed LMQ SINE test
c
c      v17e2   27sep06   Expanded the 1 and only pass to all MINBIC models
c
c      v17e1   26sep06   Retreat from 17f to remove all trimmed means in both
c                          amplitude and slopes BEFORE any significance testing
c                          in amplitude and slopes.
c
c      v17e    22sep06   Try adding back slpfail restrictions (i.e. 
c                        +/-0.09 degrees per month: 6 C/decade)
c
c      v17d    22sep06   Fix Tukey Median trim to +/-1.0 Inner-Quartile Range
c                          from Q1 and Q3 instead of +/-1.5 IQR from Median
c
c      v17c    22sep06   Resurrected ONLY TPR2 option for single pass
c
c      v17b    21sep06   Restricted number of adjustment passes from 5 to 1
c                        Modified Minbic option from ONLY TPR2 (6) to
c                        ALL MODELS (4) on the single pass
c                        Removed slpfail restriction
c
c      v17a    20sep06   Fixed problems with the pooled trim estimate
c                        and bypassed annual fitting for the unsupported
c                        metadata.
c
c      v17     14sep06   Attempt to deal with the annual periodicity in some
c                          of the paired difference segements in ESTAMT using
c                          the Levenberg-Marquardt optimizing fitting algorithm
c                          with the monthly function form:
c                           Y = (b1 + b2 * X) + b3 * sin(2pi/12 * (X - b4))
c
c      v16k    01sep06   All the following changes are in ESTAMT:
c                        1) Pool the chgpt before segment with the chgpt-1 after
c                          segment trends for the slope-fail adjustments.
c                        2) Shorten the station changepoint loop
c                        3) Remove artifacts of the Span Offset variation
c
c      v16j    31aug06   Reapply removal of average monthly temps from the 
c                          diff series BEFORE applying the undocumented
c                          changepoint tests (see setrnet)
c
c      v16i    23aug06   In the Estamt final IPASS looping - the test for 
c                          minimum number of neighs has been placed BEFORE
c                          the slope fail check
c                        Also, the minlen of difference series segments were
c                          fixed to IADJLEN (incorrectly was MINLENSHF)
c
c      v16h    15aug06   After a long absense - removed ANNUAL modifications
c
c      v16h.ann 27jun06  set up Dec,2000 as the "lock month" for 0.0 adjustment
c                          see LOCKMO is estamt
c                        !!!!!! BE AWARE - SETUP FOR ANNUAL !!!!!
c                        excessive slope threshold reduce to 1.5C/decade
c                          see slpmthmax at inhomog.parm.system.mthly.incl
c
c      v16g    14jun06   expand #passes in estamt from 2 (for testing) to 6
c                          (for production). Simplified dT-dS algorithm.
c
c      v16f    08jun06   implement the "delta Trend - delta Step" adjustment
c                          for NPASS in ESTAMT. This is fundamentally different
c                          than the chgpt adjustment because the trend is 
c                          considered as part of the change.
c
c      v16e    08jun06   fix ESTAMT median adjustment routine. The trimmed 
c                          median should be used.
c                        fix error in Slope Fail "lock down" setting so that 
c                          it can check segment 1. AND fix error in failtest.
c        NOTE: for testing while making verion 16f - skips NPASS in ESTAMT
c
c      v16d    07jun06   use slope failure amplitudes in estamt IFF the chgpt 
c                          adjustment is estimated (see v16c)
c
c      v16c    06jun06   backtracked to slopefail threshold of 1/2 paired slopes
c                        greater slpmthmax (instead of 1/3)
c
c      v16b    05jun06   changed assumption of "Diff domain" when MINBIC returns
c                        an SLR and TPR: to: assumption of "Same domain" unless
c                        and until an opposite sign TPR is returned.
c
c      v16a    02jun06   added 10 year limit to adjustment calculation to 
c                          reduce effects of non-stationarity of temp
c                        added "last pass" at minsta/minlen for skipping
c                          slope-fail (shock & recovery) segments AFTER all
c                          possible other segments are merged
c
c      v16     31may06   clamp down on "shock & recovery" segments by ignoring
c                          them in estamt and correcting past them
c                         1) added arrays to retain the entire network worth of
c                              adjustment information
c                         2) implemented multiple passes through the estamt
c                              adjustment estimation by starting with the 
c                              longest period w/highest number of neighbor
c                              comparisons. 
c                         3) refrained from stepping across segments that failed
c                              slope restrictions till the last pass.
c
c      v15h    30may06   expanded ndelete to represent the statest complementary
c                         test at the PAIRED station level to transfer info
c                         directly to the amplitude estimation routine
c
c      v15g    29may06   introduced nfail return from minbic for estamt segment
c                          removal (for estimation only - data carried through
c                          into adjusted series)
c
c      v15f    28may06   removing segment merging in ESTAMT step 2
c                          NOTE: minhits is changed from 5 to 2
c
c      v15e    26may06   1) added NULL domain to "do nothing" in the Interim
c                          chgpt when only straight lines are returned by MINBIC
c                          that is: the interim chgpts are removed, but no data
c                          is deleted (nor ndelete).
c                        2) removed segment merging in Combined Stat algorithm
c                        3) went back to MINBIC iopt=4 for ESTAMT step 2
c                              (estimating amplitude)
c
c      v15d    25may06   changed philosophy about management of chgpts in
c                          estamt for the ZERSHF & REMSHF events.
c                          Since it is very difficult to make decisions 
c                          on combining the difference paired results involving
c                          slope failure, best model but insignificant, and
c                          best model with significance for the same candidate.
c                          The changes are as follows:
c                          -------- before  ------  | ------- current -------
c                         1) list of offsets to form trimed mean:
c                          use sig & insig offsets  | only use sig offsets
c                         2) when the trimed mean includes 0 (ZERSHF)
c                          keep chgpt for brkslope  | remove as straight line
c                         3) after removing a ZERSHF or REMSHF chgpt
c                          proceed to earlier chgpt | redo previous chgpt
c
c      v15c    24may06   fixed statsubs to reinitialize rtraw difference series
c                          after stat tests and before minbic
c
c      v15b    16may06   expanded estamt to include the verification or removal
c                          of the unsupported metadata locations with full
c                          BIC model runs with slope/critvar tests (istep=1)
c                          Have to move the 4th filter (alignmoves) from the 
c                          confirmfilt routine to the estamt routine AFTER
c                          the ndelete > ndellim removal of temps
c
c      v15a    14may06   also added ndelete to remove
c                         data for the amplitude estimate but use all of the 
c                         data fields for the output (duly indicating those that
c                         failed the ndelete test). ndelete threshold is ndellim 
c
c      v15     13may06    may changes to keep minbic as the defining analysis
c                         tool for changepoints. 
c
c      v14d    11may06    changed minbic option to "2" (critvar test) for
c                           interim chgpt tests. Kept "4" (no critvar nor
c                           slope test) for estamt
c
c      v14c    10may06    expanded "interim" chgpts idea using 
c                          1) if chgpts have same sign, pick best
c                          2) else remove interim chgpts & data in paired diff
c                                then recalc minbic to earliest chgpt
c
c      v14b    09may06    added ndelete to reinforce QC
c
c      v14a    09may06    removed slope/critical value BIC test from the 
c                         amplitude estimation. Pass all estimates back.
c
c      v14     06may06    realized estamt uses minbic to calculate offset 
c                         and t/ftests. !!!! QED: no need for method 3 !!!
c                         Changed array for chgpts from nchgpt to nhits
c                         removed division of nhits by ntech for thresholds.
c                         Added matrix printout after 4th filter(alignmoves)
c
c      v13     04may06   removed ntech dimension from nfound & nspan arrays,
c                        standardized point indices for first (target) and
c                        second (pair) stations - itarg,ipair, it2pair, ip2targ
c
c      v12a    21apr06   reworked statsubs for new philosophy
c
c      v11a    17apr06   Brought neglected annual timeres back to operational
c
c      v11     06apr06   For testing - changed slpmthmax == .015 (betav6)
c
c      v10     05apr06   For testing - reinstated ability to run each 
c                          technique (-T option)
c
c      v9f2    23mar06   revived completeness check for segments (betav6)
c
c      v9f1    23mar06   set minser back to 24 (betav6) instead of 12 (betav11)
c                          in inhomog.parm.system.incl
c
c      v9f     16feb06   rearranged minseg length in estamt to only minlenshf
c                          (03mar06 minlen == 18)
c
c      v9e     09feb06   added sum-square error return from MINBIC
c                          and fixed the UNSTBL at beg/end of series
c
c      v9d     22jan06   trivial change - added inqtype=5 to minbic in estamt
c
c      v9c     21jan06   moved ntest into the MLRTest to add only when a
c                          station has data
c
c      v9b     16jan06   rearranged UndocChgPt filter in confirmfilt. Now
c                          absorbs lesser hit months into closest - higher
c                          hit month. NOT highest amp in bracket.
c
c      v9a     15jan06   attempt to fix ntest wrt number of real tests for 
c                          the "neigh-cand" pair (vs the "cand-neigh" pair)
c                          Alternative is to add extra neighbor indices to
c                          all the work arrays, and there is just not enough
c                          space nor time.
c
c      v9      13jan06   brought back model type for each chgpt from MLRTest
c                          in order to add a model based confirm filter
c
c      v8H     27dec05   for speed sake: set WM = 99 first loop
c
c      v8g     18nov05   FINAL ESTAMT BIC CALL:
c                        Good Neigh restricted to pairs with adj > inner quart
c                        and segments > minlenshf (the know chgpts tests)
c
c      v8f     11oct05   remove slope restiction on BIC models 2,4,5 for a 
c                        more robust fit
c
c              29sep05   method = 3 uses the BIC only to find the best
c                        amplitude estimate, not to determine whether there
c                        is a chgpt.
c
c      v8e     28sep05(ii) method = 3 is removed from confirmfilt and goes
c                        directly to estamt (as in the beginning!)
c
c      v8d     28sep05   rearrange the use of the SHF - read the SHF
c                        along with the alignmoves at the beginning of
c                        main, right after reading the networks data.
c                        This usable segments in the data. Do not merge
c                        the metadata with the undoc-chgpts until the 
c                        confirmfilter stage so the pairwise can do its job.
c
c      v8c     22sep05   include data gaps > minlen as "raw" chgpts
c                        fixed UCHGPT filter for all hits (error @ 1914)
c
c      v8b     20sep05   rework the estamt code to determine whether a chgpt
c                          has enough significance and the final adjustment
c
c      v8      07sep05   Reworked statistical tests to strengthen independence
c                          and completeness:
c                          1) removed the Lag1 autocorr test
c                          2) added the Xioulan constant non-zero slope model
c                          3) removed the t-test, then
c                          4) added the more general Bayesian Info Criteria
c                               to test which model is "best" by estimating 
c                               the effect of the model's complexity
c
c      v7      09aug05   major work involving "fixing" a changepoint to a
c                          given stn/yr/mth no matter the missing data and
c                          pairwise compression. 
c
c                          nfound changed to char*1 & add nspan
c                          and consolidated many iy/im indices into imo
c
c      v6k     02aug05   If the data is missing at a station history "hit" 
c                          then shift the SHFHIT to the first month with data
c
c      v6j     22jul05   Correlation weighting for the adjustment estimate
c                        just does NOT work. Try the inner-quartile range
c
c      v6i     22jul05   Changed minlen from 5 to 18
c                        Upgraded the difference from full temperature to
c                          monthly anomalies (wrt whole series)
c                        Added a toggle to examine effect of correlation of
c                          monthly anomaly on the adjustment estimation
c
c      v6h     20jul05   Added completeness check for the adjustment estimation
c                        attempting to reign in the variation from pair to pair
c
c       -------------------- 2005_beta1 ---------------------
c       v6g    11jul05   At last filter - window set to smallest size (=5)
c
c       ------------------- 2005_alpha ----------------------
c       v6f    30jun05   Several minor modifications to optimize algorithm
c                        1) FINSHF filter changed from amp > std to
c                             t-test homog/inhomog
c                        2) method = 3 threshold set to 1
c                        3) corrected imo index conversion to/from mth-yr
c
c       v6e    26may05   test to reintroduce the t-test (method == 3) 
c
c       v6d    26may05   attempt to span window gap in estamp
c
c       v6c    26may05   retracted the 95% to 90% inclusion (these are what
c             the thresholds are based upon)
c
c     ucpmonthly.v6b 20may05:
c          1) SHFHIT (pass1-2nd) found did not use std. amp. - fixed
c          2) UCHGPT (pass1-3rd) filter 
c               Use the sorted amplitude to conflate UndocChg
c          2) Use 95% inclusion for the amp vs range spread - done
c          3) Add t-test (SHAP type) check for unresolved ChgPts at 
c               Pass2 before the amplitude filter
c          4) Currently - leave amplitude avg > std threshold as last filter
c               modify with weighting wrt 1diff corr
c
c     ucpmonthly.v6a 18may05 
c             adjusts offset amp wrt avg mthly stddev for width of merge window
c
c     ucpmonthly.v6 10may05:
c             Added output summary for each candidates estimate in "estamt".
c             This effectively removes method=3 which now has no calls to 
c             the MLRTest nor confirmfilt.
c
c     ucpmonthly.v5a 24mar05 
c             The clearing algorithm for the "highest erasure filter" in
c             confirmfilt did not remove all of the paired hit increments
c             because it back indexed from the highest hit station neighbor
c             list. Rewrote to look at all of the neighbors of all the 
c             stations, because station(2) may be a neighbor to station(1),
c             but station(1) does not necessarily have to be a neighbor of
c             station(2)!
c
c             Also had to set the minimum location spread vs amplitude to
c             the minimum window size for the statistical tests. Any other
c             size amplitude will "disappear" if there is any undoc or
c             not undoc changes within the minimum window. COULD SPAN THE
c             WINDOW!
c
c     ucpmonthly.v5 16mar05 Modified method 3 from the SHAP algorithm to 
c             merging in the metadata info into the suspect changepoint
c             array and employing the pairwise technique with the T-test.
c             Still holding onto method 2 (composite reference stations) as 
c             it may be useful at some future date.....
c
c     ucpmonthly.v4 16mar05 MLR Testseg vs SHAP Stnbn vs SHAP Mrgbn chgpt 
c             comparison. Shows that a new direction with the MLR Testseg
c             output gives BETTER results. This is a locked-in version with
c             all the associated models locked-in.
c
c     ucpmonthly.v3a 07mar05 Added results of the research into parameters for
c             pairwise (including confirm_filter.v2.f). As of UCPM3aM11bS41 and
c             its associated routines, there are no more confirmation input
c             parameters. For work related to the parameter study, see version:
c             UCPM3Lm2bM11aS41.
c
c     ucpmonthly.v3  07Feb05 Removed SHAP to speed up runs for pairwise 
c             questions. See references to stnhist. Added variable monthly
c             merge depending upon the estimated amplitude of the changepoint
c             (3 mar 05)
c
c     ucpmonthly.v2  02feb05 Removed all calls to the LV package(s). Mainly
c             removed the monthly estimation after the SHAP final call and the
c             entire multi-linear-regression options - not needed in pairwise
c             algorithm. 
c             Also, removed Durbin-Watson in favor of AutoLag1 only, because
c             the threshold can go down to lower obs and no need to defend 
c             whether using "gray" area between upper and lower bounds as in DW.
c
c     ucpmonthly.v2  18jan05 add SHAP history input (including capability for
c             input from MM's meta files)
c
c     ucpmonthly.v1o 10Dec04 tests the "high hit erasure" technique
c
c     ucpmonthly.v1n 5nov04 is the next generation using MONTHLY data after
c              Annual changepoints have been identified. Major changes include
c              expansion of the nfound array to include monthly + annual. And
c              a process toggle variable itimeres ( 0=annual;1=monthly )
c
c     ucpdriver.v1n 30Jul04 NWS48 change removed restriction on output for
c               No changepoints  
c               need output file whether chgpt found or not
c
c     ucpdriver.v1 10May04 - Testing code for the "looping" questions
c              v1 revives the iterative scheme and the SHAP adjustment
c
c     octpart.v5d 26apr04 - work around the 10 station limit for MM's
c              subnets
c
c     octpart.v5c2 9apr04 - pairwise with no loops 
c SHAP
c
c     octpart.v5b2 5apr04 - Pairwise with no loops
c
c     octpart.v5b 05apr04 - Expanded nfound for hits to include the
c             "technique" dimension
c
c     octpart.v5a 12mar04 - Beginning of the formalized testing of the 
c             pairing algorithm, first by testing the individual ref/stat
c             combo's with each benchmark case.
c
c     octpart.v5 14jan03 - major re-examination of the import of the
c             separate methods (i.e. LV, PEEP, MD) as compared to the
c             separate of techniques to generate reference series (1stdiff,
c             LV station avg, MD z-score avg) and separate statistical tests
c             (Fmax, Tmax, Durbin-Watson) so that we now have the following
c             matrix:
c                           | Fmax  | Tmax | D-W  |
c                           -----------------------
c                   1stdiff |  PS   |  PT  |  PW  |
c                   LVref   |  LS   |  LT  |  LW  |
c                   MDref   |  MM   |  MD  |  MW  |
c
c             Note: PE is a conglomeration of (PS or PT) and MRPP
c                   this will not be addressed past version 'v4d'
c
c             Note: Lucie Vincent's algorithm implements the Durbin Watson 
c                test in two ways... 
c                1) LW uses the original code from L. Vincent - the
c                   reference series is generated as average anomaly
c                2) LV is developed from L.V.'s paper - the reference 
c                   is the independent stations compared with 
c                   Multi-linear Regression
c
c     v4dEP       09jan04 - parameter study on the PEEP method...
c             Original PEEP contained F(3,4),StuT,MRPP,dT(2,3)
c             Current PEEP inplementation Fmax,dT(avg)
c             Where F(3,4) is the F(3,n-4) thresholds
c                   Fmax are the Fmax thresholds
c                   StuT is the Student's T-test
c                   MRPP is the Mielke Ramdom Perturbation Test
c                   dT(2,3) is the generation of Ref_dT using 2-3 stations
c                   dT(avg) is Ref_dT from avg of all stations dT
c             Input parameters have been expanded to run tests separately
c                   or in combination.
c
c     octpart.v4d 03dec03 - lowlim = 2 for method = 1 worked for the
c             6 station benchmarks provided by MM. But, with a full
c             subnet of 20 stations, it provides far too many false
c             positives. Have returned to the idea of using a fraction
c             of the number of pair-wise "hits" to make a better 
c             measure. 
c
c            v4d also contains the first transfer and use of the 
c            METHOD type into the STATSUBS routine - just AUTOLAG currently
c
c     octpart.v4c 4nov03 - fixed bugs in random series (benchmark) code
c        for updates made since v4 began.
c
c     octpart.v4b 24oct03 - expanded Journal files to include saving and
c        checking of size parameters for state arrays that are written/read
c
c     octpart.v4a 18oct03 - incorporates LV's monthly adjustments into the
c        SHAP after the final annual adjustments are completed
c
c     octpart v4 15oct04 - ADDED RESTART ABILITY (Journal files)
c        implements the saving and reading of critical variables and arrays
c        used to restart the process at two points.....
c        METHOD_LOOP - entry point at the top of the loop defined as:
c             do while (method .lt. 3 .and. lht .ne. nht)
c          
c
c     octpart v3 26sep03 - NEXT GENERATION - IN PREPARATION FOR THE REAL WORLD
c        implements the expansion of the looping to include an entire network
c        includes first logic for the non-circularity of candidate subnets
c        (neighborhoods). Arrays and variables that drive this code are:
c
c        Pointers linking the stations and their data in each neighborhood
c          from /vapor5/USHCN/src/ushcn_corr_2000.v2.f
c        nindx(#subnets in network, #neigh in each subnet + circularity buffer) 
c
c     octpart v2b 17sep03 - calls mlrtest7, an abbreviated version of MM's mlr
c        which removes all vestiges of monthly processing and includes some 
c        subtle, but critical indexing changes for near perfect agreement with
c        MM's monthly and annual tests.
c
c     octpart v2a 16sep03 - added begin and end station input for rerunning for
c        recovery and special studies.
c
c     octpart v2 26aug03 - instituted the iexact and "last found" (lf??) arrays
c        because for_paper_I/case3b and caus_mest/3sigma got stuck in infinite loops.
c
c     octpart v1c as of 21 Aug 03 - consolidated and standardized debug output 
c             24 Jul 03 - First Concensus Version using the original 
c        Quadpart algorithms. Major changes in reponse to the Canadian Urban
c        question raised by RV from TK. Some reduction in the parameters 
c        passed and institution of a need-to-read array reduced execution
c        demands and defined the higher level looping better. 
c
c     octpart v1b 22 Jul 03 - Annual only (with month flags) defined
c        in readnet (read_write.f) - includes new versions of statistic
c        routine wrappers as well.....
c
c     octpart v1a uses the "pairing" method to find the largest chgpts, then
c        attempts the "greedy" method for finer resolution
c        (Note: original >unfinished< coding at quadpart.v2c.f)
c
c     octpart v1 16Jul 03 combines the classic statistical inhomog techniques
c        i.e. LV, EP, MD, and MM, with the mix&match reference series (Lucie's 
c        autolag and Peterson's firstdiff) with stat tests (SNIT and
c        TWOPHREG) suggested by mmenne, i.e. LS, LT, PS, and PT.
c
c     quadpart.v2b 10Jul03 implements a firmware switch (netloop) to  
c        enable(=1)/disable(=0) looping through the network to untangle 
c        chgpts in the neighborhood. Needed for benchmark statistics and
c        paper_I.
c
c     quadpart.v2a 26Jun03 Explores using only 1 network station ("pairing")
c
c     quadpart.v2 23Jun03  Explores the "round-robin" approach to determining
c         whether an inhomog is due to the candidate or network.....
c         unable to sort out attribution - abandoned approach
c
c     quadpart.v1 27May03  Added the two phase regression technique to the
c        MD suite for testing
c
c     tripart.v2 07may03   Pulled out LV's autolag subroutine so that major
c        changes can be made to the Test Statistic parsing. 
c
c     Tripart.v1 24dec02   After initial testing of the detection of the times
c       for inhomogeneity, this version 1) turns on output (for LV & PE) and
c       has modifications to run the GHCN through the PE code.
c
c     PROTO1            First program to test, head to head, the inhomogeneity 
c       techniques developed in the last 10 years since the original USHCN.
c
c     The first techniques to be tested will be:
c       PEEP: the first difference technique developed by Peterson & Easterling
c       LV: the lag correlation technique developed by Lucie Vincent
c       MLR: the maximum likelyhood ratio (formerly the standard normal homog)
c               technique developed by Alexxanderson and improved by Menne and Duchon

c     This include file contains the common parameters used by all the various routines
      INCLUDE 'inhomog.parm.mthly.incl'  
      
c     This include file contains the command line parameters
      include 'inhomog.comm.mthly.incl'
      
c     This include file contains common parameters for the MLR tests
      include 'inhomog.MDparm.mthly.incl'

c     This include file contains data required for restarts
      include 'inhomog.restart.mthly.incl'

      character*132 version
      
c     network station temporary (round robin) list (base station first)
      character*6 ntstn(maxstns), rrstn(nstns)
      
c     Raw monthly temperature data array (stn/yr/molowest)
      real orig(maxstns,begyr:endyr,13),rrtemp(2,nmo)
      character*1 rrflg(nmo) 
c     For method == 2, use(=1),not-use(=0) neighbor array
      integer mxstn(nstns)
      
c     since the network is NOT circular - but to keep down the number
c       of paired iterations - keep track of pairs done.
      integer ipaird(maxstns, nstns)
      
c     NEED TO RUN array - if, for some strange reason, there is no data for
c       a station in the meta or corr files then skip (=0) else process (=1)
c     Added intr to keep track of the subnets read in, must keep the index
c       and the stations in sync, otherwise nothing works!!!!
      integer intr(maxstns),rnunit
      
c     array of MMTS dates from read_hist first is when going to MMTS,
c         second is leaving MMTS, no MMTS is 9999. Use with comm option immts
      integer mmdates(maxstns,2)
      
      character*132 isfile
      
c     restart network (input) used to read input file up to restart point
      integer rsnet
      
C     UNITS 11, 12, 13, 14, & 15 are RESERVED!!!!! 
c      Subroutines must NOT use these units except when passed as arguments!
c     input unit for candidate - reference network (drives process)
      nnunit = 11
c     input unit for data input (used many times - open/close when using)
      idunit = 12
c     output unit
      iounit = 13
c     output plotting data (opened in main, accumulates plot data)
      ipunit = 14
c     Matt metafile for generated series
      mmunit = 15 
c     matt datafile for generated series
      rnunit = 16
c     Output file for compressed metafile from Matt Gen series for skill scores
      mounit = 17
      
c     netloop enables(=1)/disables(=0) network looping after candidate chgpts identified
      netloop = 0
      
c     Counter for random series pseudo-numbers for stations
      irancnt = 0

c     Initialize begin and end stations
      begstn = 0
      endstn = 999999

c     Initialize begin and end networks
      firstnet = 0
      lastnet = 999999
      
c     ---------------------------------------------------------------------
c     Read the command line options
      call commline()
      
c     begin with annual
      if(itimeres .eq. 0) then
        lomth = 13
        himth = 13
      else  
c       monthly
        lomth = 1
        himth = 12
      endif  
      print *,'lo-hi mth set: ', lomth, himth

c     open Matt's metafile for generated series
      if(irandom .ne. 0 .or. igra .ne. 0) then
        open (unit=mmunit,file=mattmeta, err = 215)
c        print *,' Opened: ',mmunit,' File: ', mattmeta
        isfile = mattmeta(1:lnblnk(mattmeta)) // '.compress'
        open (unit=mounit,file=isfile) 
      endif

c     ---------------------------------------------------------------------
c     Open Candidate stations meta file to process
c     -- CHANGES BY DANIEL ROTHENBERG <darothen@mit.edu> -- BEGIN
c        -- gfortran compiler expects second named argument to be called
c        -- 'file', not 'name'. 
      open(nnunit, file=netfile, err=200)
      if(mattdata .ne. '') then
        open(rnunit, file=mattdata, err = 210)  
c     -- CHANGES -- END
      endif  
      
      print *,'  Decisions on inhomogeneity of the Candidate are made'
      print *,'    with the threshold limit after the statistic value'
      print *,'    > indicates the statistic value GREATER than the'
      print *,'    limit is inhomogeneous and < for LESS than'
      print *,''
      print *,'  The following network dependent parmeters are set'
      print *,'  by the inhomog.parm.mthly.incl source file'
      print *,'   project = ', project
      print *,'   begyr = ', begyr
      print *,'   endyr = ', endyr
      print *,'   nstns = ', nstns
      print *,'   maxstns = ', maxstns
      print *,'   maxfound = ', maxfound
      print *,''
      print *,' Definitions from inhomog.parm.system.mthly.incl'
      print *,'  numyr = ', numyr
      print *,'  nmo = ', nmo
      print *,'  minser = ', minser
      print *,'  minann = ', minann
      print *,'  minlenshf = ', minlenshf
      print *,'  slpmthmax = ', slpmthmax
      print *,'  minlen = ', minlen
      print *,'  mincomp = ', mincomp
      print *,'  compt = ', compt
      print *,'  minslpknt = ', minslpknt
      print *,'  minsta = ', minsta
      print *,'  minhits = ', minhits

c     if the journal file is not blank then this is a restart run
c     restart is not used for random series
c     currently (V52) only restart = 3 is enabled
      rsnet = 0
      rentry = 0
      if(jrnlfile.ne.'' .and. irandom.eq.0)
     *  call read_restart(rsnet)
     
c     The input loop used to be (prior to v3) the major loop through the 
c     process. The philosophy used to be that each sub-net was a complete
c     network, but in the real world, it is not. A network is a collection
c     of sub-nets (neighborhoods) centered in a correlated sense, around
c     each station in the network.
c     In the real world there is only one network (usually)
      ieof = 0
      network = 0
      do while(ieof .eq. 0)

c       do not initilize arrays read in by restart!!!!
   10   if(rentry .eq. 0) then
          do imo = 1,nmo
            do itarg = 1, maxstns
              nhits(itarg,imo) = 0
            enddo
          enddo  
        
c         initialize paired network index 
          do itarg = 1, maxstns
            do istn = 1, nstns
              nindx(itarg, istn) = 0
            enddo
          enddo
        endif  
              
c       Read in one complete network - all candidates and their neighborhoods
c       nhits is used to set a chgpt at any gap in the data > minlen 
c         contiguous months.                                  22Sep05
        call readnet(mmunit,nnunit,idunit,rnunit,orig,intr,ieof,ntstn,
     *    tflg)
        if(ieof .ne. 0) goto 999
        network = network + 1
c       if restart has been initiated read input up to the restart network
        if(rsnet .gt. network) go to 10
c       only process networks between first and last
        if(firstnet.ne.0 .or. lastnet.ne.0) print *,'Network:',network
        if(network .lt. firstnet) go to 10
        if(network .gt. lastnet) go to 999

c       If past the restart network, then initialize restart arrays
        if(rsnet .lt. network) then

c         copy the original data to a temp array for finding inhomog
          do itarg = 1, maxstns
            do iy = begyr, endyr
              do im = 1, 13
                temp(itarg, iy, im) = orig(itarg, iy, im)
              enddo
            enddo
          enddo

c         initialize arrays for SHP adjustment "history"
          do inh = 1, ninh
            do itarg = 1, maxstns
              sahist(itarg,inh) = 0
            enddo
          enddo
          
c         initialize mmts dates
          do idate = 1, 2
            do itarg = 1, maxstns
              mmdates(itarg,idate) = 9999
            enddo
          enddo
          
c         fill the sahist array with station history data
c          (sahist is not initialized in read_hist)
          if(ihyear .le. 0) then
            call read_hist(mmunit,mounit,ntstn,nht,mmdates)
c           adjust for MMTS as in USHCN v1
            if(immts .eq. 1) then
              call mmts_adj(mmdates)
            endif
          endif  

c         initialize need-to-run to input
          do itarg = 1, maxstns
            ntr(itarg) = intr(itarg)
          enddo       

          nloop = 0
          method = 0

c         Initialize the total number of changepoints found (NHT) 
c                    and the number of changepoints in the last interation (LHT)
          lht = 0
c ------ need to change nht == metahistory changepoints for v2a
          nht = 0
        endif  

c       initialize round robin flags
        do imo = 1, nmo
            rrflg(imo) = ' '
        enddo

c       if restart is on and rentry == 3, skip all the way to the 
c         "Final stnhist" call
        if(rentry .eq. 3) goto 8090

c       ---------------------------------------------------------------------
c       Process for the candidate ONLY - this option is used primarily for
c         benchmark and individual stat test debugging
        if(netloop .eq. 0) then
        
          print *,' Init nfound & nspan'
          do itarg = 1, maxstns
            do ipair = 1, nstns
              ipaird(itarg, ipair) = 0
              do imo = 1, nmo
                nfound(itarg, ipair, imo) = czero
                nspan(itarg, ipair, imo) = czero
                ndelete(itarg, ipair, imo) = czero
              enddo
            enddo
          enddo

          print *,' Init ntest, schgpt, zchgpt, nchgpt'
          do itarg = 1, maxstns
            do imo = 1, nmo
              ntest(itarg, imo) = 0
              schgpt(itarg, imo) = 0.0
              zchgpt(itarg, imo) = 0.0
              nchgpt(itarg, imo) = 0
            enddo
          enddo

c         Method set to 3 because no looping
          method = 3
          
          do k = 1, nstns
            movnum = 0
            itarg = nindx(k,1)
            if(itarg .ne. 0) then
              call setrflg(itarg, rrflg, tflg)
              call setrnet(k, itarg, rrtemp, temp, lowest)
              if(k .eq. 1 .and. lowest .lt. minser) then
                print *, k,' ',ntstn(k),' not enough data - skipping'
                goto 120
              endif  
              numstn = k
            else  
              call inirnet(k, rrtemp)
            endif
          enddo
          rrstn(1) = ntstn(1)
          rrstn(2) = 'Noloop'
          call statsubs(1, 1, 0, -1, idunit, rrtemp, rrflg, rrstn)
 
c         add up the number of candidate hits for each iym
          itfound = 0
          do imo = 1,nmo
            if(nfound(1,1,imo) .ne. czero) itfound = itfound + 1
c           if 2 or more, make a suspect changepoint for SHAP
            if(itfound .ge. 2) then
              movnum = movnum + 1
              sahist(1,movnum) = imo
              print *,' Candidate Only LFOUND:',iy,im,itfound
            endif  
          enddo
          
        else
        
c         -------------------- METHOD LOOP ---------------
c         Process the network through the statistical tests until there are
c           no more significant changepoints found in any station
          do while (method .lt. 3)
 
c           In UCPM5a - sped up compute with lotech & hitech...
            lotech = 1
            hitech = ntech

c           if rentry == 2 skip to the confirmfilt call entry point
            if(rentry .eq. 2) goto 8080

c           else continue from here with current values from the
c             restart file or the regular process

            if(rentry .ne. 1) nloop = nloop + 1
c           if current loop hit total = last loop hit total then
c           increment method and initialize as needed
            if (nht .eq. lht .or. nloop .gt. 1) then
              if(rentry .ne. 1) then
                method = method + 1
c               skip method 2
                if(method .eq. 2) method = 3
                print *,' Method changed to ',method
                nloop = 1
              endif  

              if(method .ge. 3) then
                lotech = ntech + 1
                hitech = ntech + 1
              endif    

              print *,' Init ntest, schgpt, zchgpt, nchgpt'
              do itarg = 1,maxstns
                do imo = 1,nmo
                  ntest(itarg, imo) = 0
                  schgpt(itarg, imo) = 0.0
                  zchgpt(itarg, imo) = 0.0
                  nchgpt(itarg, imo) = 0
                enddo
              enddo

              print *,' Init ipaird, nfound & nspan'
              do itarg = 1,maxstns
                do it2pair = 1, nstns
                  ipaird(itarg, it2pair) = 0
                  do imo = 1, nmo
                    nfound(itarg, it2pair, imo) = czero
                    nspan(itarg, it2pair, imo) = czero
                    ndelete(itarg, it2pair, imo) = czero
                  enddo
                enddo
              enddo

c             initialize round robin array
              do is = 1, 2
                do imo = 1, nmo
                  rrtemp(is, imo) = amiss
                enddo
              enddo 
              
c             this IF skips the pairwise and CRS looping and goes
c               directly to the SHAP - commented out at UCPmonthly.v5
c              if(method .eq. 3) goto 100
            endif  
            print *,' LOOP = ', nloop

c >>>>      Write the "method/loop" journal file for restart
c            if(rentry .eq. 0 .and. unique .ne. '' .and. irandom.eq.0) 
c            if(rentry .eq. 0 .and. unique .ne. '') 
c     *        call write_restart(1)
c >>>>      CRITICAL point for restart - come here if rentry == 1
            rentry = 0

            lht = nht

c           ------------- Pairing Method (#1) and (#3 as of UCPM5)------------
c           Each station is pair with every other, concensus determines which
c           stations have which prospective chgpts and when
            if(method .eq. 1 .or. method .eq. 3) then
              
c             Use only one network station
c             keeping the year and amount of inhomogeneity for each iteration
              do k = 1, numsubs
                itarg = nindx(k,1)
c                print *,' k, itarg:',k,itarg
                if(ntr(itarg) .eq. 0) goto 55
c               setup the temporary flag array for the candidate
                call setrflg(itarg, rrflg, tflg)

                do it2pair = 2, nstns
                  ipair = nindx(itarg, it2pair)
c                 finished with neighborhood if neighbor index == 0 or
c                   there is no data for the "pair" station
                  if(ipair .eq. 0 .or. ntr(ipair) .eq. 0) goto 55
c                 set the temporary "round-robin" arrays for the candidate
                  rrstn(1) = ntstn(itarg)
                  rrstn(2) = ntstn(ipair)
                  if(itarg .gt. ncand .and. 
     *              ntstn(itarg) .gt. ntstn(ipair))then
                    rrstn(1) = ntstn(ipair)
                    rrstn(2) = ntstn(itarg)
                  endif
                  
c                  if(ntstn(itarg) .eq. '011084' .and. 
c     *              ntstn(ipair) .eq. '016640')  then
c                      print *,' Stop for a moment'
c                  endif    

                  if(idebug .ge. 1)
     *              print *,' Pair ', rrstn(1), ' with ', rrstn(2)

                  if(ipaird(itarg, it2pair) .eq. 1) then
                    if(idebug.ge.2) print *,' Preprocessed ',rrstn(1),
     *                ' with ',rrstn(2)
                    goto 50
                  endif  
                
                  call setrnet(1, itarg, rrtemp, temp, lowest)
                  if(lowest .lt. minser) then
                    print *, itarg, rrstn(1),' excess miss 1 - skipped'
                    goto 50
                  endif  
                  call setrnet(2, ipair, rrtemp, temp, lowest)
                  if(lowest .lt. minser) then
                    print *, ipair, rrstn(2),' excess miss 2 - skipped'
                    goto 50
                  endif  
                  
c                 check to see if stn1 is a neigh of stn2 (ip2targ=-1 if not)
                  ip2targ = -1
                  do kstn = 2, nstns
                    js = nindx(ipair, kstn)
c                   if there is no value in the pair array, loop out
                    if(js .eq. 0) goto 30
                    if(js .eq. itarg) then
                      ip2targ = kstn
                      if(idebug .ge. 2)
     *                 print *,' Linked ', ntstn(itarg), itarg, it2pair,
     *                    ' from ',ntstn(ipair), ipair, ip2targ
                    endif
                  enddo
c                 if not back linked, add target to pair neighbor list
   30             if(itarg .le. ncand) then
                    print *,' No link back for HCN'
                  else if(ip2targ .eq. -1) then
                    if(kstn .le. nstns) then
                      nindx(ipair,kstn) = itarg
                      ip2targ = kstn
                      print *,' Link added ', ntstn(itarg), itarg,
     *                  it2pair, ' from ',ntstn(ipair), ipair, ip2targ
                    else
                      print *,' Link full ', ntstn(ipair), ipair, 
     *                    ' at ', kstn
                    endif
                  endif  
                  
c                  print *,' 1:itarg,temp: ', itarg,temp(itarg, 1950, 1)
c                  print *,' 1:ipair, temp: ', ipair, temp(ipair, 1950, 1)
                  call statsubs(itarg, it2pair, ipair, ip2targ, idunit,
     *              rrtemp, rrflg, rrstn)

c                 set targ-neigh pair index to DONE
                  ipaird(itarg, it2pair) = 1
                  do imo = 1, nmo
                    if(rrtemp(1,imo) .gt. amiss+1.0)
     *                ntest(itarg, imo) = ntest(itarg, imo) + 1
                  enddo  

                  if(method .ne. 3) then
c                   nfound has been set for the "it2pair in itarg", now use
c                     the cross-pair info to set the "itarg in it2pair"
                    if(ip2targ .gt. 0) then
                      ipaird(ipair, ip2targ) = 1
                      do imo = 1, nmo
                        nfound(ipair, ip2targ, imo) =
     *                    nfound(itarg, it2pair, imo)
                        nspan(ipair, ip2targ, imo) = 
     *                    nspan(itarg, it2pair, imo)
                        ndelete(ipair, ip2targ, imo) = 
     *                    ndelete(itarg, it2pair, imo)
                        if(rrtemp(2,imo) .gt. amiss) 
     *                    ntest(ipair, imo) = ntest(ipair, imo) + 1
                      enddo  
                    else
c                     No link back in paired stn
                      print *,' No link back to ',ntstn(itarg),
     *                  ' from ',ntstn(ipair)
                    endif  
                  endif ! end if method != 3
   50           enddo ! end do all neighbors
   55         enddo ! end do all candidates

c    !!!!!! Daniel
c            call exit()
c

            else if(method .eq. 2) then
c           --------------- Composite Reference Method (#2) ------------------
c           ****** Abandoned as of ucpdriver.v1 - kept "just-in-case" ********
c           ---- may use in future for temporary local trend corrections ----
c             Process each candidate against its sub-network
c             keeping the year and amount of inhomogeneity for each iteration

c             initialize found array for each loop for 
c             method 2 inhomog discovery
              do itarg = 1,numsubs
                do it2pair = 1, nstns
                  do imo = 1,nmo
                    nfound(itarg, it2pair, imo) = czero
                    nspan(itarg, it2pair, imo) = czero
                  enddo
                enddo
              enddo

              do k = 1, numsubs
                if(ntr(k) .eq. 0) goto 60
                itarg = nindx(k,1)
c               setup the temporary flag array for the candidate
                call setrflg(itarg, rrflg, tflg)
              
                print *,' Candidate for merged method:',ntstn(itarg)

c               determine the best "merged network" for this candidate...
c               fill up the data arrays
                do it2pair = 1, nstns
                  ipair = nindx(itarg, it2pair)
c                 finished with neighborhood if neighbor index == 0 or
c                   there is no data for the "pair" station
                  if(ipair .ne. 0 .and. ntr(ipair) .ne. 0) then
c                   setup temporary data array
                    call setrnet(it2pair, ipair, rrtemp, temp, lowest)
                    numstn = it2pair
                  else  
c                   reiniti the temporary data array
                    call inirnet(it2pair, rrtemp)
                  endif  
                enddo  
                
c               find "best" neighbor combination using max# station-years 
                call maxnet(numstn, rrtemp, mxstn)

c               fill round-robin arrays with the best "merged network" 
c               fill up the data arrays
                it2pair = 0
                numstn = 0
                do ind = 1, nstns
c                 reiniti each stn in temporary data array
                  call inirnet(ind, rrtemp)
                  if(mxstn(ind) .eq. 1) then
                    it2pair = it2pair + 1
                    ipair = nindx(itarg, ind)
c                   refill temporary data array for mxstn
                    call setrnet(it2pair, ipair, rrtemp, temp, lowest)
                    numstn = it2pair
                  endif  
                enddo  
          
                if(numstn .gt. 0) then
                  rrstn(1) = ntstn(itarg)
                  rrstn(2) = 'Merged'
c                  print *,' 2:ipair,temp: ', ipair,temp(ipair, 1930, 13)
                  call statsubs(itarg, 1, 0, -1, idunit, rrtemp, rrflg,
     *              rrstn)
                else
                  print *,' Merged net does not meet minimum req.'
                endif

   60         enddo   
            endif  ! End of if method 1, 2 and/or 3
            
c >>>>      Write the "confirmfilt" journal file for restart
c            if(rentry .eq. 0 .and. unique .ne. '' .and.irandom .eq. 0)
C     *        call write_restart(2)
c >>>>      CRITICAL point for restart - come here if rentry == 2
 8080       if(rentry .eq. 2) print *,'---Rentry point 2---'
            rentry = 0
            
c           finished with statsubs - do not want method == 3 going
c             through confirmfilt (28sep05)
            if(method .eq. 3) goto 8085
 
c           !!!!!!!!!!!!!!!!!! test 5c2 - NO SHAP !!!!!!!!!!!!!!!!!!!
c           !!!!!!!!!!! reinstated SHAP at UCPdriver.v1 !!!!!!!!!!!!!
c           ------------------ Filter nfound into nhits ---------------
c           Revived ucpdriver.v1                             10may04

c           As of UCPM5 - confirmfilt ALSO:
c           read in history, merge with undoc and loop back around...
            call confirmfilt(ntstn,iexact,ifound,mmunit,mounit)
            
c           estamt uses minbic to calculate offset and QED t/ftests
c             no need for method 3                          06may06
            goto 8085

c           removed to retest t-test - from a comparison of v6c and v5, 
c              this is the only line to change!               26may05
c           no need to go back through all this - est. amp. & finish
c            method = 3
        
c           stop iterating if change points have not changed...
            if(iexact .eq. 1 .and. ifound .eq. 0) then
              print *,' No Changepoints or same AS LAST ITERATION ',
     *          '- SKIPPING LOWER LEVEL STNHIST'
              goto 95
            endif
              
   95     enddo ! end of (nht .ne. lht) loop - finding network chgpts
        
  100   endif ! end of NETLOOP option if
  
c >>>>  Write the "Final stnhist" journal file for restart
 8085   continue
c 8085   if(rentry.eq.0 .and. unique .ne. '' .and. irandom.eq.0)
c        if(rentry.eq.0 .and. unique .ne. '')
c     *    call write_restart(3)
c >>>>  CRITICAL point for restart - come here if rentry == 3
 8090   rentry = 0

c       using the nhits output of the first confirmfilt call, estimate
c         the amplitude of the changepoints
        call estamt(ntstn,idunit)
        
c        write out the results - 
c     BEWARE !!!!!! ONLY USED WITH read_write.mthly.v2a.f !!!!!!
c          call writsta(idunit, ipunit, orig, tflg, ntstn, nhits,
c     *      nindx)
        
c     end of all networks
c      print *,' Stopping Test'
c      goto 999
  120 enddo
  
c     close all the debug files for this station
      call closeunits()
c     close plotting data unit
      close(ipunit)
        
      go to 999

c     -----------------------------------------------------------
c     input/output errors
  200 print *,' Cannot open candidate-network file:', netfile
      stop 1

  210 print *,' Cannot open Matt data file:', mattdata
      stop 1

  215 print *,' Cannot open Matt meta file:', mattmeta
      stop 1

  999 end
     
     
c     =======================================================================
      subroutine statsubs(itarg, it2pair, ipair, ip2targ, idunit,
     *   rrtemp, rrflg, rrstn)
      INCLUDE 'inhomog.parm.mthly.incl'  
      INCLUDE 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'

c     version 2 has the 
c     Funnel the proper array subsets to the appropriate statistical
c       in the MLRTest application. Then assign the detected changepoint
c       info to the proper array.
c     parameters - 
c       itarg - station list index of the target station
c       it2pair - target's neighbor index to paired station
c       ipair - station list index of the paired station
c       ip2targ - pair's neighbor index to target station
c       idunit - output unit number for writing station data
c       rrtemp - cand & net temperature data
c       rrstn - network station list (base station first)
      character*6 rrstn(nstns)
      
c     Raw monthly temperature data array (stn/yr/molowest)
      real  rrtemp(2,nmo),rTraw(nmo)
      character*1 rrflg(nmo), tchar
      integer ichanges(nmo)

c     Returned discontinuity year/month, est. offset, model type,
c      and number inhomog from stats indate includes the month
c      of discontinuity (see imo2iym subroutine)
      integer indate1(ninh, ntech+1), indate2(ninh, ntech+1), 
     *  ndis(ntech+1), imtype(ninh, ntech+1), istnat(ninh)
      real discont(ninh, ntech+1), zscore(ninh, ntech+1)

      character*512 outstr, tmpstr, outcon, tmpcon
      real rmuq(2), rslpq(2), rsseq
      integer mknt2(2)

      integer ichgtyp(ninh), ichgimo(ninh)
      integer ichgintr(ninh), begchg, endchg
      real qx(nmo), qy(nmo)
            
      real j1intr(ninh), tintr(ninh), iqintr(ninh), imintr(ninh)

c     For the melding of the MLR into the IHT the following values are
c     used from the inhomog.parm.mthly.incl header file
c 
c     amiss - missing {real}
c     begyr - begin year 
c     endyr - end year
c     nstns - number of stations (including candidate)
        
      iedebug = 0
 
c     set writing for interim series for each test (YES=1;NO=0)
      iwrt = 0
      
      do it = 1, ntech+1
        ndis(it) = 0
      enddo 
      
      do imo = 1, nmo
        ichanges(imo) = 0
      enddo

c      !!!!   BEWARE: If any of the stat calls are changed - they must !!!!
c      !!!!    be changed in the alternate NETLOOP condition as well.  !!!!
      if(method .eq. 1 .or. method .eq. 2) then 

c       run the techniques backward, most to least complicated so the
c         ichange(imo) value is the lowest order model
c    -------------     run the split/merge routine   -------------------
        do it = hitech, lotech, -1
          if(itech(it).eq.1) then
            call splitmerge(idunit, rrtemp, rrstn, it, istnat,
     *        indate1(1,it), indate2(1,it), imtype(1,it), discont(1,it),
     *        zscore(1,it), ndis(it), iwrt, itarg, it2pair, ipair, 
     *        ip2targ, rTraw, imob, imoe)
            do ichg = 1, ndis(it)
              ichanges(indate1(ichg,it)) = it
c              print *,' combine :', it, ichg, indate1(ichg,it)
            enddo
          endif  
        enddo
        
        ichg = 0
        npts = 0
        do imo = 1, nmo
          if(ichanges(imo) .gt. 0) then
            ichg = ichg + 1
            ichgtyp(ichg) = ichanges(imo)
            ichgimo(ichg) = imo
c            print *,' ichanges: ', ichg, ichanges(imo), imo
          endif
        enddo
        nchg = ichg

c       loop back through the combined change points -
c         use the same assumptions that the single test routines use 
c         1) the chgpt is still valid (not straight line)
c         2) same type as originally detected
c       also
c         3) if 2 or more of the chgpts are within MINLEN
c            a) if the chgpt est are of the same sign then
c                 test each singly with same endpoints & keep lowest BIC Q
c            b) if not same sign then
c               opt1:  remove all data (add to ndelete) over interim chgpts
c               opt2:  remove all data (add to ndelete) over both segments
c               in either case: retain earliest chgpt
        begchg = imob
        ichg = 1
        
        do imo = 1, nmo   
          if(rrtemp(1, imo) .gt. amiss+1. .and. 
     *      rrtemp(2, imo) .gt. amiss+1. ) then
            if(imob .eq. 0) imob = imo
            imoe = imo
            rTraw(imo) = rrtemp(1,imo) - rrtemp(2,imo)
          else
            rTraw(imo) = amiss  
          endif  
        enddo  

        do while(ichg .le. nchg)
c         initialize segment series and raw diff series
          do imo = 1, nmo
            qx(imo) = amiss
            qy(imo) = amiss
          enddo

c         initialize interim changepoints
          intr = 1
          ichgintr(intr) = ichgimo(ichg)
          ilen = 0

c         if any are within minlen, put into interim chgpt array
          do while(ilen .le. minlen .and. ichg .lt. nchg)
            ilen = 0
            do imo = ichgimo(ichg)+1, ichgimo(ichg+1)
              if(rTraw(imo) .gt. amiss + 1.0) ilen = ilen + 1
            enddo 
            if(ilen .le. minlen) then 
              ichg = ichg + 1
              intr = intr + 1
              ichgintr(intr) = ichgimo(ichg)
            endif  
          enddo
          if(ichg .lt. nchg) then
            endchg = ichgimo(ichg+1)
          else
            endchg = imoe
          endif    

          nintr = intr
          call imo2iym(iyb,imb,begchg)
          call imo2iym(iye,ime,endchg)
          
          srtemp1 = 0.0
          srtemp2 = 0.0
          s2rtemp1 = 0.0
          s2rtemp2 = 0.0
          rntemp = 0.0
        
c         setup the x/y segments for testing
          jmo = 0
          do imo = begchg, endchg
            if(rTraw(imo) .gt. amiss + 1.0) then
              if(jmo .eq. 0) jmobeg = imo
              jmo = imo - jmobeg + 1
c             data series for minbic
              qy(jmo) = rTraw(imo)
              qx(jmo) = imo

c             for the standardized offset
              srtemp1 = srtemp1 + rrtemp(1,imo)
              srtemp2 = srtemp2 + rrtemp(2,imo)
              s2rtemp1 = s2rtemp1 + rrtemp(1,imo)*rrtemp(1,imo)
              s2rtemp2 = s2rtemp2 + rrtemp(2,imo)*rrtemp(2,imo)
              rntemp = rntemp + 1
            endif
          enddo
          jmoend = jmo
       
c         calculate std for each station in window.... assume the
c           lower std as the one without a chgpt and use for chgpt zscore
          stdrtemp1 = sqrt((s2rtemp1 - 
     *      (srtemp1 * srtemp1 / rntemp)) / (rntemp - 1))
          stdrtemp2 = sqrt((s2rtemp2 - 
     *      (srtemp2 * srtemp2 / rntemp)) / (rntemp - 1))
          stdrtemp = stdrtemp1
          if(stdrtemp2 .lt. stdrtemp .and. stdrtemp2 .ne. 0.0) 
     *      stdrtemp = stdrtemp2
c          print *,'stdrtemp: ',stdrtemp1, stdrtemp2, rntemp, stdrtemp

c         run BIC on all of the interior points, keeping the best BIC Q chgpt
          qlow = 9999.
          izero = 0
          do intr = 1, nintr
            do imo = 1, jmoend
              if(qx(imo) .ne. amiss .and. qx(imo) .le. ichgintr(intr))
     *          jend1 = imo
            enddo  
            call minbic(1, qx, qy, jend1, jmoend, tcrit, qstat, qmin,
     *        tadj, rmuq, rslpq, rsseq, 0, iqtype, mknt2, ifail)
            if(iedebug .gt. 0) then
              call imo2iym(iyi,imi,int(qx(jend1)))
              write(6,'("Interim chgpt: ",a,"-",a,2(i5,i2.2," - "),
     *          i5,i2.2,4f8.2,2f7.3,3i5))')
     *          rrstn(1), rrstn(2), iyb, imb, iyi, imi, iye,
     *          ime, qmin, qstat, tcrit, tadj, rslpq, iqtype, mknt2
            endif
            tintr(intr) = tadj
            iqintr(intr) = iqtype
            j1intr(intr) = jend1
            imintr(intr) = int(qx(jend1))
            if(qmin .lt. qlow) then
              minintr = intr
              minimo = int(qx(jend1))
              mintyp = iqtype
              qlow = qmin
              qadj = tadj
              qsse = rsseq
              tclow = tcrit
              qslow = qstat
            endif
c           see if all of the adjustments have the same sign (idomain = 1)
c             or not (idomain = 0)
            if(intr .eq. 1) then
              idomain = 1
              qsign = tadj
            else
c             ucpmonthly.v16a
c             changed "lt" to "le" - does this fix the Interim straight line?
c             ucpmonthly.v16b changed "le" back to "lt" 
              if(qsign * tadj .lt. 0.0) idomain = 0
            endif  
            if(tadj .eq. 0.0) then
              izero = izero + 1
            endif  
   10     enddo
          
c         if the interim chgpts are all the same sign, pick best
          if(idomain .eq. 1) then
            if(iedebug .gt. 0) then
              call imo2iym(iyi,imi,minimo)
              write(6,'("Same domain - Lowest Interim : ", a, "-", a,
     *          2(i5, i2.2, " - "), i5, i2.2, 4f8.2, 2i5))')
     *          rrstn(1), rrstn(2), iyb, imb, iyi, imi,
     *          iye, ime, qlow, qslow, tclow, qadj, mintyp, minintr
            endif
          else
c           if all of the interim chgpts came back as straight lines
c           remove chgpts - no need to remove data nor inc ndelete
            if(izero .eq. nintr) then
              call imo2iym(iyi,imi,minimo)
              write(6,'("Null domain - Earliest Interim : ", a, "-", a,
     *          2(i5, i2.2, " - "), i5, i2.2, 4f8.2, 2i5))')
     *          rrstn(1), rrstn(2), iyb, imb, iyi, imi,
     *          iye, ime, qmin, qstat, tcrit, qadj, mintyp, minintr
            else
c             if chgpts are not same sign, assume unusable data - remove
c               OPTION 1: data between interim chgpts.
              jend1 = j1intr(1)
              jmiss = jend1
c              print *,' set ndelete: ', itarg, ipair,
c     *          imintr(1)+1, ' to ', imintr(nintr)
              do imo = imintr(1)+1, imintr(nintr)
                ndelete(itarg,it2pair,imo) = 'D'
                rTraw(imo) = amiss
                jmiss = jmiss + 1
                qx(jmiss) = amiss
                qy(jmiss) = amiss
              enddo
c             segment evaluation has changed, go back through MINBIC for finals
              call minbic(1, qx, qy, jend1, jmoend, tcrit, qstat, qmin,
     *          tadj, rmuq, rslpq, rsseq, 0, iqtype, mknt2, ifail)
              minintr = 1
              qadj = tadj
              minimo = imintr(1)
              mintyp = iqtype
              if(iedebug .gt. 0) then
                call imo2iym(iyi,imi,minimo)
                write(6,'("Diff domain - Earliest Interim : ", a, "-",
     *             a,2(i5, i2.2, " - "), i5, i2.2, 4f8.2, 2i5))')
     *            rrstn(1), rrstn(2), iyb, imb, iyi, imi,
     *            iye, ime, qmin, qstat, tcrit, qadj, mintyp, minintr
              endif
            endif  
          endif

          if(mintyp .ge. 3) then
            do imo = minimo, nmo
c             Save the changepoint for the confirm filters
c             nfound is the info transport array for chgpt model type
              nfound(itarg,it2pair,imo) = achar(mintyp)
              nspan(itarg,it2pair,imo) = achar(ichg)
c             print *,' IN: ', mchg, ':', itarg,it2pair,imo

c             the changpoint is POSITIVE for the first sta of pair
              schgpt(itarg,imo) = schgpt(itarg,imo) + qadj
              zadj = abs(qadj / stdrtemp)
              zchgpt(itarg,imo) = zchgpt(itarg,imo) + zadj
              nchgpt(itarg,imo) = nchgpt(itarg,imo) + 1
c              print *,'1st:',itarg,imo,qadj,zadj,schgpt(itarg,imo),
c     *          zchgpt(itarg,imo),nchgpt(itarg,imo)

c             if the indx2 equal 0 then working with merged network
c              if(ip2targ .gt. 0) then
c             the changpoint is NEGATIVE for the second sta of pair
              schgpt(ipair,imo) = schgpt(ipair,imo) - qadj
              zchgpt(ipair,imo) = zchgpt(ipair,imo) + zadj
              nchgpt(ipair,imo) = nchgpt(ipair,imo) + 1
c              print *,'2nd:',ipair,imo,-qadj,zadj,schgpt(ipair,imo),
c     *          zchgpt(ipair,imo),nchgpt(ipair,imo)
c              endif
              if(rTraw(imo+1) .gt. amiss + 1.0) goto 20
            enddo
          else
c           for SLR push back begin of segment
c            minimo = begchg - 1  
c           DO NOT push back begin of segment
          endif  
          
c         set beginning of next test segment to best interim chgpt
  20      begchg = minimo + 1
          ichg = ichg + 1
        enddo
      else

c       set writing series to the index of the candidate (YES>0;NO=0)
        iwrt = itarg
        it = indKW
      
c       MM's Z-score Ref with T-test (the "known homogeneity test")
c       With metadata replaces SHAP - only executes the "merge" segtest
        if(itech(indKW).eq.1) call splitmerge(idunit, rrtemp, rrstn,
     *    indKW, istnat, indate1(1,indKW), indate2(1,indKW),
     *    imtype(1,indKW), discont(1,indKW), zscore(1,indKW), 
     *    ndis(indKW), iwrt, itarg, it2pair, ipair, ip2targ, rTraw, 
     *    imob, imoe)

c       By this stage, the Undoc Chgpt "hits" have been assigned to a SHF
c       event or coalesed into a solid Undoc Chgpt. The SHF events without
c       previously defined Undoc Chgpt have been appended and the output
c       includes only those that are significant with the "a priori" test.
c       So, at this point, only the validity of the unattached SHF events
c       are needed. However, all of the chgpts are checked, just for 
c       completeness.

          iocon = 0
          iostr = 0
          do in = 1, ndis(it)
            ischg = istnat(in)
            if(ischg .eq. 1) then
              isindx = itarg
            else
              isindx = ipair
            endif  
            x = discont(in,it)
c           changepoint is negative for the paired station
            if(ischg .eq. 2) x = -1. * x
            z = zscore(in,it)
            if(ischg .eq. 0) then
              if(idebug .ge. 2) print *,' WHAT!!!! ', itarg, it2pair,
     *          ipair, nchg1, nchg2, indate1(in,it), indate2(in,it)
              stop 99
            else if(ischg .gt. 2) then
c             do not use if there is more than one changepoint from one
c             or both stations within the gap
              if(idebug .ge. 2) then
                write(6,'("BIC confused: ",a,"-",a," ",a,i5,i3," ",
     *            " stn",i1,2(i5,i2.2,i5),2f7.2)') rrstn(1),rrstn(2),
     *            c2tech(it), itarg, it2pair, ischg, iy1, im1,
     *            indate1(in,it), iy2, im2, indate2(in,it), x, z
              endif  
            else
              call imo2iym(iy1,im1,indate1(in,it))
              call imo2iym(iy2,im2,indate2(in,it))
c             for each discontinuity, determine which station/chgpt to 
c             attribute it to, or whether the discont can be used at all
              if(idebug .ge. 2) then
                write(6,'("BIC chngpts: ",a,"-",a," ",a,i5," ",
     *            " stn",i1,2(i5,i2.2,i5),2f7.2)') rrstn(1),rrstn(2),
     *            c2tech(it), isindx, ischg, iy1, im1,
     *            indate1(in,it), iy2, im2, indate2(in,it), x, z
              endif
              do imo = indate1(in,it), indate2(in,it)
c               accumulate BIC chngpts for amplitude estimation (estamt)
                schgpt(isindx,imo) = schgpt(isindx,imo) + x
                zchgpt(isindx,imo) = zchgpt(isindx,imo) + z
                nchgpt(isindx,imo) = nchgpt(isindx,imo) + 1
              enddo
              print *,'3rd:',isindx,indate1(in,it),'-',indate2(in,it),
     *          x,z,schgpt(isindx,imo),zchgpt(isindx,imo),
     *          nchgpt(isindx,imo)
            endif ! end if one and only one chgpt check
   50     enddo ! end of chgpt returned from BIC loop

      endif ! end if method == 1/2 and else

  999 return
      end

c     =======================================================================
      subroutine setrflg(iorig, rrflg, origflg)
c     iorig - index in the "original" flag array 
c     rrflg - round robin flag array
c     origflg - original flag array
      INCLUDE 'inhomog.parm.mthly.incl'  
      INCLUDE 'inhomog.comm.mthly.incl'  
      character*1 origflg(maxstns,begyr:endyr,13), rrflg(nmo)
      integer iorig
      do iy = begyr+1, endyr
c       if annual time resolution (lomth,himth=13)
        do im = lomth, himth
          call iym2imo(iy,im,imo)
          rrflg(imo) = origflg(iorig, iy, im)
        enddo
      enddo
      return
      end
          
c     =======================================================================
      subroutine setrnet(irr, iorig, rrdata, origdata, lowest)
c     irr - index in the round robin array
c     iorig - index in the "original" flag array 
c     rrdata - round robin data array
c     origdata - original data array
c
c     include removal of series monthly average to make monthly anomaly series
      INCLUDE 'inhomog.parm.mthly.incl'  
      INCLUDE 'inhomog.comm.mthly.incl'  
      real origdata(maxstns,begyr:endyr,13), 
     *  rrdata(2,nmo),mthavg(13)
      integer iorig, irr, lowest
      lowest = 0
c     if annual time resolution (lomth,himth=13)
      do im = lomth, himth
        sum = 0.0
        num = 0
        do iy = begyr+1, endyr
          call iym2imo(iy,im,imo)
          if(origdata(iorig,iy,im) .gt. amiss+1) then
            sum = sum + origdata(iorig, iy, im)
            num = num + 1
          endif  
        enddo
        if(num .gt. minann) then 
          mthavg(im) = sum / num
        else
          mthavg(im) = amiss
        endif
      enddo
c      print *,'Mthavg: ', mthavg

      do iy = begyr+1, endyr
        do im = lomth, himth
          call iym2imo(iy,im,imo)
          if(origdata(iorig, iy, im) .gt. amiss+1. .and. 
     *      mthavg(im) .gt. amiss+1.) then
            rrdata(irr, imo) = origdata(iorig, iy, im) - mthavg(im)
            lowest = lowest + 1
          else
            rrdata(irr, imo) = amiss
          endif
        enddo
      enddo

      return
      end
          
c     =======================================================================
      subroutine inirnet(irr, rrdata)
c     irr - index in the round robin array
c     rrdata - round robin data array
      INCLUDE 'inhomog.parm.mthly.incl'  
      real rrdata(2,nmo)
      integer irr
      do im = 1, nmo
        rrdata(irr, imo) = amiss
      enddo
      return
      end
          
c     =======================================================================
      subroutine confirmfilt(ntstn,iexact,ifound,mmunit,mounit)

c     ntstn - station names
c     nfound - "hits" from the current technique
c     nhits - 2ofX filtered "hits" from all of the techniques

      INCLUDE 'inhomog.parm.mthly.incl'  
      include 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'
c     number of model indices in itypshow (used for model types = 3 to 8)
      parameter (ntyp = 6)
      character*6 ntstn(maxstns)
      character*6 procstr
      character*500 outstr
      character*32 tmpstr
      integer jhits(maxstns, nmo), khits(nmo),
     *  shhits(maxstns,nmo), ktests(nmo)
      real amps(maxstns, nmo), akhigh(nmo)
      integer ifndshow(nmo,maxstns), itypshow(nmo,maxstns,ntyp)

      maxmo = nmo
      if(itimeres .eq. 0) maxmo = numyr
      
c      icdebug = 2
      icdebug = 10
      numtech = (hitech - lotech) + 1

C     ---------------------------------------------------------------------
c     filter years and separate into individual stations 
c       with the "2-of-X" logic
      do itarg = 1, maxstns
        lastimo = 0
        lasthits = 0
        do imo = 1, nmo
c         keep track of the last nhit
          if(nhits(itarg,imo) .ne. 0) then
            lastimo = imo
            lasthits = nhits(itarg,imo)
            nhits(itarg,imo) = 0
          endif  
        enddo
        if(lastimo .gt. 0) then
          nhits(itarg,lastimo) = lasthits
          print *,'lastimo: ',itarg,lastimo
        endif  
      enddo 

      indx = 0
      do imo = 1, nmo
        indx = indx + 1
        do k = 1, maxstns
          do it = 1, ntech+1
            jhits(k, indx) = 0
            amps(k, indx) = 0.0
          enddo
        enddo
      enddo        
        
c      print *,'          itarg     it2pair  iy1  nfound'
c     follow all the "found" chgpts - if exactly equal to last
c       loop then increment the method
      ifound = 0
      iexact = 0
c     initialize show array
      do k = 1, maxstns
        if(ntr(k) .eq. 0) goto 10
        itarg = nindx(k,1)
        do imo = 1, nmo
          ifndshow(imo,itarg) = 0
          do ityp = 1, ntyp
            itypshow(imo,itarg,ityp) = 0
          enddo  
        enddo  
   10 enddo 
  
      print *,' Ndelete array: ', ifound
      if(icdebug .lt. 2) goto 610
      print *,' Exact wrt last: ', iexact
      do k1 = 1, numsubs, 100
        k2 = k1 + 99
        if(k2 .gt. maxstns) k2 = maxstns
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 590
          itarg = nindx(k,1)
          outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(1:3)// '|'
  590   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 595
          itarg = nindx(k,1)
          outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(4:6)// '|'
  595   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        nofound = 0
        do imo = 1, nmo
          call imo2iym(iy,im,imo)
          write(outstr, '(i4,i3,i5," |")') iy,im,imo
          mofound = 0
          do k = k1, k2
            if(ntr(k) .eq. 0) goto 600
            itarg = nindx(k,1)
            if(temp(itarg,iy,im) .eq. amiss) then
              outstr = outstr(1:lnblnk(outstr)) // '-X-'
            else
              idel = 0
              do ipair = 1, nstns
                if(ndelete(itarg,ipair,imo) .ne. czero) idel = idel + 1
              enddo
              if(idel .eq. 0) then
                outstr = outstr(1:lnblnk(outstr)) // '---'
              else
                write(tmpstr,'(i3)') idel
                outstr = outstr(1:lnblnk(outstr)) // tmpstr
                mofound = mofound + 1
              endif  
            endif  
            outstr = outstr(1:lnblnk(outstr)) // '|'
          enddo  
  600     if(mofound .gt. 0) then
            if(nofound .gt. 0) then
c              write(6,*) nofound
              nofound = 0
            endif  
            write(6,*) outstr(1:lnblnk(outstr))
          else
            nofound = nofound + 1  
          endif  
        enddo
      enddo  

c     ------------------ pre-FILTER 1 ----------------------
c     for each station and its sub-net in the network sum
c       from the paired itfound array:
c       1) all the non-zero hits into the ifndshow array
c       2) the model types into the itypshow array
  610 do k = 1, numsubs
        if(ntr(k) .eq. 0) goto 30
        itarg = nindx(k,1)
        do imo = 1, nmo
c         filter found by the paired confirm hits per technique
c         added ucpdriver.v1                           10may04
          itfound = 0
c         changed begin loop from itarg+1 to itarg, 14 Aug 03
          do it2pair = 1, nstns
            ipair = nindx(itarg, it2pair)
            if(ipair .eq. 0) goto 20
c           convert chgpt model byte to integer and remove SLR codes
            mtype = iachar(nfound(itarg,it2pair,imo)) - 2
c     -- CHANGES BY DANIEL ROTHENBERG <darothen@mit.edu> -- BEGIN
c        -- make sure mtype doesn't break things. The above call to 
c        -- iachar will retrive a '\' because the default 'nfound' value
c        -- was set to czero = '\'. iachar('\0') yields 90, so mtype above
c        -- could be set to 90, which will break the lookup in 'itypshow'
c        -- below, whose third dimension is only size 6.
            if(mtype .gt. 0 .and. mtype .ne. 90) then
c     -- CHANGES -- END
              if(icdebug .ge. 3)
     *          print *,'mtype: ',itarg,it2pair,imo,it,mtype
              itfound = itfound + 1
              itypshow(imo,itarg,mtype) = itypshow(imo,itarg,mtype)+1
            endif  
          enddo  
c         filter found by the number of separate pairs hit
   20     ifndshow(imo,itarg) = itfound
          if(itfound .gt. 0) then
            call imo2iym(iy,im,imo)
            if(icdebug .ge. 3)
     *        print *,'ifndshow:',ifndshow(imo,itarg),
     *          imo,iy,im,itarg
          endif
   25   enddo
   30 enddo
  

      if(icdebug .lt. 1) goto 38
c     print out the current working output array
      print *,' Exact wrt last: ', iexact
      do k1 = 1, numsubs, 100
        k2 = k1 + 99
        if(k2 .gt. maxstns) k2 = maxstns
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 32
          itarg = nindx(k,1)
          if(numtech .eq. 1) then
            outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(6:6)// '|'
          else
            outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(1:3)// '|'
          endif  
   32   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 34
          itarg = nindx(k,1)
          if(numtech .eq. 1) then
            outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(6:6)// '|'
          else
            outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(4:6)// '|'
          endif  
   34   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        nofound = 0
        do imo = 1, nmo
          call imo2iym(iy,im,imo)
          write(outstr, '(i4,i3,i5," |")') iy,im,imo
          mofound = 0
          do k = k1, k2
            if(ntr(k) .eq. 0) goto 36
            itarg = nindx(k,1)
            if(ifndshow(imo,itarg) .lt. 1) then
              outstr = outstr(1:lnblnk(outstr)) // '---'
            else
              write(tmpstr,'(i3)') ifndshow(imo,itarg)
              outstr = outstr(1:lnblnk(outstr)) // tmpstr
              mofound = mofound + 1
            endif  
            outstr = outstr(1:lnblnk(outstr)) // '|'
          enddo  
   36     if(mofound .gt. 0) then
            if(nofound .gt. 0) then
c              write(6,*) nofound
              nofound = 0
            endif  
            write(6,*) outstr(1:lnblnk(outstr))
          else
            nofound = nofound + 1  
          endif  
        enddo
      enddo  

      if(icdebug .lt. 1) goto 38
c     print out the current working output array
      do k1 = 1, numsubs, 100
        k2 = k1 + 99
        if(k2 .gt. maxstns) k2 = maxstns
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 320
          itarg = nindx(k,1)
          outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)// '|'
  320   enddo
        write(6,*) outstr(1:lnblnk(outstr))
        write(outstr, '("             |")')
        nofound = 0
        do imo = 1, nmo
          call imo2iym(iy,im,imo)
          write(outstr, '(i4,i3,i5," |")') iy,im,imo
          mofound = 0
          do k = k1, k2
            if(ntr(k) .eq. 0) goto 360
            itarg = nindx(k,1)
            do ityp = 1, ntyp
              if(itypshow(imo,itarg,ityp) .lt. 1) then
                outstr = outstr(1:lnblnk(outstr)) // '-'
              else if(itypshow(imo,itarg,ityp) .gt. 9) then
                outstr = outstr(1:lnblnk(outstr)) // '#'
                mofound = mofound + 1
              else
                write(tmpstr,'(i1)') itypshow(imo,itarg,ityp)
                outstr = outstr(1:lnblnk(outstr)) // tmpstr
                mofound = mofound + 1
              endif  
            enddo
            outstr = outstr(1:lnblnk(outstr)) // '|'
          enddo  
  360     if(mofound .gt. 0) then
            if(nofound .gt. 0) then
c              write(6,*) nofound
              nofound = 0
            endif  
            write(6,*) outstr(1:lnblnk(outstr))
          else
            nofound = nofound + 1  
          endif  
        enddo
      enddo  

C     ---------------------------------------------------------
c     Expanded from one-dimensional (single yr-mth for all stn/pairs)
c     to a multi-dimensional filter to clear spans of missing data
c     as implemented in the statsubs routine. These spans vary between
c     each paired comparison.
c                                                     Aug 05 cw
c     test the "higher count erasure" technique - 
c     From the paired hit array "nfound" containing "1" or "0" hits for each
c     target-pair-yr/mth-stat: final output "jhits" & "amps" contain the
c     total hits & est adj. parsed and attributed to the "Best Guess" stations
   38 if(icdebug .gt. 0) then
        print *,'  NET  STN    FILT TECH    YEAR MTH          ',
     *    '  AVG   STD NUM ITFOUND'
      endif
c     ... find the highest chgpt hit occurance for all targets/all yr-mths
      istop = 0
      do while (istop .eq. 0) 
        ihighit = 0
c       ... for each year/month
        do imo = 1, nmo
c         ... for each target
          do k = 1, numsubs
            if(ntr(k) .eq. 0) goto 40
            itarg = nindx(k,1)
            if (ifndshow(imo,itarg) .gt. ihighit) then
              ihighit = ifndshow(imo,itarg)
              ihighnet = itarg
              ihighmo = imo
            endif
   40     enddo ! end target search loop
        enddo ! end yr-mth search loop

c       keep looping until ihighit == 1!!!!
        if(ihighit .ge. 2) then
          if(icdebug .ge. 4) then
            call imo2iym(iy,im,ihighmo) 
            print *,' Yr/mth: ',ihighmo,iy,im,' ihighnet: ',ihighnet,
     *        ' ihighit: ', ihighit
          endif
c         set the current highest hit station for this y/m/t
          itarg = ihighnet
c         get the highest ifndshow out of the way...
c           ... go thru the pairs
          do it2pair = 2, nstns
            ipair = nindx(itarg, it2pair)
            if(ipair .eq. 0) goto 45
c           erase the hits for the event with the pair
            mtype = iachar(nfound(itarg,it2pair,ihighmo)) - 2
            if(mtype .gt. 0) then
c              nfound(itarg,it2pair,ihighmo) = czero
c             fetch the chgpt number of the event
              nchg = iachar(nspan(itarg,it2pair,ihighmo))
              if(nchg .eq. 0) then
                print *,'1Nfound > 0;Nspan == 0!!!!!', k, itarg, it2pair,
     *            ipair, ihighmo,nchg,' ',mtype
                stop 1
              endif  
c             go backward ... clearing out peripherals
              do imo = ihighmo-1, 1, -1
                ichg = iachar(nspan(itarg,it2pair,imo))
                if(ichg .ne. nchg) goto 42
                nfound(itarg,it2pair,imo) = czero
                nspan(itarg,it2pair,imo) = czero
                ifndshow(imo,itarg) = ifndshow(imo,itarg)-1
                ntest(itarg,imo) = ntest(itarg,imo) - 1
              enddo
c             go forward ... clearing out peripherals
   42         do imo = ihighmo+1, nmo
                ichg = iachar(nspan(itarg,it2pair,imo))
                if(ichg .ne. nchg) goto 44
                nfound(itarg,it2pair,imo) = czero
                nspan(itarg,it2pair,imo) = czero
                ifndshow(imo,itarg) = ifndshow(imo,itarg)-1
                ntest(itarg,imo) = ntest(itarg,imo) - 1
              enddo
            endif  
   44     enddo ! end of paired station loop for target

   45     ifndshow(ihighmo,itarg) = -1 * ifndshow(ihighmo,itarg)
c         go through all of the stations, skipping the highest
c           (chosen) target
          do it2pair = 1, numsubs
            ipair = nindx(it2pair, 1)
            if(ipair .eq. itarg) goto 50
c           go through the neighbors of each paired station
            do ip2targ = 2, nstns
              js = nindx(ipair, ip2targ)
c             if one of the paired station neighbors is the chosen target
              if(js .eq. itarg) then
c               zero out the chgpt for that pair!
                mtype = iachar(nfound(ipair,ip2targ,ihighmo)) - 2
                if(mtype .gt. 0) then
c                 fetch the chgpt number of the event
                  nchg = iachar(nspan(ipair,ip2targ,ihighmo))
                  if(nchg .eq. 0) then
                    print *,'2Nfound > 0;Nspan == 0!!!!!', k, itarg, 
     *                it2pair, ipair, ihighmo,nchg,' ',mtype
                    stop 2
                  endif  
c                 go backward ... clearing out peripherals
                  do imo = ihighmo-1, 1, -1
                    ichg = iachar(nspan(ipair,ip2targ,imo))
                    if(ichg .ne. nchg) goto 47
                    nfound(ipair,ip2targ,imo) = czero
                    nspan(ipair,ip2targ,imo) = czero
                    ifndshow(imo,ipair) = ifndshow(imo,ipair)-1
                    ntest(ipair,imo) = ntest(ipair,imo) - 1
                  enddo
   47             imo1 = imo + 1
c                 go forward (including ihighmo)... clearing out peripherals
                  do imo = ihighmo, nmo
                    ichg = iachar(nspan(ipair,ip2targ,imo))
                    if(ichg .ne. nchg) goto 49
                    nfound(ipair,ip2targ,imo) = czero
                    nspan(ipair,ip2targ,imo) = czero
                    ifndshow(imo,ipair) = ifndshow(imo,ipair)-1
                    ntest(ipair,imo) = ntest(ipair,imo) - 1
                  enddo
   49             imo2 = imo - 1
                  if(icdebug .ge. 4) 
     *              print *,'Zeroed STN: ',ntstn(itarg),' neigh:', 
     *              ntstn(ipair),' neigh index:',ip2targ,' range:',imo1,
     *              imo2,' fndshow:',ifndshow(imo,ipair)
                endif
                goto 50
              endif  
            enddo ! end neighbors of pair loop
            if(icdebug .ge. 4) 
     *       print *,' No link from ',ntstn(ipair),' back to ', 
     *        ntstn(itarg)
   50     enddo ! end of search for pairs to target
        else
          istop = 1  
        endif
      enddo ! end of the highest occurance search for this tech

c     all of the peripheral "hits" in the nfound array associated 
c       with each filtered changepoint should be zeroed out by this
c       point, leaving only the attributed chgpt-station data
c     go thru the year/months again.....
      do imo = 1, nmo
c       for each target ...
        do k = 1, numsubs
          if(ntr(k) .eq. 0) goto 70
          itarg = nindx(k,1)
          itfound = 0
c         ... and its sub-net in the network          
          do it2pair = 1, nstns
            ipair = nindx(itarg, it2pair)
            if(ipair .eq. 0) goto 60
c           accumulate the remaining hits
            if(nfound(itarg,it2pair,imo).ne.czero) itfound = itfound + 1
          enddo ! end of paired station loop
              
c         filter found by the number of separate pairs hit
c           using ifndshow for array printout
   60     ifndshow(imo,itarg) = itfound
          if(itfound .ge. iconfirm) then
c           generate adjustment around output (iconfirm > 2 assumed)
            sum = schgpt(itarg,imo)
            sum2 = zchgpt(itarg,imo)
            num = nchgpt(itarg,imo)
c           save for the next filter
            jhits(itarg,imo)=num
            avg = sum / num
c           zchgpt has the z-scores for the size of hit
            avgz = sum2 / num
            amps(itarg,imo) = avgz
            if(icdebug .gt. 0) then
              call imo2iym(iy,im,imo)
              procstr = 'CONFRM'
              write(6,1000) itarg,ntstn(itarg),procstr,
     *          nloop,imo, iy, im, avg, avgz, num, itfound
              ifound = ifound + 1
 1000         format(i5,1x,a6,'-',a6,' MW',i1,' at ',2i5,i3,
     *          ' AVG ADJ: ',2f6.2,2i4)
            endif
          else
            ifndshow(imo,itarg) = 0  
          endif  
   70   enddo ! end of target loop
      enddo ! end of yr/mth loop
      print *,' Confirm filter: ',ifound

      if(icdebug .lt. 2) goto 100
c     print out the current working output array
      print *,' Exact wrt last: ', iexact
      do k1 = 1, numsubs, 100
        k2 = k1 + 99
        if(k2 .gt. maxstns) k2 = maxstns
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 80
          itarg = nindx(k,1)
          if(numtech .eq. 1) then
            outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(6:6)// '|'
          else
            outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(1:3)// '|'
          endif  
   80   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 85
          itarg = nindx(k,1)
          if(numtech .eq. 1) then
            outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(6:6)// '|'
          else
            outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(4:6)// '|'
          endif  
   85   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        nofound = 0
        do imo = 1, nmo
          call imo2iym(iy,im,imo)
          write(outstr, '(i4,i3,i5," |")') iy,im,imo
          mofound = 0
          do k = k1, k2
            if(ntr(k) .eq. 0) goto 90
            itarg = nindx(k,1)
              if(ifndshow(imo,itarg) .lt. iconfirm) then
                outstr = outstr(1:lnblnk(outstr)) // '---'
              else
                write(tmpstr,'(i3)') ifndshow(imo,itarg)
                outstr = outstr(1:lnblnk(outstr)) // tmpstr
                mofound = mofound + 1
              endif  
            outstr = outstr(1:lnblnk(outstr)) // '|'
          enddo  
   90     if(mofound .gt. 0) then
            if(nofound .gt. 0) then
c              write(6,*) nofound
              nofound = 0
            endif  
            write(6,*) outstr(1:lnblnk(outstr))
          else
            nofound = nofound + 1  
          endif  
        enddo
      enddo  

C     ---------------------------------------------------------------------
c     now the new second filter.....
c     these are the rules:
c     TEST WITH THE STATION HISTORY DATA IN SAHIST
c     for each station
c       go through the year/months.....
c         for each technique...
c            using the est amplitude of each jhit - 
c              if any SHF is within +/- mrgyr then
c                merge jhit into SHF location
c         (all "jhits"/"amps" used in accum are zeroed out)
c         Output accum SHF hits in the "nhits"/"ntests" arrays

  100 do imo = 1, nmo
        do is = 1, numsubs
          itarg = nindx(is,1)
          shhits(itarg, imo) = 0
c          ntests(itarg, imo) = 0
        enddo
      enddo    
  
      do is = 1, numsubs
        itarg = nindx(is,1)
        if(itarg .eq. 0) goto 105
        do inh = 1, ninh
          if(sahist(is,inh) .eq. 0) goto 105
          imo = sahist(is,inh)
c         stn hist is locked into position by SHF (and UCP in filter 3)
c          print *,' SHHITS set:', itarg,imo
          shhits(itarg,imo) = 1
          nhits(itarg,imo) = 100
        enddo
  105 enddo        

c     for each target station...
      if(icdebug .ge. 2) then
        print *,'  NET  STN    FILT TECH     YR MTH JSUM  AVG   STD  ',
     *    'IMO JINC RNGE NHITS'
      endif
      ifound = 0
      do k = 1, numsubs
        if(ntr(k) .eq. 0) goto 120
        itarg = nindx(k,1)
c       average of the monthly standard deviations
        call stdmth(itarg,stdk)
c       ....for all the year/months
        do imo = 1,nmo
c         ....find the hits from statsubs
          if(jhits(itarg,imo) .gt. 0) then
c	      use amps array to define mth range
            ahigh = amps(itarg,imo)
            jsum = jhits(itarg,imo)
c           standardize the size of the offset
            astd = abs(ahigh)
            do irange = 1,nrange
	        if(astd .lt. arange(irange)) goto 110
	      enddo
  110       if(irange .gt. nrange) irange = nrange

c           make the search bracket twice the range plus the current yr/mth
            ibracket = mrgyr(irange)*2 + 1
            if(icdebug .ge. 3) then
              call imo2iym(iy,im,imo)
              print *,' ASTD: ',itarg,ntstn(itarg),nloop,iy,im,astd,
     *          ahigh, stdk, jsum, ibracket
            endif

c           go through the bracket of months, look for SHF record
c           start at the hit-point index and expand outward 
c           in the series 0, -1, 1, -2, 2, etc.
c           this will merge the hits with the closest SHF record
c           keep track of missing values in both directions
            jupmiss = 0
            jdnmiss = 0
            do iinc = 1, ibracket
c             jinc is the expansion index
              jinc = iinc/2
              if(mod(iinc,2).eq.1) then
  111           jmo = imo + jinc + jupmiss
                if(jmo .gt. maxmo) goto 114
                call imo2iym(iyj,imj,jmo)
                if(temp(itarg,iyj,imj) .lt. amiss+1.) then
c                 this month is missing, goto the next
                  jupmiss = jupmiss + 1
                  goto 111
                endif  
              else 
  112           jmo = imo - jinc - jdnmiss
                if(jmo .lt. 1) goto 114
                call imo2iym(iyj,imj,jmo)
                if(temp(itarg,iyj,imj) .lt. amiss+1.) then
c                 this month is missing, goto the next
                  jdnmiss = jdnmiss + 1
                  goto 112
                endif  
              endif

c             If there is an SHF record (at jmo) absorb the UCP hits (from imo)
              if(shhits(itarg,jmo).gt.0)then
                nhits(itarg,jmo) = nhits(itarg,jmo) + jsum
c                ntests(itarg,jmo) = ntests(itarg,jmo) + ntest(itarg,jmo)
c               zero test array block for next iter
c                print *,'Zero jhits:',itarg,imo
                jhits(itarg,imo) = 0
                amps(itarg,imo) = 0.0
c               let the world know
                if(icdebug .ge. 2) then
                  call imo2iym(iy,im,imo)
                  write(6,2000) itarg, ntstn(itarg),nloop,iy,im,jsum,
     *              ahigh,astd,imo,jmo,mrgyr(irange),
     *              nhits(itarg,jmo)
 2000             format(i5,1x,a6,'-SHFHIT KW',i1,' at ',i4,i3,i4,
     *              2f6.2,4i5)
                endif
c               found what we wanted, break out
                goto 115
              endif ! end of SHF record occurence
                
  114       enddo ! end of the bracket loop
  115       continue
          endif 
        enddo ! end of year/month loop
  120 enddo ! end of station loop
  
      print *,' Merge to SHF: ',ifound
      if(icdebug .lt. 2) goto 150
c      if(icdebug .lt. 3) goto 150
      print *,' Exact wrt last: ', iexact
      write(outstr, '("             |")')
      do k1 = 1, numsubs, 100
        k2 = k1 + 99
        if(k2 .gt. maxstns) k2 = maxstns
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 130
          itarg = nindx(k,1)
          outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(1:3)// '|'
  130   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 135
          itarg = nindx(k,1)
          outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(4:6)// '|'
  135   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        nofound = 0
        do imo = 1, nmo
          call imo2iym(iy,im,imo)
          write(outstr, '(i4,i3,i5," |")') iy,im,imo
          mofound = 0
          do k = k1, k2
            if(ntr(k) .eq. 0) goto 140
            itarg = nindx(k,1)
            if(nhits(itarg,imo) .eq. 0) then
              outstr = outstr(1:lnblnk(outstr)) // '---'
            else
              write(tmpstr,'(i3)') nhits(itarg,imo)
              outstr = outstr(1:lnblnk(outstr)) // tmpstr
              mofound = mofound + 1
            endif  
            outstr = outstr(1:lnblnk(outstr)) // '|'
          enddo  
  140     if(mofound .gt. 0) then
            if(nofound .gt. 0) then
c              write(6,*) nofound
              nofound = 0
            endif  
            write(6,*) outstr(1:lnblnk(outstr))
          else
            nofound = nofound + 1  
          endif  
        enddo
      enddo  

C     ---------------------------------------------------------------------
c     now the third filter.....(the old second filter)
c     THESE SHOULD BE UCP'S WITH NO SHF RECORDS
c     these are the rules:
c       Input is all the hits that are leftover from the SHF filter in the
c         "jhits" & "amps" arrays
c       go back through year/months.....
c         for each technique...
c           find the highest remaining hits (down to confirm)
c            accumulate all of the "jhits" & "ntests" for each month 
c              +/- mrgyr(from "amps") while skipping missing data
c         (all "jhits" & "amps" used in accum are zeroed out)
c         for all accumulations within nmrgyr 
c           when given month are greater or equal to ithres then add to 
c             the "nhits" & "ntests" arrays

c     for each station...
  150 if(icdebug .gt. 0) then
        print *,'  NET  STN    FILT TECH     YR MTH JSUM  AVG  STD ',
     *    'RNG PR THRES'
      endif
      ifound = 0
      do k = 1, numsubs
        if(ntr(k) .eq. 0) goto 180
        itarg = nindx(k,1)
c       initialize the standard deviation for the amps vs range xform
        call stdmth(itarg,stdk)
        
c       initialize interim temp arrays...
        do imo = 1, nmo
          khits(imo) = 0
          ktests(imo) = 0
          akhigh(imo) = 0.0
        enddo
                  
c       iterate until there are no more high points
        istop = 0
        do while (istop .eq. 0) 
        
          ihighit = 0
          ahigh = 0.0
c         ....find the "highest count"
          do imo = 1,nmo
c           ... sum all techniques
            isum = 0
            asum = 0.0
c           sum only those greater than iconfirm (:= 2)
            if(jhits(itarg,imo) .ge. iconfirm)then
              jhit = jhits(itarg,imo)
              isum = isum + jhit
              asum = asum + amps(itarg,imo) * jhit
            endif

c           ... find the highest chgpt hit station
c           not sure whether this should be by "hits" or "amps"....
c           this way is amps......
c            if (asum .gt. ahigh) then
c           this way is hits......
            if (isum .gt. ihighit) then
              ihighit = isum
              ihighmo = imo
              ahigh = asum / isum
            endif
          enddo ! end find the highest count loop
          
c         the highest hit value in the array is at IHIGHMO, the sum of
c           hits over all stat-tests is IHIGHIT, and the est. adj. is AHIGH

          if(idebug .ge. 4) 
     *     print *,'----itarg,ihighit,ihighmo,ahigh,stdk',
     *      itarg,ihighit,ihighmo,ahigh,stdk
     
c         keep going until there are no more hits
          if(ihighit .gt. 0) then
c           bracket the highest hit +/- nmrgyr (from command line)
            if(nmrgyr .ne. -2) then
              ibracket = nmrgyr * 2 + 1
            else
c	      else bracket using ampitude of chgpt to define mth range
              if(method .eq. 3) then
                irange = nrange
              else
c               base the range on the standardized amplitude est.
c               standardize the size of the offset
                astd = abs(ahigh)
                do irange = 1,nrange
                  if(astd .lt. arange(irange)) goto 160
                enddo
c                print *,' ASTD: ', astd, ahigh, stdk
              endif  
  160         if(irange .gt. nrange) irange = nrange
              ibracket = mrgyr(irange) * 2 + 1
            endif ! end of defining the search bracket

c           go through the bracket, look for the highest hits already found
c           start at the hit-point index and expand outward 
c           in the series 0, -1, 1, -2, 2, etc.
            ihighmo = ihighmo
c           keep track of missing values in both directions
            jupmiss = 0
            jdnmiss = 0
            do iinc = 1, ibracket
c             jinc is the expansion index
              jinc = iinc/2
              if(mod(iinc,2).eq.1) then
  161           jmo = ihighmo + jinc + jupmiss
                if(jmo .gt. maxmo) goto 164
                call imo2iym(iyj,imj,jmo)
                if(temp(itarg,iyj,imj) .lt. amiss+1.) then
c                 this month is missing, goto the next
                  jupmiss = jupmiss + 1
                  goto 161
                endif  
              else 
  162           jmo = ihighmo - jinc - jdnmiss
                if(jmo .lt. 1) goto 164
                call imo2iym(iyj,imj,jmo)
                if(temp(itarg,iyj,imj) .lt. amiss+1.) then
c                 this month is missing, goto the next
                  jdnmiss = jdnmiss + 1
                  goto 162
                endif  
              endif ! end of finding the next non-missing datum

c             absorb lesser hit into closest higher hit
              if(khits(jmo) .gt. 0) then
                khits(jmo) = khits(jmo) + ihighit
                akhigh(jmo) = akhigh(jmo) + ahigh*ihighit
                if(idebug .ge. 3)
     *            print *,' Absorb hit: ', itarg, ihighmo, ' to ', 
     *              jmo, khits(jmo), ktests(jmo), 
     *              akhigh(jmo)/khits(jmo)
                goto 165
              endif  
  164       enddo ! end of bracket loop
  
c           if no hits found, setup new hit
            khits(ihighmo) = ihighit
            ktests(ihighmo) = ntest(itarg,ihighmo)
            akhigh(ihighmo) = ahigh * ihighit
            if(idebug .ge. 3)
     *        print *,' New CHG hit: ', itarg, ihighmo, khits(ihighmo),
     *          ktests(ihighmo), akhigh(ihighmo)/khits(ihighmo)

c           zero test array block for next iter
  165       jhits(itarg,ihighmo) = 0
            amps(itarg,ihighmo) = 0.0
          else
c           no hits >= confirm: stop
            istop = 1
          endif ! end of working with this UCP
        enddo ! end of finding UCP for this station
            
c       examine the interim khits array for station's filtered chgpts
        do imo = 1, nmo
c         ... if highest hits >= ithres(npair) (from param study) then save
c         fetch the number of pairs tested
c         WARNING - THE NUMBER OF PAIRWISE TESTS ARE SUMMED FOR EACH
C           TECHNIQUE. FOR THE "PARMSET" VERSION THIS IS ASSUMED TO
C           BE THE TOTAL NUMBER OF STATISTICAL TECHNIQUES AVAILABLE.
C           IF THIS IS NOT THE CASE, NPAIR MUST BE ADJUSTED APPROPRIATELY!
C                                                      10 MARCH 2005 CW
          if(khits(imo) .gt. 0) then
            call imo2iym(iy,im,imo)
            npair = ktests(imo)
            jsum = khits(imo)
            if(npair .le. 5 .and. npair .lt. nstns - 2) then
              print *,'NT ',ntstn(itarg),itarg,imo,iy,im,
     *          ktests(imo),npair
            endif  
            do ihss = 1, nhss-1
              if(npair .le. nhtstns(ihss)) goto 170
            enddo  
  170       ihthres = nhthres(ihss)

c           TEST CODE - for FINSHF threshold test
            if(method .eq. 3) ihthres = 1                
            if(idebug .ge. 3) 
     *      print *,'itarg,imo,iy,im,npair,jsum,ihthres,stdk',
     *        itarg,imo,iy,im,npair,jsum,ihthres,stdk
            if(jsum .ge. ihthres) then
c             passed threshold test - put interim into final
              nhits(itarg,imo) = nhits(itarg,imo) + jsum
              ifound = ifound + 1
              if(icdebug .gt. 0) then
                ahigh = akhigh(imo)/khits(imo)
                astd = ahigh
                if(method .ne. 3) then
                  procstr = 'UCHGPT'
                  write(6,3000) itarg, ntstn(itarg), procstr, nloop, iy,
     *              im, jsum, ahigh, astd, ibracket, npair, ihthres
                else
                  procstr = 'FINSHF'
                  write(6,3000) itarg, ntstn(itarg), procstr, nloop, iy,
     *              im, jsum, ahigh, astd, ibracket, npair, ihthres
                endif    
 3000           format(i5,1x,a6,'-',a6,' KW',i1,' at ',i4,i3,i4,2f6.2,
     *            3i3)
              endif
            endif  
          endif ! end of using interim hit
        enddo ! end if examining interim hit array  
  180 enddo
  
      if(method .ne. 3) then
        print *,' Undoc filter: ',ifound
      else
        print *,' BIC final filter: ', ifound
      endif    
      if(icdebug .lt. 2) goto 210
      print *,' Exact wrt last: ', iexact
      do k1 = 1, numsubs, 100
        k2 = k1 + 99
        if(k2 .gt. maxstns) k2 = maxstns
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 190
          itarg = nindx(k,1)
          outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(1:3)// '|'
  190   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 195
          itarg = nindx(k,1)
          outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(4:6)// '|'
  195   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        nofound = 0
        do imo = 1, nmo
          call imo2iym(iy,im,imo)
          write(outstr, '(i4,i3,i5," |")') iy,im,imo
          mofound = 0
          do k = k1, k2
            if(ntr(k) .eq. 0) goto 200
            itarg = nindx(k,1)
            if(nhits(itarg,imo) .eq. 0) then
              outstr = outstr(1:lnblnk(outstr)) // '---'
            else
              write(tmpstr,'(i3)') nhits(itarg,imo)
              outstr = outstr(1:lnblnk(outstr)) // tmpstr
              mofound = mofound + 1
            endif  
            outstr = outstr(1:lnblnk(outstr)) // '|'
          enddo  
  200     if(mofound .gt. 0) then
            if(nofound .gt. 0) then
c              write(6,*) nofound
              nofound = 0
            endif  
            write(6,*) outstr(1:lnblnk(outstr))
          else
            nofound = nofound + 1  
          endif  
        enddo
      enddo
      
c     4th filter (alignmoves) has migrated to the estamt routine

  210 return
      end

c     =======================================================================
      subroutine maxnet(numstn, rrdata, mxstn)

      INCLUDE 'inhomog.parm.mthly.incl'  

c         maxnet determines the combination of neighbors that
c         maximizes the number of station-years with at least
c         minyr years and minstn stations in the subnet
c         returns a mxstn station array giving the indices in
c         the nindx(itarg, ...) of the maximized subnet.

c     Subroutine arguments -
c     rrdata - input neighborhood data array
c     mxstn - output array of rrtemp stations used for best network
c               (1 = use; 2 = reject)

c     Local variables - 
c     mxcomb - combination array with max sta-yr where all sta have each yr
c     mxdat - max station X year
c     scomb - current combination being summarized
c     ndat - current combo number of data vals
c     nyr - current combo number of years
c     ncomb - current combo number of stations
c     istop - indicates when all combinations have been evaluated (=1)
c     iny - number of data vals in current year

      real rrdata(nstns,nmo)
      integer mxstn(nstns)
      integer mxcomb(nstns), scomb(nstns)
      
c     current values for minimum number of years and stations
c     LV's DW thresholds define the minyr = 15
c     No need to do less than two stations (one stn already done in pairs)
      minyr = 15
      minstn = 2
      
      istop = 0
      
      do ind = 1, nstns 
        mxstn(ind) = 0     
        mxcomb(ind) = 0
        scomb(ind) = 0
      enddo  
      
      mxdat = 0
      mxyr = 0
      lcomb = minstn - 1
      mxnum = lcomb
      
      do while (istop .eq. 0)
        call combine(numstn, scomb, istop)
        if(istop .eq. 1) goto 100
        
c       find the last station, if not the candidate, then skip
         do ip = 1, numstn
           if(scomb(ip) .eq. 0) then
             ncp = ip-1
             goto 10
           endif  
         enddo
         ncp = ip -1
   10    if(scomb(ncp) .ne. 1) goto 50

c       no need to look at a combo until the minstn > candidate
        if(scomb(minstn) .gt. 1) then
          ndat = 0
          nyr = 0
          ncomb = 0
          do imo = 1, nmo
            iny = 0
c           no need to count if the candidate data is missing
            if(rrdata(1,imo).gt.amiss+1.) then
              do ind = 1, numstn
c               candidate is the last stn in the current combo
                if(scomb(ind) .eq. 1) goto 20
                if(rrdata(scomb(ind), imo).gt.amiss+1.) then
                  iny = iny + 1
                else
c                 if any station in current combo is missing, 
c                 do not use any data for that year
                  iny = 0
                  goto 40
                endif
              enddo ! end of station loop

c             accumulate number of data vals and number of years
   20         ndat = ndat + iny
              nyr = nyr + 1
              ncomb = ind - 1
            endif
          
   40     enddo ! end of year loop

c         to cut down on the number of combinations tested, if
c           there was no network combo found with ncomb-1 stations
c           with nyr ge maxyr then there won't be any in higher
c           combo's either.
          if(ncomb .gt. lcomb) then
            if(mxnum .ne. lcomb) goto 100
            lcomb = ncomb
          endif  
   
c         save this combination if the number of sta-yr is highest so far
c         putting in a priority wrt length of series
          if(nyr.ge.minyr .and. ncomb.ge.minstn .and. 
     *      ndat.ge.mxdat) then
c           for a given combination order (ncomb = mxnum) 
c             later combinations have less correlation
c           each replaced network must have at least as 
c             many years as the last best
            if((ncomb .gt. mxnum .and. nyr .ge. mxyr) .or.
     *        (ncomb .eq. mxnum .and. nyr .gt. mxyr)) then
              mxdat = ndat
              mxyr = nyr
              mxnum = ncomb
              do ind = 1, nstns
                mxcomb(ind) = scomb(ind)
              enddo
c            else
c              print *,' Dismissed: ', nyr, ncomb, ndat  
            endif  
          endif    
        endif 
   50 enddo ! end combination loop
      
c     generate output array
  100 do ind = 1, nstns
        if(mxcomb(ind) .gt. 0) mxstn(mxcomb(ind)) = 1
      enddo
      
      print *,'Max Merged net: ',mxnum,' stations by ', mxyr,' years'   

      return
      end

C***********************************************************************
c      COMBINE - generate successive combinations from 1 to NP
c        integers into the S array. When ISTOP == 1 all of the
c        combinations have been generated
C***********************************************************************
       subroutine combine(np, s, istop)
       INCLUDE 'inhomog.parm.mthly.incl'
       integer s(nstns)
       istop = 0
         
       s(1) = s(1) + 1
       if(s(1) .gt. np) then
         do ind = 2, np
           if(s(ind) .lt. np-ind+1) then
             s(ind) = s(ind) + 1
             do ii = ind - 1, 1, -1
               s(ii) = s(ii+1) + 1
             enddo
             goto 10
           endif
         enddo
c        if we get here, the last permutation has been made
         istop = 1
         return
       endif

   10  continue
c       print *, (s(ind),ind= np,1,-1)
c         pause
       return
       end
                 
C***********************************************************************
C End of Subroutine COMBINE.                                         *
C***********************************************************************

c     convert month indx to year,month
      subroutine imo2iym(iy,im,imo)
      INCLUDE 'inhomog.parm.mthly.incl'
      INCLUDE 'inhomog.comm.mthly.incl'
      if(itimeres .eq. 0) then
c       annual
c        iy = imo + begyr - 1
        iy = imo + begyr
        im = 13
      else
c       monthly
        iy = int((imo-1)/12) + (begyr + 1)
        im = mod(imo, 12)
        if(im .eq. 0) im = 12
      endif  
      return
      end
      
c     convert year,month to month indx
      subroutine iym2imo(iy,im,imo)
      INCLUDE 'inhomog.parm.mthly.incl'
      INCLUDE 'inhomog.comm.mthly.incl'
      if(itimeres .eq. 0) then
c       annual
c        imo = iy - begyr + 1
        imo = iy - begyr
      else
c       monthly      
        imo = (iy-begyr-1)*12 + im
      endif  
      return
      end
      
c =======================================================================

      subroutine estamt(ntstn,idunit)
      
c     A major rework was accomplished on 07 Sept. 06
c       The major steps in determining the best adjustment value
c       for each station/changepoint. Entire network undergoes each of
c       the following processes. In order:
c     1) Remove unusable data. Align moves with respect to non-missing
c        data and compress out changes that are too close AND the data
c        between them.
c     2) ISTEP=2 processing begins the adjustment process by removing
c         the Non-significant changepoints to lengthen segments.
c     3) NPASS (:= ISTEP=3) finishes the adjustment process by testing
c         for the minimum number of months in a segment and number
c         of neighbors with which the difference series can be examined.
c     4) Final adjusted output is written
      
c     WARNING: the changepoints in the nhits array are assumed to have
c       been processed through the alignmoves routine while going through
c       confirmfilter. That is: all of the chgpt position alignment and 
c       data adjustment for too short segments has been done.
      
      INCLUDE 'inhomog.parm.mthly.incl'
      INCLUDE 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'
      
c      ----- NPASS includes the "Not sig" step as of v21 (18dec06) ----------
      parameter (npass = 2)
c      parameter (npass = 1)
c      ----------- post-One Pass Test ----------
c      parameter (npass = 4)
      parameter (mlsx3 = minlen * 3)
      parameter (mlsx2 = minlen * 2)
      
c     temporary output file name
      character*132 fname
      
      character*6 ntstn(maxstns), procstr

c     adjustment arrays for output - one station at a time
      real outtemp(begyr:endyr, 13),adjtemp(begyr:endyr, 13), 
     *   contemp(begyr:endyr, 13)
      character*1 monflag(begyr:endyr, 13), adjflag(begyr:endyr, 13),
     *   conflag(begyr:endyr, 13), aflg
     
      character*2 otag/'WM'/
      character*500 outstr,tmpstr
      character*1 use
      character*66 chgptstr

c     arrays for alignmoves
      integer move(ninh), mday(ninh)
      real amt(ninh)
      
c     for monthly anomaly calc
      real avgmth(maxstns, 13), smth(13)
      integer nmth(13)
      
c     holding array for original data
      real filtemp(maxstns, begyr:endyr, 13)

c     begin & end month indices for segment 1 & 2
      integer beg1, end1, beg2, end2
      
c     for the final stat test (including BIC) routine
      real diff(nmo), qx(nmo), qy(nmo)
c     1st diff series 
      real qx1d(nmo), qy1d(nmo)
      integer idiff(nmo), iqtype(nstns), mknt2(2)
      real rmuq(2), rslpq(2, nstns), ampq(2), phiq(2)
      real trend(nstns, ninh), tradj(nstns, ninh),
     *  trcor(nstns, ninh)
      real tqmin(nstns, ninh)
      integer itrqt(nstns, ninh), spanob(nstns, ninh)
      character*1 isegused(nstns, ninh)

c     set the "D" flag for all station/neigh/months
      character*1 cedit(maxstns,nstns,nmo)
      character*2 outid/'  '/

c     work arrays for spans and chgpt removal
      integer ispanq(maxstns,ninh), iremchg(maxstns,ninh)
c     trimmed adj arrays for model >= 3
      real tchgs(nstns), fchgs(nstns), tccor(nstns)
      integer iqchgs(nstns), ichgs(nstns), iseg1(nstns), iseg2(nstns),
     *  itind(nstns)

c     list of unstable dates per station and number of unstable chgpts
      integer iunstbl(maxstns,ninh), nunstbl(maxstns), nundat(ninh)

      character*1 iused(nstns)
      real tseg1(nstns), tseg2(nstns)
      
c     for inner-quartile range
      real x(0:nmo),tadj(nstns), tcorr(nstns)
      real qtarg(nmo), qpair(nmo)
      
c     beg/end period of record
      integer ifirstyr(maxstns), ilastyr(maxstns)
      
c     network changepoint arrays
c      ichgmo - month index of the changepoints for the entire network
c      nchgin - number of changepoints in each target station
c      PRECEEDING ARRAYS DO NOT CHANGE IN ADJUSTMENT LOOPS
c      ichgout - status of each changepoint per pass/step
c                0 - to be evaluated
c                +1 to +7 - MINBIC type model solution
c       num - number of targ-neigh pairs used to form adjustment estimate
c       adj - adjustment value (trimmed median)
c       std - error of adjustment (1.5 * inner-quartile)
      integer ichgmo(maxstns, ninh), nchgin(maxstns)
      integer ichgout(maxstns, ninh), num(maxstns, ninh)
      real adj(maxstns, ninh), std(maxstns, ninh)

c     list of "station minbic hits" limits to use in each ipass
c      ----- NPASS includes the "Not sig" step as of v21 (18dec06) ----------
      integer iadjlist(npass)/3, minsta/
      integer iadjlen(npass)/minlen, minlen/

      character*13 tstr

      numrem = 0
      NLOOP = 99
      
c     correlation threshold for lpairs
      corthres = 0.0

c     debug output level
c      iedebug = 1
c      icdebug = 2
      iedebug = 0
      icdebug = 0

c     initialize temporary EDIT array for all station/neigh/months
      do itarg = 1, maxstns
        do ipair = 1, nstns
          do imo = 1, nmo
            cedit(itarg,ipair,imo) = ' '
          enddo
        enddo
      enddo      

c     initialize metadata removal array
      do itarg = 1, maxstns
        do ichg = 1, ninh
          iremchg(itarg, ichg) = 0
        enddo
      enddo    

c     filter temp array to remove "ndelete < ndellim" data
      do k = 1, numsubs
        itarg = nindx(k,1)
        do iy = begyr+1, endyr
c         if annual time resolution (lomth,himth=13)
          do im = lomth, himth
            if(temp(itarg,iy,im).gt.amiss+1.) then
              call iym2imo(iy,im,imo)
c             flag suspect data where sum(ndelete) >= ndellim for adj est
              idel = 0
              do ipair = 1, nstns
                if(ndelete(itarg,ipair,imo).ne.czero) then
                  idel = idel + 1
                  cedit(itarg,ipair,imo) = 'd'
                endif  
              enddo  
              if(idel .ge. ndellim) then
c               set to missing               
                tflg(itarg,iy,im) = 'd'
              endif  
            endif
          enddo
        enddo
      enddo  
      
c     initialize unstable date number and count
      do k = 1, numsubs
        nunstbl(k) = 1
        iunstbl(k,1) = nmo
        do inh = 2, ninh
          iunstbl(k,inh) = 0
        enddo
      enddo    

C     ---------------------------------------------------------------------
c     the fourth filter (moved from confirmfilt)

c     Since the amplitude est MUST rely upon a minimum of MINLEN (=18) months
c     to get even close to a reliable estimate at this point, it is assumed
c     that the UCP are as good as the SHF - 
c     Therefore, align moves with respect to non-missing data and compress out 
c      changes that are too close AND the data between them (i.e. less than
c      MINLEN in the "temp" array)

c     for each station .....
      do k = 1, numsubs
        if(ntr(k) .eq. 0) goto 5
        itarg = nindx(k,1)
        
        do inh = 1, ninh
          move(inh) = 0
          amt(inh) = 0.0
          mday(inh) = 0
        enddo     

        movnum = 0        
c       ... gen arrays for alignment
        do imo = 1, nmo
          if(nhits(itarg,imo) .ne. 0) then
c           ... increment segment for each move
            movnum = movnum + 1
c           ... define alignment arrays
            move(movnum) = imo
            amt(movnum) = nhits(itarg,imo)
            mday(movnum) = 31
          endif
        enddo
        
        if(movnum .ne. 0)
     *    call alignmoves(ntstn(itarg),itarg,0,movnum,move,amt,mday)
          
c       put the data back into nhits
        do imo = 1, nmo
          nhits(itarg,imo) = 0
        enddo
        
        do imove = 1, movnum
          if(iedebug .ge. 2) print *,imove, move(imove), amt(imove)
          i = move(imove)
          j = amt(imove)
          if(iedebug .ge. 2) print *,i,j
          nhits(itarg,i) = j
        enddo    

    5 enddo
  
c     move incoming temp array into filtemp for output
      do itarg = 1, numsubs
        do iy = begyr+1, endyr
          do im = 1, 13
            filtemp(itarg, iy, im) = temp(itarg, iy, im)
          enddo
        enddo  
      enddo
      
c     print the best BIC model return from confirmfilt 
      if(iedebug .gt. 3) then
        do itarg = 1, numsubs
          do imo = 1, nmo
            if(nhits(itarg,imo).gt.0) then
c             generate adjustment amound output (iconfirm > 2 assumed)
              inum = nchgpt(itarg,imo)
              avg = schgpt(itarg,imo) / inum
              avgz = zchgpt(itarg,imo) / inum
              call imo2iym(iy,im,imo)
             write(6,'(a6,i4," BICFIN:",2i5," off",2f7.2,i4," in",i4)')
     *         ntstn(itarg),itarg,iy,im,avg,avgz,inum,nhits(itarg,imo)
            endif   
          enddo
        enddo        
      endif

      print *,' Alignmoves & Ndelete filter: ', ifound
      if(icdebug .lt. 3) goto 25
      print *,' Exact wrt last: ', iexact
      do k1 = 1, numsubs, 100
        k2 = k1 + 99
        if(k2 .gt. maxstns) k2 = maxstns
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 10
          itarg = nindx(k,1)
          outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(1:3)// '|'
   10   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        write(outstr, '("             |")')
        do k = k1, k2
          if(ntr(k) .eq. 0) goto 15
          itarg = nindx(k,1)
          outstr=outstr(1:lnblnk(outstr))//ntstn(itarg)(4:6)// '|'
   15   enddo  
        write(6,*) outstr(1:lnblnk(outstr))
        do imo = 1, nmo
          call imo2iym(iy,im,imo)
          write(outstr, '(i4,i3,i5," |")') iy,im,imo
          do k = k1, k2
            if(ntr(k) .eq. 0) goto 20
            itarg = nindx(k,1)
            if(nhits(itarg,imo) .eq. 0) then
              if(temp(itarg, iy, im) .lt. amiss+1.0) then
                outstr = outstr(1:lnblnk(outstr)) // '-X-'
              else  
                outstr = outstr(1:lnblnk(outstr)) // '---'
              endif  
            else
              write(tmpstr,'(i3)') nhits(itarg,imo)
              outstr = outstr(1:lnblnk(outstr)) // tmpstr
            endif  
            outstr = outstr(1:lnblnk(outstr)) // '|'
          enddo  
   20     write(6,*) outstr(1:lnblnk(outstr))
        enddo
      enddo  

c     generate the series monthly averages for all of the stations
   25 do k = 1, numsubs
        itarg = nindx(k,1)
c       initialize the beg/end period of record
        ifirstyr(itarg) = 0
        ilastyr(itarg) = 0
c       accumulate sums and nums for monthly averages
        do im = 1, 13
          smth(im) = 0.0
          nmth(im) = 0
        enddo  
        do iy = begyr+1, endyr
c         if annual time resolution (lomth,himth=13)
          do im = lomth, himth
            if(temp(itarg,iy,im).gt.amiss+1.) then
              call iym2imo(iy,im,imo)
              if(ifirstyr(itarg) .eq. 0) ifirstyr(itarg) = imo
              ilastyr(itarg) = imo
              smth(im) = smth(im) + temp(itarg,iy,im)
              nmth(im) = nmth(im) + 1
            endif
          enddo
        enddo
c       calculate monthly averages
        do im = lomth, himth
          avgmth(itarg,im) = smth(im) / nmth(im)
        enddo
      enddo    

      do k = 1, maxstns
        nchgin(k) = 0

        do inh = 1, ninh
          ichgmo(k, inh) = 0
          adj(k, inh) = 0.0
          std(k, inh) = 0.0
          num(k, inh) = 0
          ispanq(k, inh) = 0
        enddo  
      enddo  

c     initialize the changepoint arrays AFTER passing through alignmoves
c       update in the adjloop of the step == 2 pass below
c       go through all of input chgpts (first/last are beg/end por)
c         accumulating the output chgpts & adjustments
      do k = 1, numsubs
c       itarg is the index of the sub-network (candidate)
        itarg = nindx(k,1)

c       Find the confirmfilt output changepoints in nhits 
        nchg = 1
        ichgmo(itarg, nchg) = ifirstyr(itarg)
        call imo2iym(iy,im,ichgmo(itarg, nchg))
        if(iedebug .ge. 3)
     *    write(6,'(i6," Chgpt  Begin year: ",4i6)') itarg, nchg, iy, 
     *      im, ichgmo(itarg, nchg)
        do imo = 1, nmo
          if(nhits(itarg,imo).gt.0) then
            nchg = nchg + 1
            ichgmo(itarg, nchg) = imo
            call imo2iym(iy,im,ichgmo(itarg, nchg))
            if(iedebug .ge. 3)
     *        write(6,'(i6," Confirmfilt chgpt: ",5i6)') itarg, nchg,
     *          iy,im,ichgmo(itarg, nchg),nhits(itarg, imo)
          endif
        enddo
c       inserted if around following lines to fix end series UNSTBL problem
        if(ichgmo(itarg, nchg) .ne. ilastyr(itarg)) then
          nchg = nchg + 1
          ichgmo(itarg, nchg) = ilastyr(itarg)
          call imo2iym(iy,im,ichgmo(itarg, nchg))
          if(iedebug .ge. 3)
     *      write(6,'(i6," Chgpt    End year: ",4i6)') itarg, nchg, iy,
     *        im, ichgmo(itarg, nchg)
        endif
        nchgin(itarg) = nchg
      enddo

c     The subnetwork processing became a multi-step process plus a
c      "post-process pass" to manage:
c       1) problems with the documented changepoints with NO undocumented
c             support. 
c       2) determine the best amplitude estimation for each confirmed
c             changepoint
c                                                         v15b - 16may06

c      ------------- Sig/Adj test and final solution (21dec06) -----------
c                       Note: istep == 1 is not used
      do istep = 2, 3

c      ichgout initialization relocated here v20d (18dec06)
c        used temporarily in the "steps" but kept in the "passes"
        do k = 1, maxstns
          do inh = 1, ninh
            ichgout(k,inh) = 0
          enddo  
        enddo
      
c       ------------ ipass loopback for istep == 2 -------------
   27   if(istep .eq. 1) then
          goto 87
        else if(istep .eq. 2) then
          print *,' -------------- NOT SIG REMOVAL ---------'
          tstr = 'Not sig: '
          outid = 'NS'
          ipass = 1
          iminlen = minlen
          numclim = 3
        else
          print *,' -------------- ADJUST DISCONTINUITY STEP ---------'
          tstr = 'Good Neigh:  '
          ipass = ipass + 1
          outid = 'GN'
          if(ipass .gt. npass) goto 90
          iminlen = iadjlen(ipass)
          numclim = iadjlist(ipass)
          print *,'Adjpass, iminlen, numclim: ', ipass, iminlen, 
     *      numclim
          if(ipass .eq. npass) then
            print *,' -------------------- NPASS -----------------'
            tstr = 'Dstep Dtrend: '
            outid = 'WM'
            iedebug = 1
            icdebug = 2
          endif  
        endif  

c       process each subnetwork (candidate/target) one at a time
        if(istep .eq. 2) write(6,'(a,a,a)') 
     *  '  NET  STN    FILT TECH        ------ AFTER -----   ',
     *  '    ----- BEFORE -----              AVG   STD NEIGH   CUM ',
     *  'NHITS'
     
        do k = 1, numsubs
c         itarg is the index of the sub-network (candidate)
          itarg = nindx(k,1)
          nchg = nchgin(itarg)
          iloop = 0

          do istn = 1, nstns
            do ichg = 1, ninh
              tradj(istn, ichg) = amiss
              trcor(istn, ichg) = amiss
              trend(istn, ichg) = amiss
              tqmin(istn, ichg) = 9999.
              isegused(istn,ichg) = ' '
              spanob(istn, ichg) = 0
              itrqt(istn, ichg) = 0
            enddo
          enddo
          
          minindx = nindx(k,minsta+1)
          if(minindx .eq. 0) then
            print *,' Neighbors < minsta: ', ntstn(itarg)
            go to 77
          endif  

c         List the confirmfilt output changepoints
          if(iedebug .ge. 1 .and. istep .ge. 2) then
            do ichg = 1, nchg
              call imo2iym(iy,im,ichgmo(itarg, ichg))
              if(ichgmo(itarg,ichg) .le. 0) then
                print *,' ichgmo error: ', itarg, ntstn(itarg), ichg, 
     *            iy, im, ichgmo(itarg, ichg)
                nchg = 1
                goto 30
              else
                write(6,'(i3,i5,1x,a," Estamt chgin: ",2i5,i2.2,3i6,
     *            2f7.2)') ipass, itarg, ntstn(itarg), ichg, 
     *            iy, im, ichgmo(itarg,ichg), 
     *            nhits(itarg,ichgmo(itarg,ichg)), ichgout(itarg,ichg),
     *            adj(itarg,ichg),std(itarg,ichg)
              endif
            enddo  
          endif
            
c         ---------------- station chgpt loop --------------------
c         find all of the "not done" changepoints (ichg) up to this point
c           and the "not done" changepoints (ichg1 & ichg2) that bracket it
c         Note: ichgmo(1) & ichgmo(nchg) are the begin/end por
   30     lchg = nchg
          ichg = nchg - 1
          siglast = 0.0
          do while (ichg .ge. 2)
c            print *, 'Incoming: ', ichg, lchg

c           reset current chgpt by going backward from the current ichg 
c             to the first chgpt not done
            do jchg = lchg-1, 2, -1
              if(ichgout(itarg, jchg) .eq. 0) goto 32
            enddo
   32       ichg = jchg
            if(ichg .eq. 1) goto 75
            
c           go forward to the latter chgpt or first chgpt not done (==0)
            do jchg = ichg+1, lchg
              if(ichgout(itarg, jchg) .eq. 0) goto 33
            enddo
   33       if(jchg .gt. lchg) jchg = lchg
            ichg2 = jchg
   
c           finally, keep going back to the next earlier chgpt not done (==0)
            do jchg = ichg-1, 2, -1
              if(ichgout(itarg, jchg) .eq. 0) goto 34
            enddo
   34       ichg1 = jchg
c           set last changepoint to current changepoint for next iteration
            lchg = ichg

c           expanded IF to fix beg series UNSTBL problem
            if(ichg1 .eq. 1) then
              beg1 = ichgmo(itarg, 1)
            else
              beg1 = ichgmo(itarg, ichg1) + 1
            endif
            end1 = ichgmo(itarg, ichg)
            beg2 = end1 + 1
            end2 = ichgmo(itarg, ichg2)
            print *, 'Oriented: ', ichg1, ichg, ichg2,
     *        beg1, end1, beg2, end2
            
c           when checking the significance of a metadata changepoint (istep=1),
c             only evaluate those with no supporting undoc stats (nhits=100)
            if(istep .eq. 1 .and. nhits(itarg,end1) .ne. 100) goto 70

c           get all of the yr/mth for print
            call imo2iym(iyb1,imb1,beg1)
            call imo2iym(iye1,ime1,end1)
            call imo2iym(iyb2,imb2,beg2)
            call imo2iym(iye2,ime2,end2)
            write(chgptstr,'("  Win1: ", 2(i5, i2.2, i5)," to ",
     *         "Win2: ", 2(i5, i2.2, i5))')
     *         iyb1, imb1, beg1, iye1, ime1, end1,
     *         iyb2, imb2, beg2, iye2, ime2, end2
     
c          ---------- This section accumulates target-neighbor ------------
c          ----------------------- comparisons --------------------------

c           see if there are enough homogenous data in the target 
c           first forward - process the after chgpt (second) window
            ibeg2 = 0
            iend2 = 0
            itot2 = 0
            do imo = beg2, end2
              call imo2iym(iy,im,imo)
              if(temp(itarg,iy,im).gt.amiss+1.) then
c               set the first month used
                if(ibeg2 .eq. 0) ibeg2 = imo
c               set last used
                iend2 = imo
                itot2 = itot2 + 1
              endif  
            enddo  

c           if the segment length of the targ station is too short for this
c             adj pass, skip this chgpt - for now....
c           however, if at NPASS then remove segment and attempt to span
            if(itot2 .lt. iminlen) then
              if(ipass .lt. npass) then
                write(6,'("Adjpass seg2 short ",2i4,a," ",a,i5,
     *            " -skip")') ipass,itarg,ntstn(itarg),chgptstr,itot2
                goto 70
              else  
                write(6,'("Adjpass seg2 short ",2i4,a," ",a,i5,
     *            " -remove")') ipass,itarg,ntstn(itarg),chgptstr,itot2
                iremseg = 2
                goto 67
              endif  
            endif  

c           now go back in time - process the before chgpt (first window)
            ibeg1 = 0
            iend1 = 0
            itot1 = 0
            do imo = end1, beg1, -1
              call imo2iym(iy,im,imo)
              if(temp(itarg,iy,im).gt.amiss+1.) then
c               set the first month used
                if(ibeg1 .eq. 0) ibeg1 = imo
c               set last used
                iend1 = imo
                itot1 = itot1 + 1
              endif  
            enddo  

c           if the segment length of the targ station is too short for this
c             adj pass, skip this chgpt - for now....
c           however, if at NPASS then remove segment and attempt to span
            if(itot1 .lt. iminlen) then
              if(ipass .lt. npass) then
                write(6,'("Adjpass seg1 short ",2i4,a," ",a,i5,
     *            " -skip")') ipass,itarg,ntstn(itarg),chgptstr,itot2
                goto 70
              else  
                write(6,'("Adjpass seg1 short ",2i4,a," ",a,i5,
     *            " -remove")') ipass,itarg,ntstn(itarg),chgptstr,itot1
                iremseg = 1     
                goto 67
              endif  
            endif  

c           number of neighbor pairs tested with adjustments
            numc = 0

c           go through all of the neighbors
            do it2pair = 2, nstns 
c             initialize paired diff array
              do imo = 1, nmo
                diff(imo) = amiss
                qtarg(imo) = amiss
                qpair(imo) = amiss
                qx(imo) = amiss
                qy(imo) = amiss
                qx1d(imo) = amiss
                qy1d(imo) = amiss
              enddo  
c             ipair - the index of the "it2pair" station in the "itarg" subnet
              ipair = nindx(itarg, it2pair)
c             if there are no more neighbors, skip out
              if(ipair .eq. 0) goto 55
          
c             and the neighbor to make the before/after windows around the
c             chgpt location
                  
c             first forward - process the after chgpt (second) window
              num2 = 0
              imo2b = 0
              imo2e = 0
c             accumulate paired data for the most recent window(2)
              ncedit = 0
              do imo = beg2, end2
                call imo2iym(iy,im,imo)
c               sum and count the differences
c               the target-pair completeness is defined for every non-missing
c               target month there is a pair month
                if(cedit(itarg,it2pair,imo) .ne. ' ') ncedit = ncedit+1
                if(temp(itarg,iy,im).gt.amiss+1. .and.
     *            temp(ipair,iy,im).gt.amiss+1. .and. 
     *            cedit(itarg,it2pair,imo) .eq. ' ') then
c                 set the first month used
                  if(imo2b .eq. 0) imo2b = imo
c                 set last used
                  imo2e = imo
c                 gen & use monthly anomalies for paired adjustment est
                  qtarg(imo) = temp(itarg,iy,im) - avgmth(itarg,im)
                  qpair(imo) = temp(ipair,iy,im) - avgmth(ipair, im)
                  diff(imo) = qtarg(imo) - qpair(imo)
                  num2 = num2 + 1
                endif
c               loop out if neighbor hits another chgpt
c               this test is done after this month is used because when
c               going forward the chgpt location is the last in this segment!
c               version v20c ------ 15dec06
                if(nhits(ipair,imo).gt.0) then 
                  if(iedebug .gt. 1)
     *              write(6,'("CHG2: ",a,i3," num,edit,2b,2e,imo,nhits",
     *                6i5)') ntstn(ipair),it2pair,num2,ncedit,
     *                imo2b,imo2e,imo,nhits(ipair,imo)
c                  goto 50
                  goto 35
                endif  
              enddo

   35         continue
              if(num2 .lt. iminlen) then
                if(iedebug .gt. 3)
     *            write(6,'("Low2: ",a,i3," num,edit,2b,2e,imo,nhits",
     *              6i5)') ntstn(ipair),it2pair,num2,ncedit,
     *              imo2b,imo2e,imo,nhits(ipair,imo)
                goto 50
              endif
              
c             now go back in time - process the before chgpt (first window)
              num1 = 0
              imo1b = 0
              imo1e = 0
              ncedit = 0
              do imo = end1, beg1, -1
c               loop out if either of the stations hit another chgpt
                if(imo .ne. end1 .and. nhits(itarg,imo).gt.0) goto 40
c               Or neighbor has changepoint in target segment
c               if(nhits(ipair,imo).gt.0) goto 40
c               version v20c ------ 15dec06
                if(nhits(ipair,imo).gt.0) then 
                  if(iedebug .gt. 1)
     *              write(6,'("CHG1: ",a,i3," num,edit,1b,1e,imo,nhits",
     *                6i5)') ntstn(ipair),it2pair,num1,ncedit,
     *                imo1b,imo1e,imo,nhits(ipair,imo)
c                  goto 50
                  goto 40
                endif  
                call imo2iym(iy,im,imo)
c               sum and count the differences
                if(cedit(itarg,it2pair,imo) .ne. ' ') ncedit = ncedit+1
                if(temp(itarg,iy,im).gt.amiss+1. .and. 
     *            temp(ipair,iy,im).gt.amiss+1. .and. 
     *            cedit(itarg,it2pair,imo) .eq. ' ') then
c                 set the first month used
                  if(imo1e .eq. 0) imo1e = imo
c                 set last used
                  imo1b = imo
c                 gen & use monthly anomalies for paired adjustment est
                  qtarg(imo) = temp(itarg,iy,im) - avgmth(itarg,im)
                  qpair(imo) = temp(ipair,iy,im) - avgmth(ipair, im)
                  diff(imo) = qtarg(imo) - qpair(imo)
                  num1 = num1 + 1
                endif    
              enddo

   40         continue
              if(num1 .lt. iminlen) then
                if(iedebug .gt. 3)
     *            write(6,'("Low1: ",a,i3," num,edit,1b,1e,imo,nhits",
     *              6i5)') ntstn(ipair),it2pair,num1,ncedit,
     *              imo1b,imo1e,imo,nhits(ipair,imo)
                goto 50
              endif
   
              numc = numc + 1

c             scan the difference array, pass segment with missing data
c               into minbic
              jmo = 0
              do imo = 1, nmo
                if(diff(imo) .ne. amiss) then
                  if(jmo .eq. 0) jmobeg = imo
                  jmo = imo - jmobeg + 1
                  qy(jmo) = diff(imo)
                  qx(jmo) = imo
                  if(imo .le. end2) jend2 = jmo
                  if(imo .le. end1) jend1 = jmo
                endif
              enddo
            
c             ONLY TPR0 from minbic but return slopes anyway (for study).
              call minbic(6, qx, qy, jend1, jend2, tcrit, tchgpt, 
     *          qmin, tadj(numc), rmuq, rslpq(1,numc), rsseq, 5, 
     *          iqtype(numc), mknt2, ifail)
     
c             let's see what first difference correlation's come out for pairs
              call frstdif(nmo, qtarg, qpair, amiss, jmo, qx1d, qy1d)
              call correl(jmo, qx1d, qy1d, tcorr(numc), amiss)

              if(iedebug .gt. 0) then
                write(6,'(2a, "-", 2a, 4f7.2, 2f7.3, 4i5)')
     *            tstr, ntstn(itarg), ntstn(ipair), chgptstr, tcrit,
     *            tchgpt, tadj(numc), tcorr(numc), rslpq(1,numc),  
     *            rslpq(2,numc), iqtype(numc), mknt2, numc
              endif
     
c             keep adjustment for each neighbor/segment
c             set/reset trend for each neighbor/segment
c               the first segment is the before (ichg) segment
c               the second segment is the after (ichg+1) segment
c                 examining previous debug output indicates that since:
c                 the difference series is same, the output is same
              tradj(it2pair,ichg) = tadj(numc)
              trcor(it2pair,ichg) = tcorr(numc)
              tqmin(it2pair,ichg) = qmin
              itrqt(it2pair,ichg) = iqtype(numc)
              trend(it2pair,ichg) = rslpq(1,numc)
              spanob(it2pair,ichg) = mknt2(1)
c             Segment 2 has already been run as Segment 1 for the last changepoint
c               if slopes are different use the one calculated from the largest # obs
              if(trend(it2pair,ichg2) .eq. amiss) then
                trend(it2pair,ichg2) = rslpq(2,numc)
                spanob(it2pair,ichg2) = mknt2(2)
              else if(trend(it2pair,ichg2) .ne. rslpq(2,numc)) then
                if(iedebug .ge. 2)
     *            write(6,'(" Seg2 diff: ",2i4," old:",f7.2,i4," new:",
     *            f7.2,i4)') it2pair, ichg2, trend(it2pair, ichg2), 
     *            spanob(it2pair,ichg2), rslpq(2,numc), mknt2(2)
c               if current count greater than this time- replace older with new slope
                if(mknt2(2) .gt. spanob(it2pair,ichg2)) then
                  trend(it2pair,ichg2) = rslpq(2,numc)
                  spanob(it2pair,ichg2) = mknt2(2)
                endif  
              endif
              if(icdebug .gt. 2) 
     *          print *,'itarg,ipair,ichg,numc,iqt,adj,trends:',
     *          itarg,it2pair,ichg,numc, 
     *          itrqt(it2pair,ichg),tradj(it2pair,ichg),
     *          trend(it2pair,ichg),trend(it2pair,ichg2)

   50       enddo ! end of the neighbors loop

   55       continue
   
c           ---------------- End of accumulating paired chgpt --------------
c           ------------------- comparisons for network --------------------
c           ----------------------------------------------------------------

c           ----------------------------------------------------------------
c           ------------------- This section determines  -------------------
c           -------------------------- adjustment --------------------------

            npairs = 0

c           go through all of the neighbors
            do it2pair = 2, nstns 
c             ipair - the index of the "it2pair" station in the "itarg" subnet
              ipair = nindx(itarg, it2pair)
c             if there are no more neighbors, skip out
              if(ipair .eq. 0) goto 57

c             setup trends and offsets for median calculations
              if(trend(it2pair,ichg) .eq. amiss .or. 
     *          trend(it2pair,ichg2) .eq. amiss) then
                iused(it2pair) = 'M'  
              else
                npairs = npairs + 1
                tchgs(npairs) = tradj(it2pair,ichg)
                tccor(npairs) = trcor(it2pair,ichg)
                tseg1(npairs) = trend(it2pair,ichg)
                tseg2(npairs) = trend(it2pair,ichg2)
                ichgs(npairs) = it2pair
                iseg1(npairs) = it2pair
                iseg2(npairs) = it2pair
                itind(npairs) = npairs
                iused(it2pair) = 'U'
              endif
            enddo

c           in step == 3, iadjlist(ipass) determines the number of 
c             "Good Neigh" estimates needed for this pass to process the
c             amplitude estimation for this chgpt.
   57       if(npairs .lt. numclim) then
              if(ipass .lt. npass) then
                write(6,'("Adjpass numc low ",2i4,a," ",4i5)')
     *            ipass, itarg, ntstn(itarg), ichgmo(itarg,ichg-1),
     *            ichgmo(itarg,ichg), ichgmo(itarg,ichg2), npairs
                goto 70
              else 
c               else there were not enough neighbors to make an evaluation
c               THEREFORE THE TARGET HAS AN UNSTABLE NETWORK AND
c                 CANNOT BE COMPLETED
                write(6,2000) itarg, ntstn(itarg), outid, nloop, 
     *            chgptstr, numc, tadj(1)
 2000           format(i5, 1x, a6, '-UNSTBL ', a, i2.2, a,
     *            ' Neigh: ', i2, ' Only est: ', f7.2)
c               keep track of unstable dates
                nunstbl(itarg) = nunstbl(itarg) + 1
                iunstbl(itarg,nunstbl(itarg)) = end1
c                iremseg = 1
c                goto 67
                goto 70 
              endif
            endif
                
c           1st - remove both adjustment and trend outliers 
c           2nd - calculate median adjustment
c           filter around the inner-quartile range.....
c             aqscale is used for trimming outliers
c             sqscale is used for significance testing
            if(ihtag .ne. '') then
c             for HOFN - use a closer (95%) clipping for outliers
c             v50d (v19d)
              aqscale = 1.0 * qscale
              sqscale = 1.0 * qscale
            else
c             for normal output - use (99%) clipping for outliers
c             v20d - qscale from command line restored
              aqscale = 1.5 * qscale
              sqscale = 1.0 * qscale
            endif  

c           sort the changepoint adjustments first
            call sort(npairs,tchgs,itind)
            call tukey_med(npairs, tchgs, pct25, pct50, pct75)
            rng = pct75 - pct25
            rnglo = pct25 - (pct50 - pct25) * aqscale
            rnghi = pct75 + (pct75 - pct50) * aqscale
            write(6,'(" TRIM p25, p75, pct50, rng, lo, hi:", 
     *          6f7.2)') pct25,pct75,pct50,rng,rnglo,rnghi
            
c           trim outliers from adjustments
            do i = 1, npairs
              if(tchgs(i) .lt. rnglo .or. tchgs(i) .gt. rnghi) then
                iused(ichgs(itind(i))) = 'X'
              endif
            enddo

            npairs = 0
            do it2pair = 2, nstns 
c             ipair - the index of the "it2pair" station in the "itarg" subnet
              ipair = nindx(itarg, it2pair)
c             if there are no more neighbors, skip out
              if(ipair .eq. 0) goto 59
              if(iused(it2pair) .eq. 'U') then
                npairs = npairs + 1
                tchgs(npairs) = tradj(it2pair,ichg)
                tccor(npairs) = trcor(it2pair,ichg)
                ichgs(npairs) = it2pair
                isegused(it2pair,ichg) = 'U'
                isegused(it2pair,ichg2) = 'U'
              endif
              
              if(iedebug .gt. 0 .and. iused(it2pair) .ne. 'M')
     *          write(6,'(a1,i4,f7.2,2f8.4,f7.2)') iused(it2pair),
     *            it2pair, tradj(it2pair,ichg), trend(it2pair,ichg),
     *            trend(it2pair,ichg2), trcor(it2pair,ichg)
            enddo
          
   59       if(npairs .lt. numclim) then
              if(ipass .lt. npass) then
                write(6,'("Insuff trimmed mean ",2i4,a," ",4i5)')
     *            ipass, itarg, ntstn(itarg), ichgmo(itarg,ichg-1),
     *            ichgmo(itarg,ichg), ichgmo(itarg,ichg2), npairs
                goto 70
              else 
c               else there were not enough neighbors to make an evaluation
c               THEREFORE THE TARGET HAS AN UNSTABLE NETWORK AND
c                 CANNOT BE COMPLETED
                write(6,2001) itarg, ntstn(itarg), outid, nloop, 
     *            chgptstr, numc, tadj(1)
 2001           format(i5, 1x, a6, '-UNSTB2 ', a, i2.2, a,
     *            ' Neigh: ', i2, ' Only est: ', f7.2)
                nunstbl(itarg) = nunstbl(itarg) + 1
                iunstbl(itarg,nunstbl(itarg)) = end1
                iremseg = 1
                goto 67
              endif
            endif  

            arng = 999.
            lpairs = npairs

c           sort wrt the correlation
            call sortall(0, lpairs, tccor, tchgs, ichgs)
c           remove all changepoints with correlations lower than THRESCOR
            threscor = 0.0
            do ipair = 1, lpairs
              if(tccor(ipair) .le. threscor) then
                lpairs = ipair - 1
                goto 61
              endif
            enddo    

c           sort adjustments to get 25-50-75 percentiles to
c           generate adjustment and 95% conf using trimmed median array
   61       call sortall(1, lpairs, tchgs, tccor, ichgs)
            call tukey_med(lpairs, tchgs, spct25, spct50, spct75)
            srng = spct75 - spct25
            srnglo = spct25 - (spct50 - spct25) * sqscale
            srnghi = spct75 + (spct75 - spct50) * sqscale
            p00 = tchgs(1)
            p100 = tchgs(lpairs)
            if(iedebug .gt. 0) then
              do i = 1, lpairs
                write(6,'("Adj: ", i4, f7.3, i4, f7.2)')
     *             i,  tchgs(i), ichgs(i), tccor(i)
              enddo
            endif  
            if(iedebug .ge. 1)
     *        write(6,'("Amp ", i5," p25, p75, pct50, rng, lo, hi:", 
     *           8f7.2)') lpairs, spct25, spct75, spct50, srng, srnglo,
     *           srnghi, p00, p100
c           hold on to the lowest inner quartile range
c            if(srng .lt. arng) then
              pct25 = spct25
              pct75 = spct75
              tpct50 = spct50
              pct00 = p00
              pct100 = p100
              arng = srng
              rnglo = srnglo
              rnghi = srnghi
              npairs = lpairs
c            endif  
            
c           test to see if the lowest corr can be removed and loop around
            if(lpairs .gt. numclim) then
c             sort wrt the correlation
              call sortall(0, lpairs, tccor, tchgs, ichgs)
c             remove all corr less than threshold (0.0)
              if(tccor(lpairs) .lt. corthres) then
                do i = lpairs, numclim, -1
                  if(tccor(i) .gt. corthres) goto 61
                  lpairs = lpairs - 1
                enddo
              endif 
            endif
            numadj = npairs
            
c           ensure offset is significant for chgpt adj
            if(rnglo * rnghi .gt. 0.0) then
              procstr = 'CONSHF'
              sigadj = tpct50
            else
c             too close to zero
              sigadj = 0.0
              procstr = 'ZERSHF'
            endif  
            if(istep .gt. 1) ichgout(itarg, ichg) = ipass
            adj(itarg, ichg) = sigadj
            std(itarg, ichg) = arng * sqscale
            num(itarg, ichg) = numadj
            ispanq(itarg, ichg2) = iend2 - ibeg2 + 1
            write(6,1100) itarg, ntstn(itarg), procstr, outid, nloop,
     *        chgptstr, sigadj, sthreslo, sthreshi, npairs, 
     *        ichgout(itarg, ichg), nhits(itarg,ichgmo(itarg,ichg)),
     *        ispanq(itarg,ichg2)
 1100       format(i5, 1x, a6, '-', a6, ' ',a, i2.2, a,
     *        ' AVG ADJ: ', 3f7.2, 3i4, i5)
            siglast = tpct50
            goto 70

c           remove the indicated segment and continue
c             1) target segment too short
c             2) Unstable segment (too few trimmed mean)
   67       print *,ntstn(itarg),' ',chgptstr,' Remove chgpt inloop: ',
     *        iremseg, ichgmo(itarg,ichg)
            if(iremseg .eq. 1) then
              irem1 = beg1
              irem2 = end1
            else
              irem1 = beg2
              irem2 = end2
            endif
            do irem = irem1, irem2
              call imo2iym(iy,im,irem)
              temp(itarg,iy,im) = amiss
            enddo  

c           remove changepoint in loop:
c             1) target segment too short
c             2) Unstable segment (too few trimmed mean)
            nhits(itarg, end1) = -1 * nhits(itarg, end1)
            do ic = ichg, nchg - 1
              ichgout(itarg,ic) = ichgout(itarg,ic+1)
              ichgmo(itarg,ic) = ichgmo(itarg,ic+1)
              num(itarg,ic) = num(itarg,ic+1)
              adj(itarg,ic) = adj(itarg,ic+1)
              std(itarg,ic) = std(itarg,ic+1)

              if(ic .le. nchg-2) then
                ispanq(itarg,ic+1) = ispanq(itarg,ic+2)
              endif  

              do istn = 1, nstns
                if(ic .le. nchg-2) then
                  trend(istn,ic+1) = trend(istn,ic+2)
                  spanob(istn,ic+1) = spanob(istn,ic+2)
                endif  
                tradj(istn,ic) = tradj(istn,ic+1)
                tqmin(istn,ic) = tqmin(istn,ic+1)
                itrqt(istn,ic) = itrqt(istn,ic+1)
                isegused(istn,ic) = isegused(istn,ic+1)
              enddo  

            enddo
            nchgin(itarg) = nchgin(itarg)-1
            nchg = nchg - 1 
            goto 70
            
c           set removal of chgpt ONLY - AFTER loop is finished
c             1) non-significant changepoint
c             2) unsupported metadata
   69       iremchg(itarg,ichg) = 1
            print *,ntstn(itarg),' ',chgptstr,' Set Remove chgpt ',
     *        ichgmo(itarg,ichg)

   70       ichg = ichg - 1
          enddo ! ICHG loop - end of station chgpt list

c         No slope work after v20d

c         List the confirmfilt output changepoints
   75     if(iedebug .ge. 1 .and. istep .ge. 2) then
            do ichg = 1, nchg
              if(ichgmo(itarg,ichg) .le. 0) then
                print *,' ichgmo error: ', itarg, ntstn(itarg), ichg, 
     *            iy, im, ichgmo(itarg, ichg)
                nchg = 1
                goto 77
              else
                call imo2iym(iy,im,ichgmo(itarg, ichg))
                write(6,'(i3,i5,1x,a," Estamt chgout: ",2i5,i2.2,3i6,
     *            2f7.2)') ipass, itarg, ntstn(itarg), ichg, 
     *            iy, im, ichgmo(itarg,ichg), 
     *            nhits(itarg,ichgmo(itarg,ichg)), ichgout(itarg,ichg),
     *            adj(itarg,ichg),std(itarg,ichg)
              endif
            enddo  
          endif
          
   77   enddo ! end of station loop  
        
c       Remove accumulated non-sig chgpts & unsupported metadata here
c           remove changepoint in loop:
c             1) target segment too short
c             2) Unstable segment (too few trimmed mean)
        do k = 1, numsubs
c         itarg is the index of the sub-network (candidate)
          itarg = nindx(k,1)
          nchg = nchgin(itarg)
          lrem = 0
          do ichg = nchg-1,2,-1
            irem = 0
            if(istep .eq. 1) then
c             for metadata loop check the removal array
              if(iremchg(itarg,ichg).eq.1) irem = 1
            else if(adj(itarg,ichg) .eq. 0.0 .and. 
     *        ichgout(itarg,ichg) .gt. 0) then
c             for all others, the adj must be 0.0
              irem = 1
            endif  

            if(lrem .eq. 1) then
              ichgout(itarg,ichg) = 0
              lrem = 0
            endif

            if(irem .eq. 1) then
              lrem = 1
              iremchg(itarg, ichg) = 0
              end1=ichgmo(itarg,ichg)
              call imo2iym(iy,im,end1)
              print *,ntstn(itarg),itarg,' Remove chgpt ',ichg,iy,im

              nhits(itarg, end1) = -1 * nhits(itarg, end1)
c             reset the chgpt for recomputation
              ichgout(itarg,ichg) = 0
              do ic = ichg, nchg - 1
                ichgmo(itarg,ic) = ichgmo(itarg,ic+1)
                num(itarg,ic) = num(itarg,ic+1)
                adj(itarg,ic) = adj(itarg,ic+1)
                std(itarg,ic) = std(itarg,ic+1)
                if(ic .le. nchg-2) then
                  ispanq(itarg,ic+1) = ispanq(itarg,ic+2)
                endif  
              enddo
              nchgin(itarg) = nchgin(itarg)-1
              nchg = nchg - 1
            endif
            
          enddo  
        enddo  

c       loopback for "npass" testing
        if(istep .eq. 3) goto 27
   87 enddo ! end of istep loop
      
c          ------------------ End of section determining ------------------
c          -------------------------- adjustment  -------------------------
c          ----------------------------------------------------------------
          
   90 ipass = npass + 1
          
c       ----------- End of section determining whether -----------------
c       ---------------- an adjustment can be made ---------------------
c       ----------------------------------------------------------------
          
c   95 enddo ! end of all candidates (the network)

      print *,' ------------ Output Adjustments ------------'
      
c     if HOFN output is enabled, open the various output files
      if(ihtag .ne. '') then
c       Cand-Net monthly anom data
        fname = odir(1:lnblnk(odir)) // 'chgpnt_series_' // 
     *    ihtag(1:lnblnk(ihtag)) // icelem // '.dat'
        open(36, file=fname, status='unknown', err=110)
c       Cand Chgpnt date/offset data
        fname = odir(1:lnblnk(odir)) // 'chgpnt_date_' // 
     *    ihtag(1:lnblnk(ihtag)) // icelem // '.dat'
        open(39, file=fname, status='unknown', err=110)
c       Cand-Net mean segment normalized diff data
        fname = odir(1:lnblnk(odir)) // 'chgpnt_mean_' // 
     *    ihtag(1:lnblnk(ihtag)) // icelem // '.dat'
        open(46, file=fname, status='unknown', err=110)
c       Cand-Net raw temp data
        fname = odir(1:lnblnk(odir)) // 'chgpnt_temp_' // 
     *    ihtag(1:lnblnk(ihtag)) // icelem // '.dat'
        open(47, file=fname, status='unknown', err=110)
      endif  
      
c     finished with entire network adj passes - write out results
      do k = 1, numsubs
c       itarg is the index of the sub-network (candidate)
        itarg = nindx(k,1)
        nchg = nchgin(itarg)
 
c       initialize output arrays
        iflg = iachar('A')
        do iy = begyr, endyr
          do im = 1,13
            adjflag(iy,im) = ' '
            outtemp(iy,im) = amiss
            adjtemp(iy,im) = amiss
            contemp(iy,im) = amiss
          enddo
        enddo

        if(nunstbl(itarg) .gt. 1) then
          do iun = 1, ninh
            nundat(iun) = 0
          enddo  
          nunstbl(itarg) = nunstbl(itarg) + 1
          iunstbl(itarg,nunstbl(itarg)) = 1
c         figure out which segment to remove (the one with less data)
c            add 'g' (gap) to the sement that will be removed in Fillin
          do iun = 1, nunstbl(itarg)-1
c           find the number of months in each segment
            do imo = iunstbl(itarg,iun), iunstbl(itarg, iun+1), -1
              call imo2iym(iy,im,imo)
              if(filtemp(itarg,iy,im).gt.amiss+1.) then
                nundat(iun) = nundat(iun) + 1
              endif
            enddo
          enddo
          print *,itarg, ' UNSTABLE imo nval'
          lonun = 9999
          do  iun = 1, nunstbl(itarg)-1
            write(6,'(3i5)') iun,iunstbl(itarg,iun),nundat(iun)
            if(nundat(iun) .lt. lonun) then
              loun = iun
              lonun = nundat(iun)
            endif  
          enddo
          print *,' Remove segment: ', loun, iunstbl(itarg,loun), 
     *      iunstbl(itarg, loun+1)
          do imo = iunstbl(itarg,loun), iunstbl(itarg, loun+1), -1
            call imo2iym(iy,im,imo)
            tflg(itarg,iy,im) = 'g'
          enddo
        endif  
              
c       Adjust the candidate series output use filtemp
c          - it has all of incoming data
        sumchg = 0.0
        sumcon = 0.0
c       count the adjusted chgpt - going backwards
        jadj = 0
        aflg = ' '
c       set the end of the first segment at the end POR
        end2 = ichgmo(itarg, nchg)
        ichg = nchg-1
        do while(ichg .ge. 1)
          jadj = jadj + 1
          if(ichg .eq. 1) then
            beg2 = ichgmo(itarg,ichg)
          else  
            beg2 = ichgmo(itarg,ichg) + 1
          endif
          print *, 'imo: ', itarg, ichg, beg2, end2
          delseg = 0
          do imo = end2, beg2, -1
            call imo2iym(iy,im,imo)
            if(filtemp(itarg,iy,im).gt.amiss+1.) then
              totchg = sumchg + delseg
              adjtemp(iy,im) = totchg
              contemp(iy,im) = sumcon
              outtemp(iy,im) = filtemp(itarg,iy,im) - totchg
              adjflag(iy,im) = aflg
c              print *, imo, end2-imo, delseg, sumchg, totchg
            endif
          enddo
          if (itimeres .eq. 0) then
c           for annual data - put adjustments into monthly data
            iy1 = beg2 + begyr
            iy2 = end2 + begyr
            print *, 'mth: ', itarg, ichg, iy1, iy2
            delseg = 0
            do iy = iy2, iy1, -1
              do im = 12, 1, -1
                if(filtemp(itarg,iy,im).gt.amiss+1.) then
                  yrchg = sumchg + delseg 
                  adjtemp(iy,im) = yrchg
                  contemp(iy,im) = sumcon
                  outtemp(iy,im) = filtemp(itarg,iy,im) - yrchg
                  adjflag(iy,im) = aflg
c                  print *, iy, im, delseg, sumchg, yrchg
                endif
              enddo
            enddo
          endif
          
          write(6,'("Adj write:",a,5i5,2f7.2,a)') ntstn(itarg), nchg,
     *       beg2, end2, ichg, jadj, sumchg, sumcon, aflg
                
          aflg = achar(iflg + jadj - 1)
          sumchg = totchg + adj(itarg,ichg)
          sumcon=sqrt(sumcon*sumcon+std(itarg,ichg)*std(itarg,ichg))
          end2 = ichgmo(itarg, ichg)
          ichg = ichg - 1  
        enddo
        
        if(ihtag .ne. '') then
          call hofnout(itarg, ntstn, nchg, ichgmo)
        else  
          call writsta(itarg, ntstn(itarg), outtemp, adjtemp,
     *      contemp, adjflag, otag, idunit)
        endif
      enddo ! end of station write output loop
      
      goto 120

c     Error for HOFN files
  110 call perror(' Cannot open HOFN data file: ' // fname)
c     make sure to close the HOFN units!!!
  120 if(ihtag .ne. '') then
        close(36)
        close(39)
        close(46)
        close(47)
      endif
     
      return
      end

c =======================================================================

      subroutine stdmth(itarg, avgstd)
      INCLUDE 'inhomog.parm.mthly.incl'
      INCLUDE 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'
      
      real rSum(13), rSum2(13), stddev(13), rNum(13)
      
      do im = 1, 13
        rSum(im) = 0.0
        rSum2(im) = 0.0
        rNum(im) = 0.0
      enddo 
      
      do iy = begyr, endyr
c       if annual time resolution (lomth,himth=13)
        do im = lomth, himth
          if(temp(itarg,iy,im).gt.amiss+1.) then
            rSum(im) = rSum(im) + temp(itarg,iy,im)
            rSum2(im) = rSum2(im) + temp(itarg,iy,im)*temp(itarg,iy,im)
            rNum(im) = rNum(im)+1.0
          endif
        enddo
      enddo
        
      sumstd = 0.0  
      numstd = 0
c     if annual time resolution (lomth,himth=13)
      do im = lomth, himth
        if(rNum(im) .gt. 1.0) then
          stddev(im) = sqrt((rsum2(im) - 
     *      (rsum(im) * rsum(im) / rnum(im))) / (rnum(im) - 1))
          sumstd = sumstd + stddev(im)
          numstd = numstd + 1
c          print *,'STDDEV: ',im,stddev(im),rsum(im),rsum2(im),rnum(im)
        else
          stddev(im) = amiss
        endif
      enddo

      if(numstd .gt. 0) then
        avgstd = sumstd / numstd
c        print *,' AVGSTD: ', avgstd, numstd
      else
        avgstd = amiss
      endif
      return
      end

C     ************************************************************

      subroutine sort(n,arr, ndx)

c     Shell's method.

c     n is the number of elements to sort in array ARR, 
c     ndx is an associated array that is sorted along with ARR
      real arr(n)
      integer ndx(n)

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
            nvit=ndx(i)

            arr(i)=arr(l)
            ndx(i)=ndx(l)

            arr(l)=t
            ndx(l)=nvit

            i=i-m
            if(i.ge.1)go to 553
          endif
        enddo
      enddo

 900  continue
      return

      end
C     ************************************************************

      subroutine sortall(iscend, n, arr1, arr2, ndx)

c     Shell's method.

c     iscend is the option for direction of sort
c         0 = descending
c         1 = ascending
c     n is the number of elements to sort in array ARR1
c     arr2, and ndx are sorted along with ARR1
      real arr1(n), arr2(n)
      integer ndx(n)

      if(iscend .ne. 0 .and. iscend .ne. 1) then
        print *,' Unknown sort direction: ',iscend
        stop
      endif  
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
c         descending sort (iscend == 0)
c         ascending sort (iscend == 1)
          if((iscend .eq. 0 .and. arr1(l).gt.arr1(i)) .or.
     *       (iscend .eq. 1 .and. arr1(l).lt.arr1(i)))then
            t1=arr1(i)
            t2=arr2(i)
            nvit=ndx(i)

            arr1(i)=arr1(l)
            arr2(i)=arr2(l)
            ndx(i)=ndx(l)

            arr1(l)=t1
            arr2(l)=t2
            ndx(l)=nvit

            i=i-m
            if(i.ge.1)go to 553
          endif
        enddo
      enddo

 900  continue
      return
      end
        
c     =======================================================================

      subroutine correl(nx,x,y,ccoef,amiss)
c********************************************************************
c
c     subroutine correl computes the pearson correlation coefficent
c     between a candidate station and its paired neighbor.
c
c     the mean, stnd dev, and cross-products are computed only
c     for those months when both candidate and neighbor have good data.  
c
c     each correlation coefficient applies to all months of the year.
c
c                               -on input-
c
c       nx      number of months of good and missing data in each 
c               time series of monthly departures.
c
c       x(nx)   the array of monthly departures for the candidate
c
c       y(nx)   the array of monthly departures for the paired neighbor
c
c                               -on output-
c
c       ccoef   the (ns-1 x 1) array of candidate-neighbor pearson correlation
c               coefficients
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
        if (x(im).gt.amiss+1..and.y(im).gt.amiss+1.) then
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
        
c     =======================================================================

      subroutine tukey_med(nchgs, tchgs, pct25, pct50, pct75)
c     trimmed averages using Tukey's median method

      INCLUDE 'inhomog.parm.mthly.incl'
      real tchgs(nstns), tempchg(nstns)
      
      do i = 1, nchgs
        tempchg(i) = tchgs(i)
      enddo
      ntemp = nchgs  
      
c     calculate the median
      if(mod(ntemp,2) .eq. 1) then
        iy50 = ntemp/2 + 1
        pct50 = tempchg(iy50)
      else
        iy50 = ntemp/2  
        pct50=(tempchg(iy50)+tempchg(iy50+1))/2.
c       if median is an average, add back into the obs
        do i = ntemp, iy50+1, -1
          tempchg(i+1) = tempchg(i)
        enddo
        tempchg(iy50+1) = pct50
        ntemp = ntemp + 1
        iy50 = iy50 + 1
      endif
c     calculate the lower quartile (include median)
      if(mod(iy50,2) .eq. 1) then
        iy25 = iy50/2 + 1
        pct25 = tempchg(iy25)
      else
        iy25 = iy50/2
        pct25=(tempchg(iy25)+tempchg(iy25+1))/2.
      endif    
c     calculate the upper quartile (include median)
      nhigh = ntemp - iy50 + 1
      if(mod(nhigh,2) .eq. 1) then
        iy75 = nhigh/2 + iy50
        pct75 = tempchg(iy75)
      else
        iy75 = nhigh/2 + iy50 - 1
        pct75=(tempchg(iy75)+tempchg(iy75+1))/2.
      endif   
      return
      end    
 
c     =======================================================================

      subroutine frstdif (nx,x,y,amiss,nout,x1d,y1d)
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
c       x       the (nx) monthly mean temperatures
c
c       amiss   defined missing value 
c
c                               -on output-
c
c       y       the (nx) first differences
c
c******************************************************************
      real x(nx), y(nx), amiss
      real x1d(nx), y1d(nx)
      kount = 0
      nout = 0

c     initialize
      do im = 1, nx
        x1d(im) = amiss
        y1d(im) = amiss
      enddo

c
c     do first difference filter accounting for missing data
c
      do im = 1,nx
        if (x(im).gt.amiss+1 .and. y(im).gt.amiss) then
          kount = kount + 1
c          print *, im, kount, x(im), y(im)
          if (kount.gt.1) then
            nout = nout + 1
            x1d(nout) = (x(im) - x(lim))/2
            y1d(nout) = (y(im) - y(lim))/2
c            write(6,'(3i5,6f7.2)') im,kount,nout,x1d(nout),x(im),
c     *        x(lim),y1d(nout),y(im),y(lim)
          endif  
          lim = im
        endif
      enddo  
 99   return
      end     
        
c     =======================================================================

      subroutine mmts_adj(mmdates)
      
      INCLUDE 'inhomog.parm.mthly.incl'
      include 'inhomog.comm.mthly.incl'
      include 'inhomog.restart.mthly.incl'
      integer mmdates(maxstns,2)
      delta = 0.0      
            
C     MMTS IMO INDICES ARE IN MONTHLY RESOLUTION
      ltimeres = itimeres
      itimeres = 1
      do itarg = 1, numsubs
        if(mmdates(itarg,1) .ne. 9999) then
          imo1 = mmdates(itarg,1)
c         This stn has MMTS - adjust till not MMTS
          imo2 = mmdates(itarg,2) - 1
          if(imo2 .eq. 9998) imo2 = numyr * 12
          call imo2iym(ly,lm,imo1)
          do imo = imo1, imo2
            call imo2iym(iy,im,imo)

            if(iy .ne. ly) then
              if(ndelta .eq. 12) then
                ydelta = sdelta / ndelta
                temp(itarg,ly,13) = temp(itarg,ly,13) + ydelta
c                print *,'Stn yr MMTS:',itarg,ly,ydelta
              endif
              sdelta = 0.0
              ndelta = 0
              ly = iy
            endif  
              
            if(temp(itarg,iy,im) .gt. amiss+1.0) then
c             if max temp then for spring +.6 else +.7
              if(inel .eq. 1) then
                if(im .ge. 3 .and. im .le. 5) then
                  delta = 0.6
                else
                  delta = 0.7
                endif
              else if(inel .eq. 2) then    
c             if min temp then for summer/fall -.6 else -.4
                if(im .ge. 6 .and. im .le. 11) then
                  delta = -0.6
                else
                  delta = -0.4
                endif
              endif
              temp(itarg,iy,im) = temp(itarg,iy,im) + delta
              sdelta = sdelta + delta
              ndelta = ndelta + 1
            endif
          enddo
        endif    
      enddo
c     Set time resolution back to incoming
      itimeres = ltimeres
      return
      end
