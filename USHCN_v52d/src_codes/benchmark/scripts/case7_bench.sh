#!/bin/sh

# set number of required parameters
reqparm=`expr 1`

if [ $# -lt $reqparm ]
then
  echo "Not all command line input defined. Please enter"
  echo '  $1 is the current UCP routine (i.e. UCPM24pM21qC6pL1Si6pR5p.lf95)'
  echo '     FAST.MLY.21_21 version - FOR USE WITH BENCHMARK CASE ONLY'
  exit
fi

runid=$1.DBUG.MLY.21_21
# runid=$1.DBUG.MLY.USHCN
# runid=$1.FAST.MLY.21_21

inelem="avg"
ielem=3
iproc=raw

# directory with APH executable binary RUNID
bin=$USHCNBASE/bin

# base directory with benchmark data
base=$USHCNBASE/src_codes/benchmark

# directory with corr/meta/input data files
corrdir=$base/corr
  
# use the benchmark metadata that has been randomly clipped (50%)
metafile=$corrdir/case7.meta.clip50
# 1st difference correlated network file
corrfile=$corrdir/case7.corr.ann
# benchmark formatted data file
mattdata=$corrdir/case7.data.txt
  
# directory with output log files
logdir=$base/output
# directory with input/output data directories
datadir=$base/monthly

# input data directory
indata=$datadir/$iproc

oproc=52d
# ensure the output directories are created
outaph=$datadir/WMs.$oproc
if [ ! -d $outaph ]; then mkdir -p $outaph; fi 

dtag=`date +%y%m%d%H%M`
      
logfile=$logdir/case7_bench.$dtag.log
echo "CASE7_BENCH FULL UCP: $runid oproc:$oproc dtag:$dtag"
nice $bin/$runid -Q 1.46 -S 18 -s 5 -P -o $oproc -d 1 -T 100 -l -e 3 \
  -R -21 -p $iproc -q $iproc -c 0 -t 1 -n $corrfile -C $datadir $datadir \
  -r $mattdata -m $metafile -N $datadir $datadir -O . > $logfile

