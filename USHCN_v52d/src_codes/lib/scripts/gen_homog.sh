#!/bin/sh

# set number of required parameters
reqparm=`expr 6`

if [ $# -lt $reqparm ]
then
  echo "Not all command line input defined. Please enter"
  echo '  $1 is the current UCP routine, for example:'
  echo '     UCPM24aM21fC6bL1Si6cR5aH1.FAST.MLY.USHCN'
  echo '  $2 is the Element ( max min avg )'
  echo '  $3 is the input process id (i.e. raw)'
  echo '  $4 is the output process id (i.e. 52d)'
  echo '  $5 is metafile name in realworld/meta'
  echo '  $6 is the USHCN Base directory'
  exit
fi

# set the date tag for this run
dtag=`date +%Y%m%d%H%M`

ucp=$1
elem=$2
iproc=$3
oproc=$4
mfile=$5
USHCNBASE=$6

if [ $elem == "max" ]
then
  export iel=1
  export gel="tmax"
  echo "MAX"
elif [ $elem == "min" ]
then
  export iel=2
  export gel="tmin"
  echo "MIN"
elif [ $elem == "pcp" ]
then
  export iel=4
  export gel="prcp"
  echo "PCP"
else
  echo "Unknown element: $elem"
  exit
fi      

# MODIFY - define the metadir path for preparatory files
metadir=$USHCNBASE/realworld/meta
# set station list file (input)
metafile=$metadir/$mfile
# set station distance neighborhood file (input)
distfile=$metadir/$mfile.dist
# set final station neighborhood file (output)
corrfile=$metadir/$mfile.corr
 
# MODIFY - define the input raw data file directory path
datadir=$USHCNBASE/realworld/monthly

# ensure the output directories are created
if [ ! -d $datadir/WMs.$oproc ]; then mkdir -p $datadir/WMs.$oproc; fi 

# MODIFY - define the log directory path for all text output
logdir=$USHCNBASE/realworld/log

# MODIFY - define the bin directory path to the compiled codes
bin=$USHCNBASE/bin

echo "Generate UCP series using: $ucp.$dtag.$elem.$oproc"
nice $bin/$ucp -Q 1.46 -S 18 -s 5 -P -d 1 -T 100 -l -c 0 -t 1 \
  -o $oproc -e $iel -p $iproc -q $iproc -n $corrfile -m $metafile \
  -C $datadir $datadir -N $datadir $datadir -O . \
  -u $oproc.$dtag.$elem > $logdir/$ucp.$dtag.$elem.$oproc.out
