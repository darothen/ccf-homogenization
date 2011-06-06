#!/bin/sh

# set number of required parameters
reqparm=`expr 4`

if [ $# -lt $reqparm ]
then
  echo "Not all command line input defined. Please enter"
  echo '  $1 - Elem ( max min avg )'
  echo '  $2 is the output process id (i.e. 52d)'
  echo '  $3 is metafile name in realworld/meta'
  echo '  $4 is the USHCN Base directory'
  exit
fi

# set the date tag for this run
dtag=`date +%Y%m%d%H%M`

# set meteorological element (max, min, avg - temperature: pcp - precipitation)
elem=$1
ver=$2
mfile=$3
USHCNBASE=$4

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
if [ ! -d $datadir/FLs.$ver ]; then mkdir -p $datadir/FLs.$ver; fi 

# MODIFY - define the log directory path for all text output
logdir=$USHCNBASE/realworld/log

# MODIFY - define the bin directory path to the compiled codes
bin=$USHCNBASE/bin

echo "Fillin UCP series for $elem.$dtag.$ver"  
nice $bin/ushcn_fill_2004.v4p -d 0 -u 0 -e $iel \
  -p WMs.$ver -q WMs.$ver -i WMc.$ver -o FLs.$ver -j FLc.$ver \
  -n $corrfile -c $metafile  -C $datadir -N $datadir \
  > $logdir/fill2004.$elem.$dtag.$ver.out
