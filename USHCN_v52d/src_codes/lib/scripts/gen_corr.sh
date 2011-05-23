#!/bin/sh

# set number of required parameters
reqparm=`expr 3`

if [ $# -lt $reqparm ]
then
  echo "Not all command line input defined. Please enter"
  echo '   $1 - Elem ( max min avg )'
  echo '   $2 - metafile name in realworld/meta'
  echo '   $3 - USHCN Base directory'
  exit
fi

# set the date tag for this run
dtag=`date +%Y%m%d%H%M`

# set meteorological element (max, min, avg - temperature: pcp - precipitation)
elem=$1
mfile=$2
USHCNBASE=$3

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
rm -f $corrfile

# MODIFY - define the input raw data file directory path
datadir=$USHCNBASE/realworld/monthly

# MODIFY - define the log directory path for all text output
logdir=$USHCNBASE/realworld/log

# MODIFY - define the bin directory path to the compiled codes
bin=$USHCNBASE/bin

echo "Generate the correlation network file"
$bin/ushcn_corr_2004.v3 2009 $iel -u 0 -p raw -d $distfile -m $metafile \
  -o $corrfile -C $datadir -N $datadir > $logdir/corr.$elem.$dtag.log
