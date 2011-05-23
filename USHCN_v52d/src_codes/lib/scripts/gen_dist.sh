#!/bin/sh

# set number of required parameters
reqparm=`expr 3`

# if the number of required parameters are not entered, display usage message
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

# MODIFY - define the metadir path for preparatory files
metadir=$USHCNBASE/realworld/meta
# set station list file (input)
metafile=$metadir/$mfile
# set station distance neighborhood file (output)
distfile=$metadir/$mfile.dist

# MODIFY - define the log directory path for all text output
logdir=$USHCNBASE/realworld/log

# MODIFY - define the bin directory path to the compiled codes
bin=$USHCNBASE/bin

echo "Generate the distance network file"
$bin/ushcn_dist_2004.v3 -u 0 -n $metafile -h -c $metafile -o $distfile \
  > $logdir/dist.$elem.$dtag.log
