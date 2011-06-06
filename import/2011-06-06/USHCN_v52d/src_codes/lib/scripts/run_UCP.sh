#!/bin/sh

# set number of required parameters
reqparm=`expr 6`

if [ $# -lt $reqparm ]
then
  echo "Not all command line input defined. Please enter"
  echo '  $1 Conus time stamp directory to use'
  echo '  $2 is the current UCP routine (i.e. UCPM17gM21dC5iL1Si4eR3f)'
  echo '     FAST.MLY.USHCN compile version assumed'
  echo '  $3 is the element (max, min)'
  echo '  $4 is the input process id (i.e. tob)'
  echo '  $5 is the output process id (i.e. 52d)'
  echo '  $6 <optional> is the Journal process id (i.e. g31)'
  exit
fi 

conus=$1
# runid=$1.DBUG.MLY.USHCN
runid=$2.FAST.MLY.USHCN
inelem=$3
iproc=$4
oproc=$5
if [ $# -eq $reqparm ]
then
  jproc="X"
else  
  jproc=$6
fi  

# define element index from input element
if [ $inelem == "max" ] 
then
  ielem=1
elif [ $inelem == "min" ]
then
  ielem=2
else
  echo "Unknown element: " $inelem
  exit
fi

basefile=~/USHCN/v2/conus/auto_update/$conus

dtag=`date +%y%m%d%H%M`

# Correlation network file
corrfile=$basefile/meta/ushcn_cdmp.corr.us
# Station list
metafile=$basefile/meta/ushcn_cdmp.meta.us

# data directories
incand=$basefile/monthly_coop
outcand=$incand
incoop=$basefile/monthly_coop
outcoop=$incoop

# ensure the output directories are created
if [ ! -d $outcoop/WMs.$oproc ]; then mkdir -p $outcoop/WMs.$oproc; fi 
if [ ! -d $outcoop/FLs.$oproc ]; then mkdir -p $outcoop/FLs.$oproc; fi 

outfile=$basefile/output/$runid.$dtag.$inelem.$ireg.$oproc.out
if [ $jproc == "X" ]
then
  echo "FULL UCP: $runid oproc:$oproc i:$ielem elem:$inelem dtag:$dtag"
    nice /home/cwilliam/USHCN/inhomog/source_expand/$runid -Q 1.46 -S 18 -s 5 -P -o $oproc -d 1 \
      -T 100 -l -e $ielem -p $iproc -c 0 -t 1 -q $iproc  -n $corrfile -m  $metafile -O . \
      -C $incand $outcand -N $incoop $outcoop -u $oproc.$dtag.$ielem.$ireg > $outfile
  else
    list=`ls $jproc.*.$ielem.$ireg.0001_1_01_3 2> err`
    for jrnl in $list
    do
      echo "JRNL UCP: $runid oproc:$oproc i:$ielem elem:$inelem dtag:$dtag jrnl:$jrnl"
      nice /home/cwilliam/USHCN/inhomog/source_expand/$runid -Q 1.46 -S 18 -s 5 -P -o $oproc -d 1 \
        -T 100 -l -e $ielem -p $iproc -c 0 -t 1 -q $iproc -n $corrfile -m  $metafile -O . \
        -C $incand $outcand -N $incoop $outcoop -O . -j $jrnl > $outfile
    done
  fi  
  echo "FULL FILL: $runid oproc:$oproc i:$ielem elem:$inelem dtag:$dtag"
  nice /home/cwilliam/USHCN/ushcn_lib/src/ushcn_fill_2004.v4p -d 1 -u 0 -e $ielem -p WMs.$oproc \
    -q WMs.$oproc -i WMc.$oproc -o FLs.$oproc -j FLc.$oproc -n $corrfile -c $metafile  -C $outcand \
    -N $outcoop > $basefile/output/fill2004.$ireg.$inelem.$dtag.$oproc.out
done    
