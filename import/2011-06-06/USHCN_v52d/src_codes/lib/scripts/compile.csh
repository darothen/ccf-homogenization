#!/bin/csh

# upgraded parameter list to include ANNUAL as well as MONTHLY - see inhomog.parm.system.....incl

@ reqparm = 2

if($#argv < $reqparm) then
  echo "Compile the USHCN inhomog application"
  echo "  USAGE: compile.csh comptype net_def"
  echo "     where"
  echo "           comptype is the compilation options for debugging or fast"
  echo "               Current comptype are:"
  echo "                DBUG - for symbolic links for debugging"
  echo "                FAST - for fast execution options"
  echo "           net_def is the network name used to define "
  echo "               the Network dependent array settings"
  echo "               Current net_def are:"
  echo "                21_21"
  echo "                USHCN"
  echo ""
  exit
endif  

echo "Entering :" $0 $*

set comptype = $1
set net_def = $2
set USHCNBASE = $3
setenv $USHCNBASE

# Update for each USHCN v2 update
set version = UCPM24aM21fC6bL1Si6cR5aH1
set time_def = MLY
set SCRIPTS = $USHCNBASE/src_codes/lib/scripts
set SRC = $USHCNBASE/src_codes/lib/pha_src
set BIN = $USHCNBASE/bin

set cwd = `pwd`
cd $SRC

rm inhomog.parm.mthly.incl >& /dev/null
ln -s ./parm_includes/inhomog.parm.MTHLY.${net_def}.incl inhomog.parm.mthly.incl

set binfile = $BIN/${version}.${comptype}.${time_def}.${net_def}

echo "Begin $binfile Compilation"
$SCRIPTS/$version.csh $binfile $comptype
echo "End $binfile Compilation"

echo "Leaving :" $0
cd $cwd

