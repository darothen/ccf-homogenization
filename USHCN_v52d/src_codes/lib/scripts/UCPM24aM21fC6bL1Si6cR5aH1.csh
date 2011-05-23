#!/bin/csh

@ reqparm = 2

if($#argv < $reqparm) then
  echo "Compile the USHCN inhomog application"
  echo "  USAGE: current_version.csh binfile comptype"
  echo "    where binfile is the output binary executable"
  echo "          comptype is FAST or DBUG"
  echo ""
  exit
endif  

echo "Entering :" $0 $*

set binfile = $1
set comptype = $2
set SCRIPTS = $USHCNBASE/src_codes/lib/scripts
set SRC = $USHCNBASE/src_codes/lib/pha_src

set cwd = `pwd`
cd $SRC

# expand the source code with the include code
rm *.for
foreach file (*.f)
 gawk -f $SCRIPTS/expand_source.awk $file > ${file}or
 echo $file
end

# set the option for FAST or DBUG compilation
if($comptype == "DBUG") then
  set copt = "-g -fbounds-check -fflatten-arrays"
else if($comptype == "FAST") then
#  set copt = "-O2 -funroll-loops"
  set copt = " "
else
  echo "Invalid comtype: $comptype - STOPPING"
  exit  
endif
  
f77 ucpmonthly.v24a.for splitmerge.v21f.for chgptmodels.v6b.for \
  SHAPinp.v6c.for read_write.mthly.v5a.for acovf.for lmbic.for \
  lmdiff.f hofn.for $copt -o $binfile

echo "Leaving :" $0
cd $cwd
