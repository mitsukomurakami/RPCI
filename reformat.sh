#!/usr/bin/bash
# reformat.sh

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./reformat.sh <patientID>"
    exit 1
fi

pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ023532-Paragh

JOBDIR=$PREFIX/projects/$PID/$pid
DATADIR=$JOBDIR/output/variants2 

SLIST=`ls $DATADIR`

for sid in $SLIST; do

echo Processing: $sid
ifile=`ls $DATADIR/$sid/*bambino`
ofile=${ifile}.sj

if [ -f "$ifile" ]; then
  perl reformat_bambino2SJ.pl -i $ifile -o $ofile -sample $sid 
else
  echo $sid: bambino file does not exist
fi

head $ofile

ls $DATADIR/$sid

done
