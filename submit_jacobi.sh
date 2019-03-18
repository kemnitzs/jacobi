#!/bin/bash
BASE=$(pwd)

echo $BASE
RUNFOLDER=run/
mkdir $RUNFOLDER 
cd $RUNFOLDER

cp $BASE/bin_jacobi .
# run script
./bin_jacobi

cd $BASE
# plot results
STARTTIME=$(date +%s)
python3 plot_out.py $RUNFOLDER
ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete to create all images..."
