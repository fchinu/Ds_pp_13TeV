#!/bin/bash

export SETUPFIL="PbPbcharm4alice.cmnd"
export SCRIPTN="mainAA"
export NCORES=40
export NEVT=5
export FIRSTJOB=1290
export LASTJOB=1809
export OUTFILE="charmHadr.txt"

if [ -f $SETUPFIL ]; then
    echo "Generate jobs with configuration from file: " $SETUPFIL
else
    echo "Error: file" $SETUPFIL "does not exist, exiting"
    exit
fi

if [ -f $SCRIPTN ]; then
    echo "Use script: " $SCRIPTN
else
    echo "Error: file" $SCRIPTN "does not exist, exiting"
    exit
fi


for JOB in `seq $FIRSTJOB $LASTJOB`
do
    echo Job $JOB
    OUTDIR=JobPbPb$JOB
    if [ ! -d "$OUTDIR" ]; then
	echo "Output dir does not exist, create it"
	mkdir -p $OUTDIR
    fi

    cd $OUTDIR
    echo "Now in directory:" 
    pwd
    cp ../$SETUPFIL .
    SEED=$((10000000 + $JOB))
    sed -i "s/Main:numberOfEvents = 50/Main:numberOfEvents = $NEVT/g" $SETUPFIL
    sed -i "s/Random:seed = 1000000/Random:seed = $SEED/g" $SETUPFIL
    cp ../$SCRIPTN .
    # if [ -f submitJob ]; then
	# rm submitJob
    # fi
    # touch submitJob
    # echo "Universe   = vanilla" > submitJob
    # echo "" >> submitJob
    # echo "getenv = True" >> submitJob
    # echo "" >> submitJob
    # echo "Executable = /home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/pbpb/"$OUTDIR"/"$SCRIPTN  >> submitJob
    # echo "" >> submitJob
    # echo "Log        = thisJob.log" >> submitJob
    # echo "Output     = thisJob.out" >> submitJob
    # echo "Error      = thisJob.error" >> submitJob
    # echo "Arguments  = /home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/pbpb/"$OUTDIR"/"$SETUPFIL >> submitJob
    # echo "" >> submitJob
    # echo "Queue" >> submitJob
    # #condor_submit submitJob
    cd ..
done

parallel -j $NCORES 'nice -n 5  ./JobPbPb{1}/'$SCRIPTN' ./JobPbPb{1}/'$SETUPFIL' ./JobPbPb{1}/'$OUTFILE' &> ./JobPbPb{1}/job.log' ::: $(seq $FIRSTJOB $LASTJOB)
