#!/bin/bash

export CONFIL=$1
echo "Configuration file  is:" $CONFIL
if [ ! -f "$CONFIL" ]; then
    echo "Configuration file does not exist"
    exit
fi
ls -l $CONFIL
pwd
cp $CONFIL .
ls -l

/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/pbpb/mainAA > pythiaAA.log
