#!/bin/bash

JOB_ID=$1

export SETUPFIL_ORIG="PbPbcharm4alice.cmnd"
export SETUPFIL="PbPbcharm4alice_seed.cmnd"
export SCRIPTN="mainAA.cc"
export NEVT=5
export OUTFILE="charmHadr.txt"

# #mkdir testDirCreate
# #man tar >> outSeed.log
tar xzvf pythia8311.tar.gz >> outSeed.log
ls >> outSeed.log
 export PYTHIA8="$PWD/pythia8311"
 export PYTHIA8DATA="${PYTHIA8}/share/Pythia8/xmldoc"
cd pythia8311
./configure --with-root-bin=$ROOTSYS/bin/ --with-root-include=$ROOTSYS/include/ --with-root-lib=$ROOTSYS/lib/ &> outConfig.log
# --with-root-bin=/usr/bin/
# --with-root-include=/usr/include/root
# --with-root-lib=/usr/lib64/root
echo "pythia8: $PYTHIA8" > envVar.log
echo "pythia8data: $PYTHIA8DATA" >> envVar.log
echo "ROOTSYS: $ROOTSYS" >> envVar.log
echo "ROOT_USE: $ROOT_USE" >> envVar.log

make &> outMake.log
cd examples 
make mainAA &> outMakeMacro.log

SEED=$((50000000 + $JOB_ID))
cp $SETUPFIL_ORIG $SETUPFIL
sed -i "s/Main:numberOfEvents = 50/Main:numberOfEvents = $NEVT/g" $SETUPFIL
sed -i "s/Random:seed = 1000000/Random:seed = $SEED/g" $SETUPFIL

./mainAA ./$SETUPFIL ./$OUTFILE &> ./job.log

cd ../..

mv pythia8311/examples/$OUTFILE .
mv pythia8311/examples/job.log .
mv pythia8311/outConfig.log .
mv pythia8311/outMake.log .
mv pythia8311/examples/outMakeMacro.log .
cp pythia8311/Makefile.inc .
cp pythia8311/envVar.log .
rm -r pythia8311


