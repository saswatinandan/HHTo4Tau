#!/bin/bash
#
echo "input parameters: cluster, process, run path, out path, sample name" $1 $2 $3 $4
CLUSTER=$1
PROCESS=$2
RUNPATH=$3
OUTPATH=$4


#########For running locally uncommnet this line
#CLUSTER=""
#PROCESS=5
#RUNPATH="/uscms_data/d3/abdollah/Analysis/LQ2016/CMSSW_8_0_11/src/Skim_ggNtuple"
#OUTPATH=""


echo ""
echo "CMSSW on Condor"
echo ""

START_TIME=`/bin/date`
echo "started at $START_TIME"

echo ""
echo "parameter set:"
echo "CLUSTER: $CLUSTER"
echo "PROCESS: $PROCESS"
echo "RUNPATH: $RUNPATH"
echo "OUTPATH: $OUTPATH"

source /cvmfs/cms.cern.ch/cmsset_default.sh


#echo $PWD
#cd $RUNPATH
#eval `scram runtime -sh`
#echo $PWD , "for cmsenv"
#ls 
echo "check whether eos is working here"
echo ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}
echo $PWD , "for job running"

setenv SCRAM_ARCH slc6_amd64_gcc530

xrdcp root://cmseos.fnal.gov//store/user/snandan/CMSSW_8026.tar   .
tar -xf CMSSW_8_0_12.tar
rm CMSSW_8_0_12.tar
cd CMSSW_8_0_12/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh` 
ls
make

## cmsenv is an alias not on the workers

 

###### Copy the neccessary files for running the Skim
#cp $RUNPATH/Skimmer.cc .
#cp $RUNPATH/Skimmer.h .
#cp $RUNPATH/Makefile .
#cp $RUNPATH/InputSamples.txt .





#########  Smaple/Job splitting
SplitingNumber=10
DataSetArray=($(cat InputSamples.txt)) # array of the input datadets
echo ${DataSetArray[$PROCESS / $SplitingNumber]}
DataSetName=${DataSetArray[$PROCESS / $SplitingNumber]}
rootNumber=$(($PROCESS % $SplitingNumber))
#DataSetName=DataSetName_.replace("/store/user/abdollah/Moriond17/","")

########### complie the Skimmer
make


########### loop over all root file in a dataset directory
xrdfs root://cmseos.fnal.gov ls "/eos/uscms/store/user/abdollah/Moriond17/"$DataSetName | grep $rootNumber.root | while read FullDataSetName

############  Here is where the Skimmer is running     ############
do
 file=`echo $FullDataSetName`
 ShortName=${file##*Moriond17}  # This removes all the string before Moriond17 (including Moriond17)
 echo "Here is the short Name   ------>" $ShortName
 ./Skimmer  $ShortName
done
############  Here is where the Skimmer ends          ############



IFS="/"
set $DataSetName
OutName=$4$6$rootNumber".root"  # this makes the 4th and 6th pieces of the
echo $OutName
hadd -f $OutName "skimed_"*.root




##########  remove the unneccesat files
#rm skim*root  Skimmer.cc  Skimmer.h  Makefile  InputSamples.txt
echo "Done execution ..."
xrdcp $OutName  root://cmseos.fnal.gov//store/user/abdollah/Condor/$OutName



END_TIME=`/bin/date`
echo "finished at ${END_TIME}"

exit $exitcode
echo "DONE!"

