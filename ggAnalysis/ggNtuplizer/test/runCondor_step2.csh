#!/bin/tcsh                                                                                                                                     

foreach x (`seq 0 422`)

/bin/rm -f condor_W_${x}
/bin/rm -f submit_W_${x}.csh

cat > submit_W_${x}.csh << +EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
cd /uscms_data/d3/snandan/CMSSW_8_0_26_patch1/src/ggAnalysis/ggNtuplizer/test/
setenv SCRAM_ARCH slc6_amd64_gcc530
eval `scram runtime -csh`
cmsenv
cd \${_CONDOR_SCRATCH_DIR}
cp /uscms_data/d3/snandan/CMSSW_8_0_26_patch1/src/ggAnalysis/ggNtuplizer/test/run_mc80X_W_${x}.py .
cp /uscms_data/d3/snandan/CMSSW_8_0_26_patch1/src/ggAnalysis/ggNtuplizer/test/Summer16_23Sep2016V4_MC.db .
cmsRun run_mc80X_W_${x}.py
rm run_mc80X_W_${x}.py
rm Summer16_23Sep2016V4_MC.db
echo "Done"
xrdcp -f ggtree_mc_${x}.root root://cmseos.fnal.gov//store/user/snandan/W
rm ggtree_mc_${x}.root
+EOF

chmod 775 submit_W_${x}.csh

cat > condor_W_${x} << +EOF
Universe       = vanilla
Executable     = /uscms_data/d3/snandan/CMSSW_8_0_26_patch1/src/ggAnalysis/ggNtuplizer/test/submit_W_${x}.csh
use_x509userproxy = true
Requirements   =  OpSys == "LINUX" && (Arch =="INTEL" || Arch =="x86_64")
request_memory = 10000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
output  = rungg_W_${x}.out
error   = rungg_W_${x}.error       
Log     = rungg_W_${x}.log
Queue 1

+EOF

condor_submit condor_W_${x}

end
