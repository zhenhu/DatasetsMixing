# DatasetsMixing
Eventing mixing with two different datasets

To submit the CMSSW jobs through condor system at **FNAL LPC**, please follow _Example B_ in the instruction below:
https://uscms.org/uscms_at_work/computing/setup/batch_systems.shtml

**Step 1**: compile and test your code locally

**Step 2**: tar your CMSSW working area and copy it to your personal EOS area
```
tar -zcvf CMSSW940.tgz CMSSW_9_4_0 --exclude='CMSSW_9_4_0/src/DatasetsMixing/Onia/test/condor*' --exclude='*.root' --exclude='CMSSW_9_4_0/src/DatasetsMixing/.git' --exclude='CMSSW_9_4_0/src/.git' --exclude='CMSSW_9_4_0/tmp'
cp CMSSW940.tgz /eos/uscms/store/user/zhenhu/ 
```
**Step 3**: submit jobs
```
 cd CMSSW_9_4_0/src/DatasetsMixing/Onia/test/condor_3plus1/
 cmsenv
 voms-proxy-init -voms cms
 ./submit.sh 
```
