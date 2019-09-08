#!/bin/bash
#cmsswDir="/store/user/cdozen/FourMuon_Analysis/NTuples/2017/ZeroBias/"
cmsswDir="/uscms_data/d3/cdozen/CMSSW_8_0_29/src/"
#cmsswDir="/uscms_data/d3/lpcbphy/zhenhu/Upsilon/fourmuon/CMSSW_9_4_0/src/"
inputFiles=""
inputFiles2=""

myeospath1="/store/user/cdozen/FourMuon_Analysis/MuOnia/2016_v2/MuOnia/BPHSkim--Run2016B-07Aug17_ver2-v1/190620_113034/0000"
myeospath2="/store/user/cdozen/FourMuon_Analysis/MiniBias/2016_v2/ZeroBias/BPHSkim--Run2016B-07Aug17_ver2-v1/190620_113508/0000"

j=0
files2=(`ls /eos/uscms${myeospath1} | grep root`)
for files in `ls /eos/uscms${myeospath2} | grep root`
do 
	inputFiles2=${files2[j]};
	inputFiles=$files;
	echo $inputFiles;
	echo $inputFiles2;
	jobNb=${j};
	let j=${j}+1;
	#anaHeader="myntuple_${jobNb}.h";
	#anaCScript="myntuple_${jobNb}.C";
	jobCScript="runMixingRootupler_${jobNb}.py";
	scriptName="Run_${jobNb}.csh";
	condorScriptName="runOnCondor_${jobNb}";
	#cat myntuple.h | sed "s/NUMBER/${jobNb}/g" > ${anaHeader};
	#cat myntuple.C | sed "s/NUMBER/${jobNb}/g" > ${anaCScript};
	cat runMixingRootupler.py | sed "s:INPUTPATH1:${myeospath1}:" | sed "s-INPUTFILE-${inputFiles2}-" | sed "s/NUMBER/${jobNb}/g" | sed "s:MIXINPUTPATH:${myeospath2}:" | sed "s-MIXFILEINPUT-${inputFiles}-" > ${jobCScript} ;
	cat Run.csh | sed "s-FILENAME-${jobCScript}-" > ${scriptName};
	chmod +x ${scriptName}
	cat runOnCondor | sed "s/SCRIPT/${scriptName}/" | sed "s/JOBC/${jobCScript}/" > ${condorScriptName}
	condor_submit ${condorScriptName}
done
