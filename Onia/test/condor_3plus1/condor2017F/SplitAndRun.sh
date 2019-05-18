#!/bin/bash

cmsswDir="/uscms_data/d3/lpcbphy/zhenhu/Upsilon/fourmuon/CMSSW_9_4_0/src/"
inputFiles=""

myeospath1="/store/group/l1upgrades/Run2017/fourmuon/MuOnia/BPHSkim-v4-Run2017F-17Nov2017-v1/180314_061742/0000/"
myeospath2="/store/group/l1upgrades/Run2017/fourmuon/ZeroBias/BPHSkim-v6-Run2017F-17Nov2017-v1/180321_073042/0000/"

j=0
for files in `ls /eos/uscms${myeospath1} | grep root`
do 
	inputFiles=$files;
	echo $inputFiles;
	jobNb=${j};
	let j=${j}+1;
	#anaHeader="myntuple_${jobNb}.h";
	#anaCScript="myntuple_${jobNb}.C";
	jobCScript="runMixingRootupler_${jobNb}.py";
	scriptName="Run_${jobNb}.csh";
	condorScriptName="runOnCondor_${jobNb}";
	#cat myntuple.h | sed "s/NUMBER/${jobNb}/g" > ${anaHeader};
	#cat myntuple.C | sed "s/NUMBER/${jobNb}/g" > ${anaCScript};
	cat runMixingRootupler.py | sed "s:INPUTPATH1:${myeospath1}:" | sed "s-INPUTFILE-${inputFiles}-" | sed "s/NUMBER/${jobNb}/g" | sed "s:MIXINPUTPATH:${myeospath2}:" > ${jobCScript} ;
	cat Run.csh | sed "s-FILENAME-${jobCScript}-" > ${scriptName};
	chmod +x ${scriptName}
	cat runOnCondor | sed "s/SCRIPT/${scriptName}/" | sed "s/JOBC/${jobCScript}/" > ${condorScriptName}
	condor_submit ${condorScriptName}
done
