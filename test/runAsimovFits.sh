#!/bin/bash

ext=$1
mergeBs=$2

declare -a procMap=$(head -n1 ./outdir_differential_$ext/proc_cat_names_gen*.txt | tail -1)

n=0;
for proc in ${procMap//\,/\ }
do
    let "n=$n+1"
done
let "n=$n-1"

# paramStr="rOut=1"
for i in $(seq 0 $n)
do
    paramStr=$paramStr,r${i}=1
done

paramStr=${paramStr:1}

for i in $(seq 0 $n)
do
    paramStrNoOne=$paramStrNoOne,r${i}
done

paramStrNoOne=${paramStrNoOne:1}

for bRepl in ${mergeBs//\,/\ }
do
    paramStr=$(echo $paramStr | sed -e "s/"${bRepl}"=1//" )
    paramStr=$(echo $paramStr | sed -e "s/\,\,/\,/")
    paramStrNoOne=$(echo $paramStrNoOne | sed -e "s/"${bRepl}"//" )
    paramStrNoOne=$(echo $paramStrNoOne | sed -e "s/\,\,/\,/")
done

if [[ ${paramStr: -1} == "," ]]; then
    paramStr=${paramStr::-1}
fi

if [[ ${paramStrNoOne: -1} == "," ]]; then
    paramStrNoOne=${paramStrNoOne::-1}
fi

echo $paramStr
echo $paramStrNoOne

cd outdir_differential_$ext

combine -M GenerateOnly -m 125.38 --setParameters $paramStr -n AsimovPreFit --saveToys -t -1 Datacard_13TeV_differential_${ext}.root

combine Datacard_13TeV_differential_${ext}.root -t -1 --toysFile higgsCombineAsimovPreFit.GenerateOnly.mH125.38.123456.root -M MultiDimFit -m 125.38 --saveWorkspace -n AsimovBestFit --X-rtd MINIMIZER_freezeDisassociatedParams --cminDefaultMinimizerStrategy 0 --setParameters $paramStr --redefineSignalPOIs ${paramStrNoOne}

combine -M GenerateOnly -m 125.38 --setParameters $paramStr -n AsimovPostFit --saveToys --saveWorkspace --snapshotName MultiDimFit -t -1  higgsCombineAsimovBestFit.MultiDimFit.mH125.38.root

for param in ${paramStrNoOne//\,/\ }
do
    combineTool.py higgsCombineAsimovPostFit.GenerateOnly.mH125.38.123456.root -M MultiDimFit -m 125.38 --snapshotName "MultiDimFit" -t -1 -v -1 -P ${param} --floatOtherPOIs 1 -n AsimovPostFitScanFit_${param} --X-rtd MINIMIZER_freezeDisassociatedParams --cminDefaultMinimizerStrategy 0 --algo grid --points 30 --squareDistPoiStep --setParameters $paramStr --redefineSignalPOIs ${paramStrNoOne} --toysFile higgsCombineAsimovPostFit.GenerateOnly.mH125.38.123456.root --split-points 1 --job-mode condor --task-name AsimovScan_${param} --sub-opts='+JobFlavour = "microcentury"'
    
    combineTool.py higgsCombineAsimovPostFit.GenerateOnly.mH125.38.123456.root -M MultiDimFit -m 125.38 --snapshotName "MultiDimFit" -t -1 -v -1 -P ${param} --floatOtherPOIs 1 -n AsimovPostFitScanStat_${param}  --X-rtd MINIMIZER_freezeDisassociatedParams --cminDefaultMinimizerStrategy 0 --algo grid --points 30 --squareDistPoiStep --setParameters $paramStr --redefineSignalPOIs ${paramStrNoOne} --toysFile higgsCombineAsimovPostFit.GenerateOnly.mH125.38.123456.root --freezeParameters allConstrainedNuisances --split-points 1 --job-mode condor --task-name AsimovScanStat_${param} --sub-opts='+JobFlavour = "microcentury"'
done

