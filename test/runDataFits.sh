#!/bin/bash

ext=$1
baseDir=$2
mergeBs=$3

declare -a procMap=$(head -n1 ${baseDir}/outdir_differential_$ext/proc_cat_names_gen*.txt | tail -1)

n=0;
for proc in ${procMap//\,/\ }
do
    let "n=$n+1"
done
let "n=$n-1"

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

cd $baseDir
cd outdir_differential_$ext

combine Datacard_13TeV_differential_${ext}.root -M MultiDimFit -m 125.38 --saveWorkspace -n DataBestFit --setParameters $paramStr --redefineSignalPOIs $paramStrNoOne --X-rtd MINIMIZER_freezeDisassociatedParams --cminDefaultMinimizerStrategy 1 # --setRobustFitStrategy=2 --robustFit=1 #  --robustHesse=1

for param in ${paramStrNoOne//\,/\ }
do
    combineTool.py higgsCombineDataBestFit.MultiDimFit.mH125.38.root -M MultiDimFit -m 125.38 --snapshotName "MultiDimFit" -v -1 -P ${param} --floatOtherPOIs 1 -n DataScan_${param} --redefineSignalPOIs $paramStrNoOne --X-rtd MINIMIZER_freezeDisassociatedParams --cminDefaultMinimizerStrategy 1 --algo grid --points 30 --squareDistPoiStep --split-points 1 --job-mode condor --task-name DataScan_${param} --sub-opts='+JobFlavour = "microcentury"' # --setRobustFitStrategy=2 --robustFit=1
    combineTool.py higgsCombineDataBestFit.MultiDimFit.mH125.38.root -M MultiDimFit -m 125.38 --snapshotName "MultiDimFit" -v -1 -P ${param} --floatOtherPOIs 1 -n DataScanStat_${param} --redefineSignalPOIs $paramStrNoOne --X-rtd MINIMIZER_freezeDisassociatedParams --cminDefaultMinimizerStrategy 1 --algo grid --points 30 --squareDistPoiStep --freezeParameters allConstrainedNuisances --split-points 1 --job-mode condor --task-name DataScanStat_${param} --sub-opts='+JobFlavour = "microcentury"' # --setRobustFitStrategy=2 --robustFit=1
done

