#!/bin/bash

ext=$1
baseDir=$2
mergeBs=$3

# paramStr="rOut=1"

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

combineTool.py higgsCombineDataBestFit.MultiDimFit.mH125.38.root --snapshotName "MultiDimFit" -M MultiDimFit -m 125.38 --saveWorkspace -n DataSMCompat --redefineSignalPOIs ${paramStrNoOne} --X-rtd MINIMIZER_freezeDisassociatedParams --cminDefaultMinimizerStrategy 1 --algo fixed --fixedPointPOIs ${paramStr} --saveNLL --skipInitialFit #--job-mode condor --task-name DataSMCompat --sub-opts='+JobFlavour = "testmatch"'
