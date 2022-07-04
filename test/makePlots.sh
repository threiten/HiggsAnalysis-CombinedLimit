#!/bin/bash

extension=$1
baseDir=$2
mergeBs=$3

declare -A varDic
varDic[Pt]="Pt"
varDic[Njets2p5]="Njets2p5"
varDic[CosThetaS]="CosThetaS"
varDic[AbsRapidity]="AbsRapidity"
varDic[AbsRapidityFine]="AbsRapidity"
varDic[Jet2p5Pt0]="Jet2p5Pt0"
varDic[Jet2p5AbsRapidity0]="Jet2p5AbsRapidity0"
varDic[AbsDeltaRapidityGgJet2p50]="AbsDeltaRapidityGgJet2p50"
varDic[AbsDeltaPhiGgJet2p50]="AbsDeltaPhiGgJet2p50"
varDic[Jet4p7Pt1]="Jet4p7Pt1"
varDic[Jet4p7Pt1FixUF]="Jet4p7Pt1"
varDic[Jet4p7AbsRapidity1]="Jet4p7AbsRapidity1"
varDic[AbsDeltaPhiGgJjJets4p7]="AbsDeltaPhiGgJjJets4p7"
varDic[AbsDeltaPhiJ1J2Jets4p7]="AbsDeltaPhiJ1J2Jets4p7"
varDic[AbsZeppenfeldEtaJets4p7]="AbsZeppenfeldEtaJets4p7"
varDic[MjjJets4p7]="MjjJets4p7"
varDic[MjjJets4p7NewBins]="MjjJets4p7"
varDic[AbsDeltaEtaJ1J2Jets4p7]="AbsDeltaEtaJ1J2Jets4p7"
varDic[Nleptons]="Nleptons"
varDic[NBjets2p5]="NBjets2p5"
varDic[PtMiss]="PtMiss"
varDic[Jet4p7Pt1VBFlike]="Jet4p7Pt1"
varDic[AbsDeltaPhiGgJjJets4p7VBFlike]="AbsDeltaPhiGgJjJets4p7"
varDic[AbsDeltaPhiJ1J2Jets4p7VBFlike]="AbsDeltaPhiJ1J2Jets4p7"
varDic[PtVBFlike]="Pt"
varDic[PtVBFlikeMergeB]="Pt"
varDic[Pt0Jets]="Pt"
varDic[Pt1Jets]="Pt"
varDic[Pt1PJets]="Pt"
varDic[TauCJet2p50]="TauCJet2p50"
varDic[TauCJets4p7]="TauCJets4p7"
varDic[AbsPhiS]="AbsPhiS"
varDic[PtTauC0]="Pt"
varDic[PtTauC1]="Pt"
varDic[PtTauC2]="Pt"
varDic[PtTauC3]="Pt"
varDic[PtInclusive]="Pt"
varDic[PtInclusiveVBF]="Pt"
varDic[PtInclusive1L1B]="Pt"
varDic[PtInclusive1LHPtM]="Pt"
varDic[PtInclusive1LLPtM]="Pt"
varDic[PtTauCJet2p50]="Pt"
varDic[PtTauCJets4p7]="Pt"
varDic[PtTauCJets4p7MergeB]="Pt"
varDic[PtNjets2p5]="Pt"

variable=${varDic[$extension]}
cd $baseDir/outdir_differential_$extension

n=0
declare -a procStr=$(head -n1 /eos/home-t/threiten/Analysis/Differentials/FinalFitsInDir/m125_${extension}_16/proc_cat_names_gen${variable}.txt | tail -1)
for proc in ${procStr//\,/\ }
do
    let "n=$n+1"
done
let "n=$n-1"

for i in $(seq 0 $n)
do
    paramStrNoOne=$paramStrNoOne,r${i}
done

paramStrNoOne=${paramStrNoOne:1}

for bRepl in ${mergeBs//\,/\ }
do
    paramStrNoOne=$(echo $paramStrNoOne | sed -e "s/"${bRepl}"//" )
    paramStrNoOne=$(echo $paramStrNoOne | sed -e "s/\,\,/\,/")
done

if [[ ${paramStrNoOne: -1} == "," ]]; then
    paramStrNoOne=${paramStrNoOne::-1}
fi

echo $paramStrNoOne

for param in ${paramStrNoOne//\,/\ }
do
    hadd -f higgsCombineAsimovPostFitScanFit_${param}.MultiDimFit.mH125.38.root higgsCombineAsimovPostFitScanFit_${param}.POINTS.*.MultiDimFit.mH125.38.root
done

for param in ${paramStrNoOne//\,/\ }
do
    hadd -f higgsCombineAsimovPostFitScanStat_${param}.MultiDimFit.mH125.38.root higgsCombineAsimovPostFitScanStat_${param}.POINTS.*.MultiDimFit.mH125.38.root
done

for param in ${paramStrNoOne//\,/\ }
do
    plot1DScan.py --POI ${param} -o ${param}Scan --main-label="Expected stat+syst unc." higgsCombineAsimovPostFitScanFit_${param}.MultiDimFit.mH125.38.root --other higgsCombineAsimovPostFitScanStat_${param}.MultiDimFit.mH125.38.root:"Expected stat unc.":2
done
