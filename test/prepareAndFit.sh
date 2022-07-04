#!/bin/bash

ext=$1
baseDir=$2

declare -A ranges
ranges[Pt]="\[1,-1,3\]"
ranges[Njets2p5]="\[1,-1,3\]"
ranges[Njets2p5_testNames]="\[1,0,2\]"
ranges[CosThetaS]="\[1,0,2\]"
ranges[AbsRapidity]="\[1,0,2\]"
ranges[AbsRapidityFine]="\[1,-1,3\]"
ranges[Jet2p5Pt0]="\[1,-1,3\]"
ranges[Jet2p5AbsRapidity0]="\[1,-1,3\]"
ranges[AbsDeltaRapidityGgJet2p50]="\[1,-1,3\]"
ranges[AbsDeltaPhiGgJet2p50]="\[1,-1,3\]"
ranges[Jet4p7Pt1]="\[1,-2,4\]"
ranges[Jet4p7AbsRapidity1]="\[1,-1,3\]"
ranges[AbsDeltaPhiGgJjJets4p7]="\[1,-2,4\]"
ranges[AbsDeltaPhiJ1J2Jets4p7]="\[1,-1,3\]"
ranges[AbsZeppenfeldEtaJets4p7]="\[1,-1,3\]"
ranges[MjjJets4p7]="\[1,-1,3\]"
ranges[MjjJets4p7NewBins]="\[1,-1,3\]"
ranges[AbsDeltaEtaJ1J2Jets4p7]="\[1,-1,3\]"
ranges[Nleptons]="\[1,-1,3\]"
ranges[NBjets2p5]="\[1,-2,4\]"
ranges[PtMiss]="\[1,-1,3\]"
ranges[Jet4p7Pt1VBFlike]="\[1,-3,5\]"
ranges[AbsDeltaPhiGgJjJets4p7VBFlike]="\[1,-3,5\]"
ranges[AbsDeltaPhiJ1J2Jets4p7VBFlike]="\[1,-3,5\]"
ranges[PtVBFlike]="\[1,-3,5\]"
ranges[Pt0Jets]="\[1,0,2\]"
ranges[Pt1Jets]="\[1,-1,3\]"
ranges[Pt1PJets]="\[1,-1,3\]"
ranges[TauCJet2p50]="\[1,-1,3\]"
ranges[AbsPhiS]="\[1,0,2\]"
ranges[PtTauC0]="\[1,-3,5\]"
ranges[PtTauC1]="\[1,-3,5\]"
ranges[PtTauC2]="\[1,-3,5\]"
ranges[PtTauC3]="\[1,-3,5\]"
ranges[PtInclusive]="\[1,0,2\]"
ranges[PtInclusive1L1B]="\[1,-1,3\]"
ranges[PtInclusive1LHPtM]="\[1,-1,3\]"
ranges[PtInclusive1LLPtM]="\[1,-1,3\]"
ranges[PtTauCJet2p50]="\[1,-1,3\]"
ranges[PtTauCJets4p7]="\[1,-1,3\]"
ranges[PtNjets2p5]="\[1,-1,3\]"

declare -A varDic
varDic[Pt]="Pt"
varDic[Njets2p5]="Njets2p5"
varDic[Njets2p5_testNames]="Njets2p5"
varDic[CosThetaS]="CosThetaS"
varDic[AbsRapidity]="AbsRapidity"
varDic[AbsRapidityFine]="AbsRapidity"
varDic[Jet2p5Pt0]="Jet2p5Pt0"
varDic[Jet2p5AbsRapidity0]="Jet2p5AbsRapidity0"
varDic[AbsDeltaRapidityGgJet2p50]="AbsDeltaRapidityGgJet2p50"
varDic[AbsDeltaPhiGgJet2p50]="AbsDeltaPhiGgJet2p50"
varDic[Jet4p7Pt1]="Jet4p7Pt1"
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
varDic[Pt0Jets]="Pt"
varDic[Pt1Jets]="Pt"
varDic[Pt1PJets]="Pt"
varDic[TauCJet2p50]="TauCJet2p50"
varDic[AbsPhiS]="AbsPhiS"
varDic[PtTauC0]="Pt"
varDic[PtTauC1]="Pt"
varDic[PtTauC2]="Pt"
varDic[PtTauC3]="Pt"
varDic[PtInclusive]="Pt"
varDic[PtInclusive1L1B]="Pt"
varDic[PtInclusive1LHPtM]="Pt"
varDic[PtInclusive1LLPtM]="Pt"
varDic[PtTauCJet2p50]="Pt"
varDic[PtTauCJets4p7]="Pt"
varDic[PtNjets2p5]="Pt"

var=${varDic[$ext]}
echo $var

cd $baseDir

mkdir outdir_differential_$ext
cp /afs/cern.ch/work/t/threiten/Hgg/Differentials/CMSSW_10_2_13/src/flashggFinalFit/Datacard/Datacard_13TeV_differential_$ext.txt ./outdir_differential_$ext/

python /afs/cern.ch/work/t/threiten/Hgg/Differentials/CMSSW_10_2_13/src/flashggFinalFit/Datacard/cleanDatacard.py -d ./outdir_differential_${ext}/Datacard_13TeV_differential_${ext}.txt -o ./outdir_differential_${ext}/Datacard_13TeV_differential_${ext}_cleaned.txt --removeDoubleSided --factor=1.5
mv ./outdir_differential_${ext}/Datacard_13TeV_differential_${ext}.txt ./outdir_differential_${ext}/Datacard_13TeV_differential_${ext}_orig.txt
mv ./outdir_differential_${ext}/Datacard_13TeV_differential_${ext}_cleaned.txt ./outdir_differential_${ext}/Datacard_13TeV_differential_${ext}.txt

declare -a procMap=$(head -n1 /eos/home-t/threiten/Analysis/Differentials/FinalFitsInDir/m125_${ext}_16/proc_cat_names_gen${var}.txt | tail -1)
declare -a procMapO=$(head -n1 /eos/home-t/threiten/Analysis/Differentials/FinalFitsInDir/m125_${ext}_16/proc_cat_names_gen${var}.txt | tail -1),OutsideAcceptance

for yr in 16 17 18
do
    for proc in ${procMapO//\,/\ }
    do
        expr="s/CMS-HGG_sigfit_differential_"$ext"_"$proc"_hgg.*"$yr".root/CMS-HGG_sigfit_"$proc"_"$yr".root/"
	echo $expr
	sed -i -e $expr outdir_differential_$ext/Datacard_13TeV_differential_$ext.txt
    done
done

for yr in 16 17 18
do
    sed -i "s/CMS-HGG_differential_"$ext"_13TeV_multipdf.*"$yr".root/CMS-HGG_multipdf_differential_"$ext"_"$yr".root/" outdir_differential_$ext/Datacard_13TeV_differential_$ext.txt
done

for yr in 16 17 18
do
    find /afs/cern.ch/work/t/threiten/Hgg/Differentials/CMSSW_10_2_13/src/flashggFinalFit/Signal/outdir_differential_${ext}_${yr}/ \( -name CMS-HGG_sigfit_smH\*.root -o -name CMS-HGG_sigfit_OutsideAcceptance\*.root \) \! -name CMS-HGG_sigfit_differential\*.root -exec cp '{}' ./outdir_differential_$ext/ \;
    exp="s/.root/_$yr.root/"
    for proc in ${procMapO//\,/\ }
    do
	mv outdir_differential_$ext/CMS-HGG_sigfit_$proc.root $(echo outdir_differential_$ext/CMS-HGG_sigfit_$proc.root | sed -e $exp)
    done
done

for yr in 16 17 18
do
    cp /afs/cern.ch/work/t/threiten/Hgg/Differentials/CMSSW_10_2_13/src/flashggFinalFit/Background/CMS-HGG_multipdf_differential_${ext}_${yr}.root outdir_differential_$ext/
done

i=0;
for proc in ${procMap//\,/\ }
do
    if [[ $proc == *"m1000"* ]]; then
	mapStr=$(echo $mapStr $(echo $proc | sed -e "s/smH/--PO map=\.\*\/smH/" | sed -e "s/\$/:r$i\[1,0,2\]'/"))
    else
	mapStr=$(echo $mapStr $(echo $proc | sed -e "s/smH/--PO map=\.\*\/smH/" | sed -e "s/\$/:r$i${ranges[$ext]}'/"))
    fi
    let "i=$i+1"
done

# mapStr="--PO map=.*/OutsideAcceptance:rOut[1,-1,3] "$mapStr

echo $mapStr
echo $i

text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose -m 125 --PO 'higgsMassRange=123,127' $(echo $mapStr) -o outdir_differential_$ext/Datacard_13TeV_differential_${ext}.root outdir_differential_$ext/Datacard_13TeV_differential_${ext}.txt

/afs/cern.ch/work/t/threiten/Hgg/Differentials/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/test/runAsimovFits.sh $ext $i
