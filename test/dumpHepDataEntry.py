from plotBinnedSigStr import varDic, xlabels
import hepdata_lib as hdl
import argparse
import pickle as pkl
import json
import oyaml as yaml


# nameDic = {
    
# }

xlabels  = {
    'Pt': (r'p_{\mathrm{T}}^{\gamma\gamma}', r'GeV'),
    'Njets2p5': r'n_{\mathrm{jets}}',
    'CosThetaS': r'\left|\cos\theta^{\ast}\right|',
    'AbsRapidity': r'\left|y^{\gamma\gamma}\right|',
    'Jet2p5AbsRapidity0': r'\left|y^{j_{1}}\right|',
    'Jet2p5Pt0': (r'p_{\mathrm{T}}^{j_{1}}', r'GeV'),
    'TauCJets4p7': (r'\tau_{\mathrm{C}}^{j}', r'GeV'),
    'AbsDeltaRapidityGgJet2p50': r'\left|\Delta y_{\gamma\gamma,j_{1}}\right|',
    'AbsDeltaPhiGgJet2p50': r'\left|\Delta\phi_{\gamma\gamma,j_{1}}\right|',
    'Jet4p7Pt1': (r'p_{\mathrm{T}}^{j_{2}}', r'GeV'),
    'Jet4p7AbsRapidity1': r'\left|y^{j_{2}}\right|',
    'AbsDeltaPhiGgJjJets4p7': r'|\Delta\phi_{\gamma\gamma,j_{1}j_{2}}|',
    'AbsDeltaPhiJ1J2Jets4p7': r'|\Delta\phi_{j_{1},j_{2}}|',
    'AbsZeppenfeldEtaJets4p7': r'|\bar{\eta}_{j_{1},j_{2}}-\eta_{\gamma\gamma}|',
    'MjjJets4p7': (r'm_{\mathrm{jj}}', 'GeV'),
    'AbsDeltaEtaJ1J2Jets4p7': r'|\Delta\eta_{j_{1},j_{2}}|',
    'Nleptons': r'n_{\mathrm{leptons}}',
    'NBjets2p5': r'n_{\mathrm{b-jets}}',
    'PtMiss': (r'p_{\mathrm{T}}^{\mathrm{miss}}', r'GeV'),
    'AbsPhiS': r'\left|\phi_{\eta}^{\ast}\right|',
}

hepDataLabels = {
    'Pt': 'PT(2GAMMA)',
    'Njets2p5': 'N(JET)',
    'CosThetaS': 'COS(THETA*)',
    'AbsRapidity': 'ABS(YRAP(2GAMMA))',
    'Jet2p5AbsRapidity0': 'YRAP(JET)',
    'Jet2p5Pt0': 'PT(JET)',
    'TauCJets4p7': 'TAUC(JET)',
    'AbsDeltaRapidityGgJet2p50': 'ABS(DELTAYRAP(2GAMMA,JET))',
    'AbsDeltaPhiGgJet2p50': 'ABS(DELTAPHI(2GAMMA,JET))',
    'Jet4p7Pt1': 'PT(JET2)',
    'Jet4p7AbsRapidity1': 'ABS(YRAP(JET2))',
    'AbsDeltaPhiGgJjJets4p7': 'ABS(DELTAPHI(2GAMMA,2JET))',
    'AbsDeltaPhiJ1J2Jets4p7': 'ABS(DELTAPHI(JET1,JET2))',
    'AbsZeppenfeldEtaJets4p7': 'ABS(DELTAETA(2GAMMA,2JET))',
    'MjjJets4p7': 'M(2JET)',
    'AbsDeltaEtaJ1J2Jets4p7': 'ABS(DELTAETA(JET1,JET2))',
    'Nleptons': 'N(LEPTON)',
    'NBjets2p5': 'N(BJET)',
    'PtMiss': 'MET',
    'AbsPhiS': 'ABS(PHI*(2GAMMA))'
}

locDic = {
    'Pt' : 'Fig. 10 upper left',
    'Njets2p5' : 'Fig. 10 upper right',
    'AbsRapidityFine' : 'Fig. 10 lower left',
    'CosThetaS' : 'Fig. 10 lower right',
    'Jet2p5Pt0': 'Fig. 11 lower left',
    'Jet2p5AbsRapidity0': 'Fig. 11 lower right',
    'AbsDeltaRapidityGgJet2p50': 'Fig. 12 upper left',
    'AbsDeltaPhiGgJet2p50': 'Fig. 12 upper right',
    'Jet4p7Pt1': 'Fig. 12 lower left',
    'Jet4p7AbsRapidity1': 'Fig. 12 lower right',
    'AbsDeltaPhiGgJjJets4p7': 'Fig. 13 upper left',
    'AbsDeltaPhiJ1J2Jets4p7': 'Fig. 13 upper right',
    'AbsZeppenfeldEtaJets4p7': 'Fig. 13 lower left',
    'MjjJets4p7NewBins': 'Fig. 13 lower right',
    'AbsDeltaEtaJ1J2Jets4p7': 'Fig. 14 upper left',
    'Nleptons': 'Fig. 14 lower left',
    'NBjets2p5': 'Fig. 14 lower right',
    'PtMiss': 'Fig. 14 upper right',
    'Jet4p7Pt1VBFlike': 'Fig. 15 upper left',
    'AbsDeltaPhiGgJjJets4p7VBFlike': 'Fig. 15 upper right',
    'AbsDeltaPhiJ1J2Jets4p7VBFlike': 'Fig. 15 lower left',
    'PtVBFlikeMergeB': 'Fig. 15 lower right',
    'TauCJets4p7': 'Fig. 11 upper right',
    'AbsPhiS': 'Fig. 11 upper left',
    'PtTauCJets4p7MergeB': 'Fig. 17',
    'PtNjets2p5': 'Fig. 16'
}

logList = [
    'Pt',
    'Njets2p5',
    'AbsPhiS',
    'TauCJets4p7',
    'AbsDeltaPhiGgJet2p50',
    'AbsDeltaPhiGgJjJets4p7',
    'MjjJets4p7NewBins',
    'PtMiss',
    'Nleptons',
    'NBjets2p5',
    'PtTauCJets4p7MergeB',
    'PtNjets2p5'
]

overflowList = [
    'Pt',
    'AbsPhiS',
    'TauCJets4p7',
    'Jet2p5Pt0',
    'AbsDeltaRapidityGgJet2p50',
    'Jet4p7Pt1',
    'AbsZeppenfeldEtaJets4p7',
    'MjjJets4p7NewBins',
    'AbsDeltaEtaJ1J2Jets4p7',
    'PtMiss',
    'Jet4p7Pt1VBFlike',
    'PtVBFlikeMergeB',
    'PtNjets2p5',
    'PtTauCJets4p7MergeB'
]

fidSel = {
    'fid' : {
        'PT(GAMMA1)/M(2GAMMA)' : ('>0.33', ''),
        'PT(GAMMA2)/M(2GAMMA)' : ('>0.25', ''),
        'ISOGEN(GAMMA)' : ('<10', 'GeV'),
        'ABS(ETA(GAMMA))' : ('<2.5', '')
     },
    'jet2p5' : {
        'PT(JET)' : ('>30', 'GeV'),
        'ABS(ETA(JET))' : ('<2.5', '')
    },
    'jet4p7' : {
        'PT(JET)' : ('>30', 'GeV'),
        'ABS(ETA(JET))' : ('<4.7', '')
    },
    'vbf' : {
        'N(JET)' : ('>=2', ''),
        'M(2JET)' : ('>200', 'GeV'),
        'DELTAETA(JET1, JET2)' : ('>3.5', '')
    }
}


def makeDiffTable(ext, binns, obsed, obsedErrsStat, obsedErrsSyst, exp, expErrs, inDir, binNr=None, topVar=None):

    if len(binns) <= 2:
        return None
    
    variable = varDic[ext]
    if isinstance(xlabels[variable], tuple):
        var, dvUnits = variable, xlabels[variable][1]
    else:
        var, dvUnits = variable, None
    hdlDVar = hdl.Variable(
        hepDataLabels[variable],
        is_independent=True,
        is_binned=True,
        units = dvUnits
    )
    
    hdlDVar.values = [(binns[j], binns[j+1]) for j in range(len(binns)-1)]
    hdlXs = hdl.Variable(
        'Observed',
        is_independent = False,
        is_binned = False,
        units = 'fb' if dvUnits is None else 'fb/{}'.format(dvUnits)
    )
    hdlXs.values = obsed
    hdlXsUncStat = hdl.Uncertainty("Stat.", is_symmetric=False)
    hdlXsUncStat.values = [(-obsedErrsStat[0][k], obsedErrsStat[1][k]) for k in range(len(obsedErrsStat[0]))]
    hdlXs.add_uncertainty(hdlXsUncStat)
    hdlXsUncStat = hdl.Uncertainty("Syst.", is_symmetric=False)
    hdlXsUncStat.values = [(-obsedErrsSyst[0][k], obsedErrsSyst[1][k]) for k in range(len(obsedErrsSyst[0]))]
    hdlXs.add_uncertainty(hdlXsUncStat)
    hdlXs.add_qualifier('SQRTS', '13', 'TeV')
    hdlXs.add_qualifier('MH','125.38','GeV')
    for key, item in fidSel['fid'].items():
        hdlXs.add_qualifier(key, item[0], item[1])
    if any([stt in ext for stt in ['Jet2p5', 'jets2p5']]):
        for key, item in fidSel['jet2p5'].items():
            hdlXs.add_qualifier(key, item[0], item[1])
    if any([stt in ext for stt in ['Jets4p7', 'Jet4p7']]):
        for key, item in fidSel['jet4p7'].items():
            hdlXs.add_qualifier(key, item[0], item[1])
    if 'VBFlike' in ext:
        for key, item in fidSel['vbf'].items():
            hdlXs.add_qualifier(key, item[0], item[1])

    hdlExp = hdl.Variable(
        'Expected',
        is_independent = False,
        is_binned = False,
        units = 'fb' if dvUnits is None else 'fb/{}'.format(dvUnits)
    )
    hdlExp.values = exp
    hdlExpUnc = hdl.Uncertainty("Theory", is_symmetric=False)
    hdlExpUnc.values = [(-expErrs[0][k], expErrs[1][k]) for k in range(len(expErrs[0]))]
    hdlExp.add_uncertainty(hdlExpUnc)
    hdlExp.add_qualifier('SQRTS', '13', 'TeV')
    hdlExp.add_qualifier('MH','125.38','GeV')
    for key, item in fidSel['fid'].items():
        hdlExp.add_qualifier(key, item[0], item[1])
    if any([stt in ext for stt in ['Jet2p5', 'jets2p5']]):
        for key, item in fidSel['jet2p5'].items():
            hdlExp.add_qualifier(key, item[0], item[1])
    if any([stt in ext for stt in ['Jets4p7', 'Jet4p7']]):
        for key, item in fidSel['jet4p7'].items():
            hdlExp.add_qualifier(key, item[0], item[1])
    if 'VBFlike' in ext:
        for key, item in fidSel['vbf'].items():
            hdlExp.add_qualifier(key, item[0], item[1])

    vName = hepDataLabels[variable]
    tabName = r'Diff xs {}'.format(vName)
    if topVar is not None:
        topVarHepD = hepDataLabels[topVar]
        tabName += r' vs. {}'.format(topVarHepD)
    if binNr is not None:
        tabName += ' bin {}'.format(binNr)
    if 'VBFlike' in ext:
        tabName += ' (VBF PS)'
    table = hdl.Table(tabName)
    vNameDesc = xlabels[variable][0] if isinstance(xlabels[variable], tuple) else xlabels[variable]
    table.description = 'Differential fiducial higgs to diphoton cross section with respect to ${}$'.format(vNameDesc)
    if topVar is not None:
        topVarTex = xlabels[topVar][0] if isinstance(xlabels[topVar], tuple) else xlabels[topVar]
        table.description += r' vs. ${}$'.format(topVarTex)
    if 'VBFlike' in ext:
        table.description += ' in the VBF enriched PS'
    if ext in overflowList:
        table.description += '. The last bin in the differential observable extends to infinity and the measured fiducial cross section in this bin is devided by the given bin width'
    if topVar is not None:
        table.keywords['observables'] = ['D2SIG/{}/{}'.format(hepDataLabels[variable], topVarHepD)]
    else:
        table.keywords['observables'] = ['DSIG/{}'.format(hepDataLabels[variable])]
    table.keywords['reactions'] = ['P P --> HIGGS (--> GAMMA GAMMA) X']
    table.keywords['cmenergies'] = ['13000']
    if ext in logList:
        table.add_image('{0}/outdir_differential_{1}/sigStr_{1}_ylog_wSyst_Data_commonPredHX.pdf'.format(inDir, ext))
    else:
        table.add_image('{0}/outdir_differential_{1}/sigStr_{1}_wSyst_Data_commonPredHX.pdf'.format(inDir, ext))
    table.location = 'Data from {}'.format(locDic[ext])
    table.add_variable(hdlDVar)
    table.add_variable(hdlXs)
    table.add_variable(hdlExp)

    return table


def makeCorrMatrix(ext, cMap, inDir, translate, topVar=None):

    variable = varDic[ext]
    with open('/eos/home-t/threiten/Analysis/Differentials/FinalFitsInDir/m125_{}_16/proc_cat_names_gen{}.txt'.format(ext, variable)) as f:
        lines = f.readlines()
        f.close()
    procs = lines[0].split(',')
    if len(procs) <= 2:
        return None
    pars = ['r{}'.format(i) for i in range(len(procs))]
    parNames = [translate[pars[i]] for i in range(len(procs))]
    for i in range(len(parNames)):
        while '#' in parNames[i]:
            parNames[i] = '${}$'.format(parNames[i].replace('#', '\\'))
            print(parNames[i])
    print(parNames)
        
    vName = hepDataLabels[variable]
    tabName = r'Correlation {}'.format(vName)
    if topVar is not None:
        topVarHepD = hepDataLabels[topVar]
        tabName += r' vs. {}'.format(topVarHepD)
    if 'VBFlike' in ext:
        tabName += ' (VBF PS)'
    if len(tabName) > 64:
        tabName = tabName.replace('Correlation', 'Corr')
    print(tabName)
    table = hdl.Table(tabName)
    
    varX = hdl.Variable(
        'Fiducial cross section (x)',
        is_independent=True,
        is_binned=False
    )
    varY = hdl.Variable(
        'Fiducial cross section (y)',
        is_independent=True,
        is_binned=False
    )
    corr = hdl.Variable(
        'Observed correlation',
        is_independent=False,
        is_binned=False
    )
    corrs = []
    xVals = []
    yVals = []
    for i,iPar in enumerate(pars):
        for j,jPar in enumerate(pars):
            xVals.append(parNames[i])
            yVals.append(parNames[j])
            if (iPar, jPar) in cMap.keys():
                corrs.append(cMap[(iPar, jPar)])
            elif (jPar, iPar) in cMap.keys():
                corrs.append(cMap[(jPar, iPar)])
    varX.values = xVals
    varY.values = yVals
    corr.values = corrs
    for key, item in fidSel['fid'].items():
        corr.add_qualifier(key, item[0], item[1])
    if any([stt in ext for stt in ['Jet2p5', 'jets2p5']]):
        for key, item in fidSel['jet2p5'].items():
            corr.add_qualifier(key, item[0], item[1])
    if any([stt in ext for stt in ['Jets4p7', 'Jet4p7']]):
        for key, item in fidSel['jet4p7'].items():
            corr.add_qualifier(key, item[0], item[1])
    if 'VBFlike' in ext:
        for key, item in fidSel['vbf'].items():
            corr.add_qualifier(key, item[0], item[1])
    corr.add_qualifier('SQRTS', '13', 'TeV')
    corr.add_qualifier('MH','125.38','GeV') 

    table.keywords['observables'] = ['CORR']
    table.keywords['reactions'] = ['P P --> HIGGS (--> GAMMA GAMMA) X']
    table.keywords['cmenergies'] = ['13000']
    if ext == 'Pt':
        locStr = 'Data from Fig. 9 upper'
    elif ext == 'Njets2p5':
        locStr = 'Data from Fig. 9 lower'
    else:
        locStr = 'Data from additional material'

    vNameDesc = xlabels[variable][0] if isinstance(xlabels[variable], tuple) else xlabels[variable]
    table.description = 'Correlation between the measured fiducial cross sections in the different bins of ${}$'.format(vNameDesc)
    if topVar is not None:
        topVarTex = xlabels[topVar][0] if isinstance(xlabels[topVar], tuple) else xlabels[topVar]
        table.description += ' and ${}$'.format(topVarTex)
    table.location = locStr
    table.add_variable(varX)
    table.add_variable(varY)
    table.add_variable(corr)
    table.add_image('{0}/outdir_differential_{1}/corrMatrix_{1}_robustHesse.pdf'.format(inDir, ext))

    return table

def makeFidTable(ress, errsStat, errsSyst, expRes, expErrs, inDir):

    conf = yaml.load(open('/afs/cern.ch/work/t/threiten/Hgg/Differentials/newCombine/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/test/fidXSPlot.yaml', 'r'), Loader=yaml.FullLoader)
    res = []
    resErrStat = []
    resErrSyst = []
    exp = []
    expErr = []
    labels = []
    for key in conf.keys():
        if key == 'global':
            continue
        res.append(ress[key])
        resErrStat.append(errsStat[key])
        resErrSyst.append(errsSyst[key])
        exp.append(expRes[key])
        expErr.append(expErrs[key])
        labels.append(conf[key]['label'])

    labVar = hdl.Variable(
        'Phase space region',
        is_independent=True,
        is_binned=False
    )
    labVar.values = labels
    xsObs = hdl.Variable(
        'Observed fiducial cross section',
        is_independent=False,
        is_binned=False,
        units='fb'
    )
    xsObs.values = res
    obsUncStat = hdl.Uncertainty('Stat.', is_symmetric=False)
    obsUncStat.values = [(resErrStat[k][0], resErrStat[k][1]) for k in range(len(resErrStat))]
    xsObs.add_uncertainty(obsUncStat)
    obsUncSyst = hdl.Uncertainty('Syst.', is_symmetric=False)
    obsUncSyst.values = [(resErrSyst[k][0], resErrSyst[k][1]) for k in range(len(resErrSyst))]
    xsObs.add_uncertainty(obsUncSyst)
    for key, item in fidSel['fid'].items():
        xsObs.add_qualifier(key, item[0], item[1])
    xsObs.add_qualifier('SQRTS', '13', 'TeV')
    xsObs.add_qualifier('MH','125.38','GeV')
    xsExp = hdl.Variable(
        'Expected fiducial cross section',
        is_independent=False,
        is_binned=False,
        units='fb'
    )
    xsExp.values = exp
    expUnc = hdl.Uncertainty('Theory', is_symmetric=False)
    expUnc.values = [(-expErr[k][0], expErr[k][1]) for k in range(len(expErr))]
    xsExp.add_uncertainty(expUnc)
    for key, item in fidSel['fid'].items():
        xsExp.add_qualifier(key, item[0], item[1])
    xsExp.add_qualifier('SQRTS', '13', 'TeV')
    xsExp.add_qualifier('MH','125.38','GeV')
    
    table = hdl.Table('Fiducial cross sections summary')
    table.description = 'Cross sections measured in different regions of the fiducial phase space'
    table.location = 'Data from Fig. 8'
    table.add_variable(labVar)
    table.add_variable(xsObs)
    table.add_variable(xsExp)
    table.add_image('{0}/FidXSData/plotFidXS_flipped.pdf'.format(inDir))

    return table

def main(options):

    sub = hdl.Submission()
    sub.add_link('CDS', 'http://cds.cern.ch/record/2803740')
    sub.add_record_id(2803740, 'CDS')

    for ext in varDic.keys():
        if ext == 'PtInclusiveVBF':
            continue
        binns = pkl.load(open('{0}/outdir_differential_{1}/bins_{1}.pkl'.format(options.inDir, ext), 'rb'))
        obsed = pkl.load(open('{0}/outdir_differential_{1}/observedResult_{1}.pkl'.format(options.inDir, ext), 'rb'))
        obsedErrsStat = pkl.load(open('{0}/outdir_differential_{1}/observedUncStat_{1}.pkl'.format(options.inDir, ext), 'rb'))
        obsedErrsSyst = pkl.load(open('{0}/outdir_differential_{1}/observedUncSyst_{1}.pkl'.format(options.inDir, ext), 'rb'))
        exp = pkl.load(open('{0}/outdir_differential_{1}/expectedResult_{1}.pkl'.format(options.inDir, ext), 'rb'))
        expErrs = pkl.load(open('{0}/outdir_differential_{1}/expectedResultUnc_{1}.pkl'.format(options.inDir, ext), 'rb'))
        print(ext)
        if ext == 'PtNjets2p5':
            topVar = 'Njets2p5'
        elif ext == 'PtTauCJets4p7MergeB':
            topVar = 'TauCJets4p7'
        else:
            topVar = None
        print(topVar)
        for i in range(len(binns)):
            binNr = i if len(binns)>1 else None
            table = makeDiffTable(ext, binns[i], obsed[i], obsedErrsStat[i], obsedErrsSyst[i], exp[i], expErrs[i], inDir=options.inDir, binNr=binNr, topVar=topVar)
            if table is not None:
                sub.add_table(table)

        cMap = pkl.load(open('{0}/outdir_differential_{1}/corrMap_{1}.pkl'.format(options.inDir, ext), 'rb'))
        translate = json.load(open('/afs/cern.ch/work/t/threiten/Hgg/Differentials/CMSSW_10_2_13/src/flashggFinalFit/Plots/translate_{}.json'.format(ext), 'r'))
        cTable = makeCorrMatrix(ext, cMap, options.inDir, translate, topVar=topVar)
        if cTable is not None:
            sub.add_table(cTable)

    fidRes = pkl.load(open('{}/FidXSData/result_fidXS.pkl'.format(options.inDir), 'rb'))
    fidResUncStat = pkl.load(open('{}/FidXSData/resultUncStat_fidXS.pkl'.format(options.inDir), 'rb'))
    fidResUncSyst = pkl.load(open('{}/FidXSData/resultUncSyst_fidXS.pkl'.format(options.inDir), 'rb'))
    fidExp = pkl.load(open('{}/FidXSData/resultExp_fidXS.pkl'.format(options.inDir), 'rb'))
    fidExpUnc = pkl.load(open('{}/FidXSData/resultExpUnc_fidXS.pkl'.format(options.inDir), 'rb'))
    table = makeFidTable(fidRes, fidResUncStat, fidResUncSyst, fidExp, fidExpUnc, options.inDir)
    sub.add_table(table)

    sub.create_files(options.hepDataDir, validate=True, remove_old=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    requiredArgs = parser.add_argument_group()
    requiredArgs.add_argument(
        '--inDir', '-d', action='store', type=str, required=True)
    requiredArgs.add_argument(
        '--hepDataDir', '-p', action='store', type=str, required=True)
    options = parser.parse_args()
    main(options)

