import root_pandas
import scipy.interpolate as itr
import scipy.optimize as opt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
import hepdata_lib as hdl
import pickle as pkl
import oyaml as yaml
import argparse
import copy
from scipy.stats import chi2


xlabels  = {
    'Pt': (r'p_{\mathrm{T}}^{\PGg\PGg}', r'GeV'),
    'Njets2p5': r'n_{\mathrm{jets}}',
    'CosThetaS': r'\left|\cos\theta^{\ast}\right|',
    'AbsRapidity': r'\left|y^{\PGg\PGg}\right|',
    'Jet2p5AbsRapidity0': r'\left|y^{\mathrm{j}_{1}}\right|',
    'Jet2p5Pt0': (r'p_{\mathrm{T}}^{\mathrm{j}_{1}}', r'GeV'),
    'TauCJets4p7': (r'\tau_{\mathrm{C}}^{\mathrm{j}}', r'GeV'),
    'AbsDeltaRapidityGgJet2p50': r'\left|\Delta y_{\PGg\PGg,j_{1}}\right|',
    'AbsDeltaPhiGgJet2p50': r'\left|\Delta\phi_{\PGg\PGg,j_{1}}\right|',
    'Jet4p7Pt1': (r'p_{\mathrm{T}}^{\mathrm{j}_{2}}', r'GeV'),
    'Jet4p7AbsRapidity1': r'\left|y^{\mathrm{j}_{2}}\right|',
    'AbsDeltaPhiGgJjJets4p7': r'\left|\Delta\phi_{\PGg\PGg,\mathrm{j}_{1}\mathrm{j}_{2}}\right|',
    'AbsDeltaPhiJ1J2Jets4p7': r'\left|\Delta\phi_{\mathrm{j}_{1},\mathrm{j}_{2}}\right|',
    'AbsZeppenfeldEtaJets4p7': r'|\bar{\eta}_{\mathrm{j}_{1},\mathrm{j}_{2}}-\eta_{\PGg\PGg}|',
    'MjjJets4p7': (r'm_{\mathrm{jj}}', 'GeV'),
    'AbsDeltaEtaJ1J2Jets4p7': r'\left|\Delta\eta_{\mathrm{j}_{1},\mathrm{j}_{2}}\right|',
    'Nleptons': r'n_{\mathrm{leptons}}',
    'NBjets2p5': r'n_{\mathrm{\PQb jets}}',
    'PtMiss': (r'p_{\mathrm{T}}^{\mathrm{miss}}', r'GeV'),
    'AbsPhiS': r'\left|\phi_{\eta}^{\ast}\right|',
}

cutStr = {
    'Pt': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\end{flushleft}',
    'Njets2p5': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}}\right|<2.5$\\$p_{\mathrm{T}}^{\mathrm{j}}>30\,\textnormal{GeV}$\end{flushleft}',
    'CosThetaS': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\end{flushleft}',
    'AbsRapidity': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\end{flushleft}',
    'Jet2p5AbsRapidity0': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{1}}\right|<2.5$\\$p_{\mathrm{T}}^{\mathrm{j}_{1}}>30\,\textnormal{GeV}$\end{flushleft}',
    'Jet2p5Pt0': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{1}}\right|<2.5$\\$p_{\mathrm{T}}^{\mathrm{j}_{1}}>30\,\textnormal{GeV}$\end{flushleft}',
    'TauCJets4p7': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}}\right|<4.7$\end{flushleft}',
    'AbsDeltaRapidityGgJet2p50': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{1}}\right|<2.5$\\$p_{\mathrm{T}}^{\mathrm{j}_{1}}>30\,\textnormal{GeV}$\end{flushleft}',
    'AbsDeltaPhiGgJet2p50': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{1}}\right|<2.5$\\$p_{\mathrm{T}}^{\mathrm{j}_{1}}>30\,\textnormal{GeV}$\end{flushleft}',
    'Jet4p7Pt1': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{2}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{2}}>30\,\textnormal{GeV}$\end{flushleft}',
    'Jet4p7AbsRapidity1': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{2}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{2}}>30\,\textnormal{GeV}$\end{flushleft}',
    'AbsDeltaPhiGgJjJets4p7': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\end{flushleft}',
    'AbsDeltaPhiJ1J2Jets4p7': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\end{flushleft}',
    'AbsZeppenfeldEtaJets4p7': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\end{flushleft}',
    'MjjJets4p7': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\end{flushleft}',
    'AbsDeltaEtaJ1J2Jets4p7': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\end{flushleft}',
    'Nleptons': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\end{flushleft}',
    'NBjets2p5': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}}\right|<2.5$\\$p_{\mathrm{T}}^{\mathrm{j}}>30\,\textnormal{GeV}$\end{flushleft}',
    'PtMiss': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\end{flushleft}',
    'AbsPhiS': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\end{flushleft}',
    'PtVBFlike': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\\$n_{\mathrm{jets}}\geq2$\\$\Delta\eta_{j1,j2}>3.5$\\$m_{jj}>200\,\textnormal{GeV}$\end{flushleft}',
    'PtVBFlikeMergeB': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\\$n_{\mathrm{jets}}\geq2$\\$\Delta\eta_{j1,j2}>3.5$\\$m_{jj}>200\,\textnormal{GeV}$\end{flushleft}',
    'Jet4p7Pt1VBFlike': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\\$n_{\mathrm{jets}}\geq2$\\$\Delta\eta_{j1,j2}>3.5$\\$m_{jj}>200\,\textnormal{GeV}$\end{flushleft}',
    'AbsDeltaPhiGgJjJets4p7VBFlike': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\\$n_{\mathrm{jets}}\geq2$\\$\Delta\eta_{j1,j2}>3.5$\\$m_{jj}>200\,\textnormal{GeV}$\end{flushleft}',
    'AbsDeltaPhiJ1J2Jets4p7VBFlike': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{i}}\right|<4.7$\\$p_{\mathrm{T}}^{\mathrm{j}_{i}}>30\,\textnormal{GeV}$\\$n_{\mathrm{jets}}\geq2$\\$\Delta\eta_{j1,j2}>3.5$\\$m_{jj}>200\,\textnormal{GeV}$\end{flushleft}',
    'PtTauCJet2p50': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}_{1}}\right|<2.5$\\$p_{\mathrm{T}}^{\mathrm{j}_{1}}>30\,\textnormal{GeV}$\end{flushleft}',
    'PtTauCJets4p7': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}}\right|<4.7$\end{flushleft}',
    'PtTauCJets4p7MergeB': r'\begin{flushleft}$\boldsymbol{\PH\to\PGg\PGg}$\\$\left|\eta^{\mathrm{j}}\right|<4.7$\end{flushleft}',
}

varDic = {
    'Pt': 'Pt',
    'Njets2p5': 'Njets2p5',
    'CosThetaS': 'CosThetaS',
    'AbsRapidityFine': 'AbsRapidity',
    'Jet2p5Pt0': 'Jet2p5Pt0',
    'Jet2p5AbsRapidity0': 'Jet2p5AbsRapidity0',
    'AbsDeltaRapidityGgJet2p50': 'AbsDeltaRapidityGgJet2p50',
    'AbsDeltaPhiGgJet2p50': 'AbsDeltaPhiGgJet2p50',
    'Jet4p7Pt1': 'Jet4p7Pt1',
    'Jet4p7AbsRapidity1': 'Jet4p7AbsRapidity1',
    'AbsDeltaPhiGgJjJets4p7': 'AbsDeltaPhiGgJjJets4p7',
    'AbsDeltaPhiJ1J2Jets4p7': 'AbsDeltaPhiJ1J2Jets4p7',
    'AbsZeppenfeldEtaJets4p7': 'AbsZeppenfeldEtaJets4p7',
    'MjjJets4p7NewBins': 'MjjJets4p7',
    'AbsDeltaEtaJ1J2Jets4p7': 'AbsDeltaEtaJ1J2Jets4p7',
    'Nleptons': 'Nleptons',
    'NBjets2p5': 'NBjets2p5',
    'PtMiss': 'PtMiss',
    'Jet4p7Pt1VBFlike': 'Jet4p7Pt1',
    'AbsDeltaPhiGgJjJets4p7VBFlike': 'AbsDeltaPhiGgJjJets4p7',
    'AbsDeltaPhiJ1J2Jets4p7VBFlike': 'AbsDeltaPhiJ1J2Jets4p7',
    'PtVBFlikeMergeB': 'Pt',
    'TauCJets4p7': 'TauCJets4p7',
    'AbsPhiS': 'AbsPhiS',
    'PtInclusive': 'Pt',
    'PtInclusiveVBF': 'Pt',
    'PtInclusive1L1B': 'Pt',
    'PtInclusive1LHPtM': 'Pt',
    'PtInclusive1LLPtM': 'Pt',
    'PtTauCJets4p7MergeB': 'Pt',
    'PtNjets2p5': 'Pt'
}

lastBins = {
    'AbsRapidity': 2.5,
    'Jet2p5AbsRapidity0': 2.5,
    'AbsDeltaPhiGgJet2p50': np.pi,
    'Jet4p7AbsRapidity1': 4.7,
    'AbsDeltaPhiJ1J2Jets4p7': np.pi,
    'AbsDeltaPhiGgJjJets4p7': np.pi
}

fHDic = {
    'PtVBFlike' : 11,
    'PtVBFlikeMergeB' : 12,
    'Jet4p7Pt1VBFlike': 11,
    'AbsDeltaPhiGgJjJets4p7VBFlike': 11,
    'AbsDeltaPhiJ1J2Jets4p7VBFlike': 11.6,
    'PtNjets2p5' : 10.7,
    'PtTauCJets4p7' : 9,
    'PtTauCJets4p7MergeB' : 9,
    'AbsZeppenfeldEtaJets4p7' : 11
}

xTicksDic = {
    'Njets2p5' : ([0, 1, 2, 3, 4], [r'$0$', r'$1$', r'$2$', r'$3$', r'$\geq4$']),
    'NBjets2p5' : ([0, 1, 2], [r'$0$', r'$1$', r'$\geq2$']),
    'Nleptons' : ([0, 1, 2], [r'$0$', r'$1$', r'$\geq2$'])
}

ePredLDict = {
    'nominal' : r'\textsc{MadGraph}\xspace{}5\_a\textsc{mc@nlo}, \textsc{NNLOPS}',
    'noNNLOPS' : r'\textsc{MadGraph}\xspace{}5\_a\textsc{mc@nlo}',
    'noGGH' : r'\textsc{MadGraph}\xspace{}5\_a\textsc{mc@nlo} VBF+VH+ttH',
    'noGluGluH' : r'\textsc{MadGraph}\xspace{}5\_a\textsc{mc@nlo} VBF+VH+ttH',
    'HX' : r'\textsc{MadGraph}\xspace{}5\_a\textsc{mc@nlo} VBF+VH+ttH',
    'powheg' : r'{\textsc{powheg}}\xspace',
    'powhegVBFdipoleRecoil' : r'{\textsc{powheg}}\xspace (VBF dipole Recoil)',
    'ggH' : r'\textsc{MadGraph}\xspace{}5\_a\textsc{mc@nlo}, NNLOPS\xspace ggH',
    'noNNLOPSggH' : r'\textsc{MadGraph}\xspace{}5\_a\textsc{mc@nlo}\xspace ggH',
    'powhegggH' : r'{\textsc{powheg}}\xspace ggH',
    'powhegVBFdipoleRecoilggH' : r'{\textsc{powheg}}\xspace (VBF dipole Recoil)\xspace ggH'
}

commonPredLDict = {
    'noNNLOPS' : r'HX',
    'noGGH' : r'HX',
    'noGluGluH' : r'HX',
    'HX' : r'HX'
}

offset = {
    'nominal' : -1./2.,
    'ggH' : -1./2.,
    'noNNLOPS' : 1./3.,
    'noGGH' : 1./3.,
    'noGluGluH' : 1./3.,
    'HX' : 1./3.,
    'powheg' : 2./3.,
    'powhegVBFdipoleRecoil' : -1./3.,
    'noNNLOPSggH' : 1./3.,
    'powhegggH' : 2./3.,
    'powhegVBFdipoleRecoilggH' : -1./3.
}


def find1Sig(sStr, deltaNLL, minim=None, retSpl=False):
    TdeltaNLL = 2.* deltaNLL
    cSpl = itr.interp1d(sStr,TdeltaNLL,'cubic')
    minimum = minim if minim is not None else sStr[0]
    try:
        sigM = opt.root_scalar(lambda x: cSpl(x) - 1., bracket=(sStr.min(), sStr[0]))
        neg = sigM.root - minimum
    except ValueError:
        print('WARNING: Lower limit not found. Using placeholder -3')
        neg = -3. - minimum
    try:
        sigP = opt.root_scalar(lambda x: cSpl(x) - 1., bracket=(sStr[0], sStr.max()))
        pos = sigP.root - minimum
    except ValueError:
        print('WARNING: Upper limit not found. Using placeholder 5')
        pos = 5. - minimum

    if retSpl:
        return neg, pos, cSpl
    
    return neg, pos


def getNLL(par, df, cutoff=100.):
    dNLLRaw = df['deltaNLL'].values
    sStrRaw = df[par].values
    deltaNLL = np.hstack((dNLLRaw[0],dNLLRaw[1:][np.logical_and(dNLLRaw[1:]!=dNLLRaw[0], np.abs(dNLLRaw[1:])<cutoff)]))
    sigStr = np.hstack((sStrRaw[0],sStrRaw[1:][np.logical_and(dNLLRaw[1:]!=dNLLRaw[0], np.abs(dNLLRaw[1:])<cutoff)]))
    
    return sigStr, deltaNLL


def getPos(frac, lims, log=False):
    if log:
        return 10 ** (frac * (np.log10(lims[1]) - np.log10(lims[0])) + np.log10(lims[0]))
    
    return frac * (lims[1]-lims[0]) + lims[0]


def getBins(var, ext, mergeBins=None):

    mergeBEx = False
    if mergeBins is None:
        mergeBins = []
        mergeBEx = True
    cfg = yaml.load(open('/afs/cern.ch/work/t/threiten/Hgg/Differentials/CMSSW_10_2_13/src/flashggFinalFit/splitConfig_{}_2016.yaml'.format(ext)), Loader=yaml.FullLoader)
    topLvl = None
    topVar = None
    try:
        ret = cfg['splits']['gen{}'.format(var)]
        retTmp = copy.deepcopy(ret)
        for bnInd in mergeBins:
            retTmp.remove(ret[bnInd])
        retUnrem = [np.array(copy.deepcopy(ret))]
        ret = retTmp
        ret = np.array(ret)
        ret = [ret]
    except KeyError:
        ret = cfg['splits']['gen']['gen{}'.format(var)]
        lensRet = np.cumsum(np.array([0] + [len(subret) - 1 for subret in ret]))
        retTmp = copy.deepcopy(ret)
        for bnInd in mergeBins:
            ind = np.searchsorted(lensRet, bnInd, side='right') - 1
            retTmp[ind].remove(ret[ind][bnInd-lensRet[ind]])
        retUnrem = copy.deepcopy(ret)
        for i in range(len(retUnrem)):
            retUnrem[i] = np.array(retUnrem[i])
        ret = retTmp
        for i in range(len(ret)):
            ret[i] = np.array(ret[i])
        topVar = list(cfg['splits']['gen'].keys())[0]
        topLvl = cfg['splits']['gen'][topVar]

    return ret, topLvl, topVar, retUnrem


def getPred(arr, mass, nBins, weights=None, interPRepl=None, massRepl=None):
    masses = [120., 125., 130.]
    if weights is None:
        weights = np.ones(arr.shape[1])
    splines = []
    if arr.shape[1] == 1:
        if interPRepl is None or massRepl is None:
            raise ValueError("If only one masspoint is given, interPRepl and massRepl must be provided!")
        for i in range(nBins):
            splines.append(itr.UnivariateSpline(masses, interPRepl[i,:], w=weights, k=2))

        return np.array([splines[i](mass)-interPRepl[i, masses.index(massRepl)]+arr[i,0] for i in range(nBins)])

    for i in range(nBins):
        splines.append(itr.UnivariateSpline(masses, arr[i,:], w=weights, k=2))

    return np.array([splines[i](mass) for i in range(nBins)])

cmsTextDic = {
    'wip' : r'\textbf{CMS} \textit{Work in Progress}',
    'prelim' : r'\textbf{CMS} \textit{Preliminary}',
    'final' : r'\textbf{CMS}',
    'sim' : r'\textbf{CMS} \textit{Simulation}',
}

def drawCMSLogo(ax, opt=None, fs=22):

    ax.text(0, 1, cmsTextDic[opt], fontsize=fs, transform=ax.transAxes, va='bottom')


def drawIntLumi(ax, intL=138, fs=22):

    if intL is None:
        ax.text(1, 1, r'13\ensuremath{{\,\text{{Te\hspace{{-.08em}}V}}}}\xspace', fontsize=fs, transform=ax.transAxes, ha='right', va='bottom')
    else:
        ax.text(1, 1, r'{}\mbox{{\ensuremath{{\,\text{{fb}}^{{-1}}}}}}\xspace (13\ensuremath{{\,\text{{Te\hspace{{-.08em}}V}}}}\xspace)'.format(intL), fontsize=fs, transform=ax.transAxes, ha='right', va='bottom')


def main(options):

    rcP = {'text.usetex': True,
           'font.family': 'sans-serif',
           'font.sans-serif': ['Helvetica'],
           'pdf.fonttype': 42,
           'axes.labelsize': 20,
           'font.size': 16,
           'pgf.rcfonts': True,
           'text.latex.preamble': r"\usepackage{bm, xspace, amsmath, heppennames2}"}

    plt.rcParams.update(rcP)
    predColors = list(cm.Set1.colors)

    variable = varDic[options.extension]
    binns, topLvl, topVar, binnsUnrem = getBins(variable, options.extension, mergeBins=options.mergeBins)
    underflowExists = False
    cappedOverflow = [False] * len(binns)
    if topLvl is not None:
        topLvl = np.array(topLvl)
        topVar = topVar.replace('gen', '')
        if not options.showUnderflow and topLvl[0] == -1000:
            topLvl  = np.delete(topLvl, 0)
            underflowExists = True

    if underflowExists and len(binns)>1:
        binns.pop(0)
        
    nBins = 0
    bcen = []
    xerrs = []
    binws = []
    for k in range(len(binns)):
        if not underflowExists and not options.showUnderflow and binns[k][0] == -1000.:
            binns[k] = np.delete(binns[k], 0)
            underflowExists = True
        elif options.showUnderflow and binns[k][0] == -1000.:
            binns[k][0] = 0.
        if binns[k][-1] > 2 * binns[k][-2]:
            binns[k][-1] = 2 * binns[k][-2] - binns[k][-3] if len(binns[k])>2 else binns[k][-1]
            cappedOverflow[k] = True if len(binns[k])>2 else False
        if variable in lastBins.keys():
            binns[k][-1] = lastBins[variable] 
        if options.cap is not None:
            binns[k][-1] = options.cap
        bc = (0.5*(binns[k][1:]+binns[k][:-1]))
        bcen.append(bc)
        xerrs.append(binns[k][1:] - bc)
        for m in range(len(binns[k]) - 1):
            binws.append(binns[k][m+1] - binns[k][m])
        nBins += len(binns[k]) - 1
    binws = np.array(binws)
    
    if options.mergeBins is not None:
        binwUnrem = []
        for l in range(len(binnsUnrem)):
            if underflowExists:
                binnsUnrem[l] = np.delete(binnsUnrem[l], 0)
            elif binnsUnrem[l][0] == -1000.:
                binnsUnrem[l][0] = 0.
            if binnsUnrem[l][-1] > 2*binnsUnrem[l][-2]:
                binnsUnrem[l][-1] = 2 * binnsUnrem[l][-2] - binnsUnrem[l][-3] if len(binnsUnrem[l])>2 else binnsUnrem[l][-1]
            if variable in lastBins.keys():
                binnsUnrem[l][-1] = lastBins[variable]
            if options.cap is not None:
                binnsUnrem[l][-1] = options.cap
            for m in range(len(binnsUnrem[l]) - 1):
                binwUnrem.append(binnsUnrem[l][m+1] - binnsUnrem[l][m])
        binwUnrem = np.array(binwUnrem)

    print(nBins)
    add = 1 if underflowExists else 0
    if options.mergeBins is not None:
        add += len(options.mergeBins)
    parsToCons = list(range(nBins+add))
    if underflowExists:
        parsToCons.remove(0)
    if options.mergeBins is not None:
        for bnInd in options.mergeBins:
            parsToCons.remove(bnInd)

    fname = 'higgsCombineAsimovPostFitScanFit' if options.filename is None else options.filename
    errs = np.zeros((2,nBins))
    cval = np.zeros((nBins))
    for n in range(nBins):
        par = 'r{}'.format(parsToCons[n])
        df = root_pandas.read_root('{}/{}_{}.MultiDimFit.mH125.38.root'.format(options.directory, fname, par), 'limit')
        sStr, deltaNLL = getNLL(par, df)
        cval[n] = sStr[0]
        errs[:, n] = np.array(find1Sig(sStr, deltaNLL))

    if options.showSyst or options.showStat:
        fnameStat = 'higgsCombineAsimovPostFitScanStat' if options.filenameStat is None else options.filenameStat
        errsStat = np.zeros((2,nBins))
        cavlStat = np.zeros((nBins))
        for n in range(nBins):
            par = 'r{}'.format(parsToCons[n])
            dfStat = root_pandas.read_root('{}/{}_{}.MultiDimFit.mH125.38.root'.format(options.directory, fnameStat, par), 'limit')
            sStrStat, deltaNLLStat = getNLL(par, dfStat)
            cavlStat[n] = sStrStat[0]
            errsStat[:, n] = np.array(find1Sig(sStrStat, deltaNLLStat, cval[n]))
        errsSyst = np.sqrt(np.abs(np.power(errs,2)-np.power(errsStat,2))) #np.maximum(np.zeros_like(errs), np.power(errs,2)-np.power(errsStat,2))
        print('stat unc: ', errsStat)
        print('syst unc: ', errsSyst)

    mass = df['mh'].values[0]
    cval = np.array(cval)
    errs = np.abs(errs)

    with open(options.predFile, 'rb') as f:
        preds = pkl.load(f)
        f.close()

    extraPreds = []
    predLabels = {}
    ePredKeys = []
    predLabels['nominal'] = ePredLDict['nominal']
    if options.extraPredFiles is not None:
        for ePredF in options.extraPredFiles:
            with open(ePredF, 'rb') as f:
                extraPreds.append(pkl.load(f))
                f.close()
            ePredKeys.append(ePredF.split("_")[-1].split(".")[0])
            predLabels[ePredKeys[-1]] = ePredLDict[ePredF.split("_")[-1].split(".")[0]]

    if underflowExists:
        preds = np.delete(preds, 0, axis=0)
        for i in range(len(extraPreds)):
            extraPreds[i] = np.delete(extraPreds[i], 0, axis=0)

    print('binws: ', binws)
    
    if options.mergeBins is not None:
        # print('extraPreds before: ', extraPreds)
        print('binwUnrem: ', binwUnrem)
        indRed = 1 if underflowExists else 0
        for bnInd in options.mergeBins:
            preds[bnInd - indRed - 1] = (binwUnrem[bnInd - indRed - 1]*preds[bnInd - indRed - 1] + binwUnrem[bnInd - indRed]*preds[bnInd - indRed])
            for i in range(len(extraPreds)):
                extraPreds[i][bnInd - indRed - 1] = (binwUnrem[bnInd - indRed - 1]*extraPreds[i][bnInd - indRed - 1] + binwUnrem[bnInd - indRed]*extraPreds[i][bnInd - indRed])
        for bnInd in options.mergeBins:
            preds[bnInd - indRed - 1] /= binws[bnInd - indRed - 1]
            preds = np.delete(preds, bnInd - indRed, axis=0)
            for i in range(len(extraPreds)):
                extraPreds[i][bnInd - indRed - 1] /= binws[bnInd - indRed - 1]
                extraPreds[i] = np.delete(extraPreds[i], bnInd - indRed, axis=0)
            indRed += 1
        # print('extraPreds after: ', extraPreds)
    extraPred = []
    for ePreds in extraPreds:
        if ePreds.shape[1] == 1:
            extraPred.append(getPred(ePreds, mass, nBins, weights=[1., 2.3, 1.], interPRepl=preds, massRepl=125.))
        else:
            extraPred.append(getPred(ePreds, mass, nBins, weights=[1., 2.3, 1.]))

    pred = getPred(preds, mass, nBins)

    if options.commonPredComp is not None:
        with open(options.commonPredComp, 'rb') as cpF:
            commonPreds = pkl.load(cpF)
            cpF.close()
        if underflowExists:
            commonPreds = np.delete(commonPreds, 0, axis=0)
        if options.mergeBins is not None:
            indRed = 1 if underflowExists else 0
            for bnInd in options.mergeBins:
                commonPreds[bnInd - indRed - 1] = (binwUnrem[bnInd - indRed - 1]*commonPreds[bnInd - indRed - 1] + binwUnrem[bnInd - indRed]*commonPreds[bnInd - indRed])
            for bnInd in options.mergeBins:
                commonPreds[bnInd - indRed - 1] /= binws[bnInd - indRed - 1]
                commonPreds = np.delete(commonPreds, bnInd - indRed, axis=0)
                indRed += 1
        if commonPreds.shape[1] == 1:
            commonPred = getPred(commonPreds, mass, nBins, weights=[1., 2.3, 1.], interPRepl=preds, massRepl=125.)
        else:
            commonPred = getPred(commonPreds, mass, nBins, weights=[1., 2.3, 1.])
        commonPredKey = options.commonPredComp.split("_")[-1].split(".")[0]
        pred = pred + commonPred
        if options.predFile.split("_")[-1].split(".")[0] in ePredLDict.keys():
            predLabels['nominal'] = ePredLDict[options.predFile.split("_")[-1].split(".")[0]] + r' + ' + commonPredLDict[commonPredKey]
        else:
            predLabels['nominal'] = predLabels['nominal'] + r' + ' + commonPredLDict[commonPredKey]
        for k in range(len(extraPred)):
            if ePredKeys[k] == commonPredKey:
                predLabels[ePredKeys[k]] = commonPredLDict[commonPredKey] + r' = ' + predLabels[ePredKeys[k]]
                continue
            extraPred[k] = extraPred[k] + commonPred
            predLabels[ePredKeys[k]] = predLabels[ePredKeys[k]] + r' + ' + commonPredLDict[commonPredKey]
            
    if options.showUnderflow and binns[0][0] == 0.:
        pred[0] = pred[0]*((binns[0][1] + 1000.)/(binns[0][1] - 0.))
        for ePred in extraPred:
            ePred[0] = ePred[0]*((binns[0][1] + 1000.)/(binns[0][1] - 0.))
    obsed = pred * cval
    obsedErrs = np.array([pred, pred]).reshape(2, nBins) * errs
    if options.showSyst:
        errsSystObsed = np.array([pred, pred]).reshape(2, nBins) * errsSyst
    if options.showStat or options.dumpResults:
        errsStat = np.abs(errsStat)
        errsStatObsed = np.array([pred, pred]).reshape(2, nBins) * errsStat

    if options.showTheoryUnc:
        errsTheoryObsed = {}
        errsTheory = {}
        ePredKeysHere = copy.deepcopy(ePredKeys)
        for thSkipped in options.thErrToSkip:
            errsTheoryObsed[thSkipped] = np.zeros((2, nBins))
            ePredKeysHere.remove(thSkipped)

        if len(ePredKeysHere)+1 != len(options.thErrFiles):
            raise ValueError("Number of theory error files must be equal to number of predictions minus number of labels supplied to --thErrToSkip")
        
        for k, thErrFl in enumerate(options.thErrFiles):
            predLbl = ePredKeysHere[k-1] if k>0 else 'nominal'
            with open(thErrFl, 'rb') as fl:
                errsTheoryObsed[predLbl] = pkl.load(fl)
                fl.close()
        if underflowExists:
            for key in errsTheoryObsed.keys():
                if errsTheoryObsed[key].shape[1] > nBins:
                    errsTheoryObsed[key] = np.delete(errsTheoryObsed[key], 0, axis=1)
        print('errsTheoryObsed before: ', errsTheoryObsed)
        if options.mergeBins is not None:
            for key in errsTheoryObsed.keys():
                indRed = 1 if underflowExists else 0
                if errsTheoryObsed[key].shape[1] > nBins:
                    for bnInd in options.mergeBins:
                        print('errsTheoryPred ind-1 before: ', errsTheoryObsed[key][:, bnInd - indRed - 1])
                        print('errsTheoryPred ind before: ', errsTheoryObsed[key][:, bnInd - indRed])
                        errsTheoryObsed[key][:, bnInd - indRed - 1] = np.sqrt(np.power(binwUnrem[bnInd - indRed - 1]*errsTheoryObsed[key][:, bnInd - indRed - 1], 2) + np.power(binwUnrem[bnInd - indRed]*errsTheoryObsed[key][:, bnInd - indRed], 2))
                        print('errsTheoryPred single after: ', errsTheoryObsed[key][:, bnInd - indRed - 1])
                    for bnInd in options.mergeBins:
                        errsTheoryObsed[key][:, bnInd - indRed - 1] /= binws[bnInd - indRed - 1]
                        errsTheoryObsed[key] = np.delete(errsTheoryObsed[key], bnInd - indRed, axis=1)
                        indRed += 1
        print('errsTheoryObsed: ', errsTheoryObsed)
        print('pred: ', pred)
        errsTheory['nominal'] = errsTheoryObsed['nominal']/np.array([pred, pred]).reshape(2, nBins)
        for l, pd in enumerate(extraPred):
            ePredLbl = ePredKeys[l]
            errsTheory[ePredLbl] = errsTheoryObsed[ePredLbl]/np.array([pd, pd]).reshape(2, nBins)

    lengths = [len(binns[0]) - 1]
    for j in range(1, len(binns)):
        lengths.append(lengths[j-1] + len(binns[j]) - 1)
    cval = np.array_split(cval, lengths)
    errs = np.array_split(errs, lengths, axis=1)
    pred = np.array_split(pred, lengths)
    for m in range(len(extraPred)):
        extraPred[m] = np.array_split(extraPred[m], lengths)
    obsed = np.array_split(obsed, lengths)
    obsedErrs = np.array_split(obsedErrs, lengths, axis=1)
    if options.showSyst:
        errsSyst = np.array_split(errsSyst, lengths, axis=1)
        errsSystObsed = np.array_split(errsSystObsed, lengths, axis=1)
    if options.showStat or options.dumpResults:
        errsStat = np.array_split(errsStat, lengths, axis=1)
        errsStatObsed = np.array_split(errsStatObsed, lengths, axis=1)
    if options.showTheoryUnc:
        for key in errsTheory:
            errsTheory[key] = np.array_split(errsTheory[key], lengths, axis=1)
            errsTheoryObsed[key] = np.array_split(errsTheoryObsed[key], lengths, axis=1)

    if options.showPvalue:
        fnameSMCompat = 'higgsCombineDataSMCompat' if options.fnameSMCompat is None else options.fnameSMCompat
        dfSMCompat = root_pandas.read_root('{}/{}.MultiDimFit.mH125.38.root'.format(options.directory, fnameSMCompat), 'limit')
        deltaNLLSMCompat = dfSMCompat.loc[1, 'deltaNLL']
        # print('deltaNLLSMCompat: ', deltaNLLSMCompat)
        chi2pdf = chi2(nBins+1 if underflowExists else nBins)
        pvalSMCompat = 1 - chi2pdf.cdf(2*deltaNLLSMCompat)
        # print('pvalSMCompat: ', pvalSMCompat)
            
    figWidth = len(binns)/2 if len(binns) > 1 else len(binns)
    figHeight = fHDic[options.extension] if options.extension in fHDic.keys() else 8
    topHeight = figHeight - 2
    print(topHeight)
    fig, axes = plt.subplots(2, len(binns), figsize=(figWidth*8, figHeight), sharex='col', sharey='row', gridspec_kw={'height_ratios': [topHeight, 2]}, squeeze=False)
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.05, wspace=0.)
    tops = axes[0, :]
    bottoms = axes[1, :]

    # print('bcen', bcen)
    # print('errs', errs)
    # print('pred', pred)
    # print('xerrs', xerrs)
    # print('obsed', obsed)
    # print('obsedErrs', obsedErrs)
    
    xlb = xlabels[variable]
    if isinstance(xlb, tuple):
        xlabel = r'${{{0}}} \textnormal{{({{{1}}})}}$'.format(xlb[0], xlb[1])
    else:
        xlabel = r'${{{}}}$'.format(xlb)

    res_label = r'Data, stat$\oplus$syst unc.' if 'Data' in fname else r'Expected Result, stat$\oplus$syst unc.'

    # top.hist(bcen, bins=bins, weights=pred, histtype='step', label=r'amc@NLO, NNLOPS', linestyle='solid', linewidth=2)
    hatchstr = '/'*10
    for i in range(len(binns)):
        tops[i].errorbar(bcen[i], pred[i], xerr=xerrs[i], linestyle='none', markersize=0, capsize=0, elinewidth=2, label=predLabels['nominal'], color=predColors[0], zorder=4)
        for p, ePred in enumerate(extraPred):
            tops[i].errorbar(bcen[i], ePred[i], xerr=xerrs[i], linestyle='none', markersize=0, capsize=0, elinewidth=2, label=predLabels[ePredKeys[p]], color=predColors[p+1])

        (_,caps,_) = tops[i].errorbar(bcen[i], obsed[i], yerr=obsedErrs[i], xerr=xerrs[i], marker='.', linestyle='none',markersize=8, capsize=0, linewidth=2, label=res_label, color='black', zorder=5)
        for cap in caps:
            cap.set_markeredgewidth(1)

        if options.showSyst:
            errorboxesTop = [matplotlib.patches.Rectangle((x-xe, y-ye[0]), 2*xe, ye.sum(), zorder=2) for x, y, xe, ye in zip(bcen[i], obsed[i], xerrs[i], errsSystObsed[i].T)]
            pcTop = matplotlib.collections.PatchCollection(errorboxesTop, facecolor='lightslategray', alpha=0.5, edgecolor='None', zorder=3)
            topErrColl = tops[i].add_collection(pcTop)
            
        if options.showStat:
            errorboxesTopStat = [matplotlib.patches.Rectangle((x-xe, y-ye[0]), 2*xe, ye.sum(), zorder=2) for x, y, xe, ye in zip(bcen[i], obsed[i], xerrs[i], errsStatObsed[i].T)]
            pcTopStat = matplotlib.collections.PatchCollection(errorboxesTopStat, facecolor='lightslategray', alpha=0.85, edgecolor='None', zorder=3)
            topErrCollStat = tops[i].add_collection(pcTopStat)

        if options.showTheoryUnc:
            widthThUnc = (1./6.)*xerrs[i]
            errorboxesTopTheory = {}
            pcTopTheoryUnc = {}
            topErrCollTheoryUnc = {}
            for key, item in errsTheoryObsed.items():
                predHere = pred if key == 'nominal' else extraPred[ePredKeys.index(key)]
                cInd = 0 if key == 'nominal' else ePredKeys.index(key)+1
                errorboxesTopTheory[key] = [matplotlib.patches.Rectangle((x-xe, y-ye[0]), 2*xe, ye.sum(), zorder=2) for x, y, xe, ye in zip(bcen[i]+offset[key]*xerrs[i], predHere[i], widthThUnc, item[i].T)]
                pcTopTheoryUnc[key] = matplotlib.collections.PatchCollection(errorboxesTopTheory[key], hatch=hatchstr, edgecolor=predColors[cInd], facecolor='None', linewidth=0, zorder=3)
                topErrCollTheoryUnc[key] = tops[i].add_collection(pcTopTheoryUnc[key])

        (_,caps,_) = bottoms[i].errorbar(bcen[i], cval[i], yerr=errs[i], xerr=xerrs[i], marker='.', linestyle='none',markersize=8, capsize=0, linewidth=2, label=r'Expected Result, stat$\oplus$syst unc.', color='black', zorder=5)
        for cap in caps:
            cap.set_markeredgewidth(1)
        (_,caps,_) = bottoms[i].errorbar(bcen[i], np.ones_like(pred[i]), xerr=xerrs[i], markersize=0, linestyle='none', capsize=0, elinewidth=2, label=predLabels['nominal'], color=predColors[0], zorder=4)
        for p, ePred in enumerate(extraPred):
            if ePredKeys[p] in options.predToSkipRatio:
                continue
            (_,caps,_) = bottoms[i].errorbar(bcen[i], np.divide(ePred[i], pred[i]), xerr=xerrs[i], markersize=0, linestyle='none', capsize=0, elinewidth=2, label=predLabels[ePredKeys[p]], color=predColors[p+1])

        if options.showSyst:
            errorboxesBottom = [matplotlib.patches.Rectangle((x-xe, y-ye[0]), 2*xe, ye.sum()) for x, y, xe, ye in zip(bcen[i], cval[i], xerrs[i], errsSyst[i].T)]
            pcBottom = matplotlib.collections.PatchCollection(errorboxesBottom, facecolor=tuple(topErrColl.get_facecolor()[0]), alpha=0.5, edgecolor='None', zorder=3)
            bottoms[i].add_collection(pcBottom)

        if options.showStat:
            errorboxesBottomStat = [matplotlib.patches.Rectangle((x-xe, y-ye[0]), 2*xe, ye.sum()) for x, y, xe, ye in zip(bcen[i], cval[i], xerrs[i], errsStat[i].T)]
            pcBottomStat = matplotlib.collections.PatchCollection(errorboxesBottomStat, facecolor=tuple(topErrCollStat.get_facecolor()[0]), alpha=0.85, edgecolor='None', zorder=3)
            bottoms[i].add_collection(pcBottomStat)

        if options.showTheoryUnc:
            errorboxesBottomTheoryUnc = {}
            pcBottomTheoryUnc = {}
            for key, item in errsTheory.items():
                if key in options.predToSkipRatio:
                    continue
                predHere = pred if key == 'nominal' else extraPred[ePredKeys.index(key)]
                errorboxesBottomTheoryUnc[key] = [matplotlib.patches.Rectangle((x-xe, y-ye[0]), 2*xe, ye.sum()) for x, y, xe, ye in zip(bcen[i]+offset[key]*xerrs[i], np.divide(predHere[i], pred[i]), widthThUnc, item[i].T)]
                pcBottomTheoryUnc[key] = matplotlib.collections.PatchCollection(errorboxesBottomTheoryUnc[key], edgecolor=tuple(topErrCollTheoryUnc[key].get_edgecolor()[0]), facecolor='None', hatch=hatchstr, linewidth=0, zorder=3)
                bottoms[i].add_collection(pcBottomTheoryUnc[key])
            
        # if options.yLim is None:
        #     bottoms[i].set_ylim(-0.5, 2.5)
        if options.yLim is not None:
            bottoms[i].set_ylim(*options.yLim)
        bottoms[i].grid(linestyle='-.', color='lightslategrey', alpha=0.5)
        tops[i].grid(linestyle='-.', color='lightslategrey', alpha=0.5)

        lowBound = binns[i][0]
        if options.xlog:
            bottoms[i].set_xscale('log')
            if lowBound == 0.:
                lowBound = 0.1 * binns[i][1]
            
        if options.cap is not None:
            bottoms[i].set_xlim(lowBound, options.cap)
        else:
            bottoms[i].set_xlim(lowBound, binns[i][-1])

        bottoms[i].set_xlabel(xlabel)
        if options.extension in xTicksDic.keys():
            bottoms[i].set_xticks(xTicksDic[options.extension][0])
            bottoms[i].set_xticklabels(xTicksDic[options.extension][1])
        
        if cappedOverflow[i] and options.extension not in ['Njets2p5', 'NBjets2p5', 'Nleptons']:
            yPos = max([obsed[i][-1]+obsedErrs[i][1][-1]] + [pred[i][-1]] + [extraPred[k][i][-1] for k in range(len(extraPred))])
            if isinstance(xlb, tuple):
                fig.text(bcen[i][-1]+0.5*xerrs[i][-1], 1.1*yPos, r'$\boldsymbol{{\sigma_{{\mathrm{{fid}}}}\left({{{0}}}>{{{1}}}\,\textnormal{{{2}}}\right)/{{{3:.1f}}}\,\textnormal{{{2}}}}}$'.format(xlb[0], binns[i][-2], xlb[1], binns[i][-1] - binns[i][-2]), fontsize=10, transform=tops[i].transData, va='bottom', rotation=90) #-xerrs[i][-1]
            else:
                fig.text(bcen[i][-1]+0.5*xerrs[i][-1], 1.1*yPos, r'$\boldsymbol{{\sigma_{{\mathrm{{fid}}}}\left({{{0}}}>{{{1}}}\right)/{{{2:.1f}}}}}$'.format(xlb, binns[i][-2], binns[i][-1] - binns[i][-2]), fontsize=10, transform=tops[i].transData, va='bottom', rotation=90 )

        lowYlim = tops[i].get_ylim()[0]
        print(lowYlim)
                
        if options.commonPredComp is not None:
            cPredboxesTop = [matplotlib.patches.Rectangle((x0, y0), x1-x0, y1-y0, zorder=2) for x0, y0, x1, y1 in zip(bcen[i]-xerrs[i], lowYlim*np.ones_like(bcen[i]), bcen[i]+xerrs[i], extraPred[ePredKeys.index(commonPredKey)][i])]
            print(predColors[ePredKeys.index(commonPredKey)+1]+tuple([0.66]))
            cpbTop = matplotlib.collections.PatchCollection(cPredboxesTop, edgecolor=predColors[ePredKeys.index(commonPredKey)+1]+tuple([0.66]), facecolor='None', linestyle='None', linewidth=0, zorder=0) #, hatch='xxx'
            topcpbColl = tops[i].add_collection(cpbTop)
            cPredLines = matplotlib.collections.LineCollection([[[x,y0],[x,y1]] for x, y0, y1 in zip(bcen[i][:-1]+xerrs[i][:-1], extraPred[ePredKeys.index(commonPredKey)][i][:-1], extraPred[ePredKeys.index(commonPredKey)][i][1:])], linewidth=2, color=predColors[ePredKeys.index(commonPredKey)+1], capstyle='round')
            cPredLColl = tops[i].add_collection(cPredLines)

        if options.dumpResults:
            pkl.dump(binns, open('{}/bins_{}.pkl'.format(options.outDir, options.extension), 'wb'))
            pkl.dump(obsed, open('{}/observedResult_{}.pkl'.format(options.outDir, options.extension), 'wb'))
            pkl.dump(errsStatObsed, open('{}/observedUncStat_{}.pkl'.format(options.outDir, options.extension), 'wb'))
            pkl.dump(pred, open('{}/expectedResult_{}.pkl'.format(options.outDir, options.extension), 'wb'))
            pkl.dump(errsTheory['nominal'], open('{}/expectedResultUnc_{}.pkl'.format(options.outDir, options.extension), 'wb'))
            pkl.dump(errsSystObsed, open('{}/observedUncSyst_{}.pkl'.format(options.outDir, options.extension), 'wb'))
            # if isinstance(xlabels[variable], tuple):
            #     var, dvUnits = variable, xlabels[variable][1]
            # else:
            #     var, dvUnits = variable, None
            # hdlDVar = hdl.Variable(
            #     var,
            #     is_independent=True,
            #     is_binned=True,
            #     units = dvUnits
            # )
            # hdlDVar.values = [(binns[i][j], binns[i][j+1]) for j in range(len(binns[i])-1)]
            # hdlXs = hdl.Variable(
            #     'fiducial higgs to diphoton cross section',
            #     is_independent = False,
            #     is_binned = False,
            #     units = 'fb' if dvUnits is None else 'fb/{}'.format(dvUnits)
            # )
            # hdlXs.values = obsed[i]
            # hdlXsUnc = hdl.Uncertainty("Stat+syst uncertainty on fiducial cross section", is_symmetric=False)
            # hdlXsUnc.values = [(-obsedErrs[i][0][k], obsedErrs[i][1][k]) for k in range(len(obsedErrs[i][0]))]
            # hdlXs.add_uncertainty(hdlXsUnc)
            # table = hdl.Table('{} differential fiducial cross section at CMS'.format(options.extension))
            # table.description = "Differential fiducial higgs to diphoton cross section measured at CMS with respect to {}".format(variable)
            # table.keywords['observables'] = ['DSIG']
            # table.keywords['reactions'] = ['P P --> HIGGS --> GAMMA GAMMA']
            # table.add_image('{}/sigStr_{}_wSyst_Data_commonPredHX.pdf'.format(options.outDir, options.extension))
            # table.add_variable(hdlDVar)
            # table.add_variable(hdlXs)
            # sub = hdl.Submission()
            # sub.add_link('CDS', 'http://cds.cern.ch/record/2803740')
            # sub.add_record_id(2803740, 'CDS')
            # sub.add_table(table)
            # sub.create_files(options.hepDataDir, validate=True, remove_old=True)

        #tops[i].set_ylim(lowYlim, tops[i].get_ylim()[1])


    lowYlims = [top.get_ylim()[0] for top in tops]
    upYlims = [top.get_ylim()[1] for top in tops]

    if options.ylog:
        tops[0].set_yscale('log')
        lowYlims = [top.get_ylim()[0] for top in tops]
        upYlims = [top.get_ylim()[1] for top in tops]
        tops[0].set_ylim(max(1.e-5,min(lowYlims)), 10**(1.6)*(max(upYlims)))
    else:
        tops[0].set_ylim(min(lowYlims), 1.6*max(upYlims))

    bottoms[0].set_ylabel(r'\begin{center}Ratio to\\prediction\end{center}',fontsize=16)
    if isinstance(xlb, tuple):
        tops[0].set_ylabel(r'$\Delta\sigma_{{\mathrm{{fid}}}}/\Delta {{{0}}}$(\mbox{{\ensuremath{{\textnormal{{fb/{{{1}}}}}}}}}\xspace)'.format(xlb[0], xlb[1]))
    else:
        tops[0].set_ylabel(r'$\Delta\sigma_{{\mathrm{{fid}}}}/\Delta {{{0}}}$'.format(xlb) + r'(\mbox{\ensuremath{\textnormal{fb}}}\xspace)') 
        # bottom.set_ylabel(r'$\sigma(\PH\to\PGg\PGg)/\sigma_{SM}(\PH\to\PGg\PGg)$',fontsize=14)
    topLegHandles, topLegLabels = tops[-1].get_legend_handles_labels()
    if options.showSyst:
        legPatch = matplotlib.patches.Patch(facecolor=tuple(topErrColl.get_facecolor()[0]), edgecolor='None', alpha=0.5) 
        topLegHandles.append(legPatch)
        topLegLabels.append(r'syst unc.')
    if options.showStat:
        legPatch = matplotlib.patches.Patch(facecolor=tuple(topErrCollStat.get_facecolor()[0]), edgecolor='None', alpha=0.5) 
        topLegHandles.append(legPatch)
        topLegLabels.append(r'stat unc.')
    if options.showTheoryUnc:
        legPatchesTheoryErr = {}
        ePredKeysHere = copy.deepcopy(ePredKeys)
        for thSkipped in options.thErrToSkip:
            ePredKeysHere.remove(thSkipped)
        for theoryErr in ePredKeysHere+['nominal']:
            legPatch = matplotlib.patches.Patch(edgecolor=topErrCollTheoryUnc[theoryErr].get_edgecolor()[0], facecolor='None', hatch=hatchstr[:int(2*len(hatchstr)/3)], linewidth=0, zorder=3)
            legLine = mlines.Line2D([], [], color=topErrCollTheoryUnc[theoryErr].get_edgecolor()[0], marker='None', zorder=3, linewidth=2)
            topLegHandles[topLegLabels.index(predLabels[theoryErr])] = (legPatch, legLine)
    # if options.commonPredComp is not None:
    #     legPatch = matplotlib.patches.Patch(edgecolor=topcpbColl.get_edgecolor()[0], facecolor='None', linewidth=0, zorder=3) #hatch='xxx',
    #     topLegHandles[topLegLabels.index(predLabels[commonPredKey])] = legPatch
    tops[-1].legend(labels=topLegLabels, handles=topLegHandles, fontsize=14, framealpha=0.8, edgecolor='None')

    try:
        cutString = cutStr[options.extension]
    except KeyError:
        cutString = cutStr[variable]
        
    fig.canvas.draw()
    legdWExt = tops[-1].get_legend().get_window_extent()
    legPosFig = legdWExt.transformed(fig.transFigure.inverted())
    legPosAx = legdWExt.transformed(tops[-1].transAxes.inverted())
    figsize = fig.get_size_inches()*fig.dpi
    xPosCutTxt = 0.02 * (len(tops)/1.)
    cutText = tops[0].text(xPosCutTxt, 0.97, cutString, fontsize=16, va='top', ha='left', transform=tops[0].transAxes)
    textPos = cutText.get_window_extent().transformed(tops[0].transAxes.inverted())
    axPosZero = tops[0].get_position()
    if options.showPvalue:
        pvaltxt = fig.text(legPosFig.x0+0.01, legPosFig.y0-0.005, r'p-value(SM): $\boldsymbol{{{{{:.3f}}}}}$'.format(pvalSMCompat), fontsize=14, va='top', bbox=dict(fc='white', alpha=0.8, edgecolor='None'), transform=fig.transFigure) #(legPos.x0/figsize[0])+0.1, min(textPos.y0/figsize[1]-0.01, legPos.y0/figsize[1]-0.01)
        pvaltxtPos = pvaltxt.get_window_extent().transformed(tops[-1].transAxes.inverted())
        print('pvaltxtPos:', pvaltxtPos)
    if topLvl is not None:
        tlb = xlabels[topVar]
        for j in range(len(topLvl) - 1):
            if topVar == 'Njets2p5' or topVar == 'NBjets2p5' or topVar == 'Nleptons':
                if topLvl[j+1] == 100.:
                    binString = r'$\boldsymbol{{{{{0}}}\geq{{{1:.0f}}}}}$'.format(tlb, topLvl[j] + 0.5)
                else:
                    binString = r'$\boldsymbol{{{{{0}}}={{{1:.0f}}}}}$'.format(tlb, 0.5*(topLvl[j] + topLvl[j+1]))
            else:
                if topLvl[j+1] == 10000.:
                    if isinstance(tlb, tuple):
                        binString = r'$\boldsymbol{{{{{0}}}\,\textnormal{{{2}}}<{{{1}}}}}$'.format(topLvl[j], tlb[0], tlb[1])
                    else:
                        binString = r'$\boldsymbol{{{{{0}}}<{{{1}}}}}$'.format(topLvl[j], tlb)
                elif topLvl[j] == -1000.:
                    if isinstance(tlb, tuple):
                        binString = r'$\boldsymbol{{{{{1}}}<{{{0}}}\,\textnormal{{{2}}}}}$'.format(topLvl[j+1], tlb[0], tlb[1])
                    else:
                        binString = r'$\boldsymbol{{{{{0}}}>{{{1}}}}}$'.format(tlb, topLvl[j+1])
                else:
                    if isinstance(tlb, tuple):
                        binString = r'$\boldsymbol{{{{{0}}}\,\textnormal{{{3}}}<{{{1}}}<{{{2}}}}}\,\textnormal{{{3}}}$'.format(topLvl[j], tlb[0], topLvl[j+1], tlb[1])
                    else:
                        binString = r'$\boldsymbol{{{{{0}}}<{{{1}}}<{{{2}}}}}$'.format(topLvl[j], tlb, topLvl[j+1])
            
            # axPos = tops[j].get_position()
            # yPosBinString = min(textPos.y0/figsize[1]-0.01, legPos.y0/figsize[1]-0.01, pvaltxtPos.y0/figsize[1]-0.01) if options.showPvalue else min(textPos.y0/figsize[1]-0.01, legPos.y0/figsize[1]-0.01)
            # fig.text(axPos.x0 + 0.12/figWidth - axPosZero.x0, yPosBinString, binString, fontsize=16, verticalalignment='top')
            yPosBinString = min(textPos.y0-0.01, legPosAx.y0-0.01, pvaltxtPos.y0-0.01) if options.showPvalue else min(textPos.y0-0.01, legPosAx.y0-0.01)
            tops[j].text(0.1, yPosBinString, binString, fontsize=16, verticalalignment='top', transform=tops[j].transAxes)
    if options.showUnderflow and topLvl is None:
        yPos = max([obsed[i][0]+obsedErrs[i][1][0]] + [pred[i][0]] + [extraPred[k][i][0] for k in range(len(extraPred))])
        if isinstance(xlb, tuple):
            fig.text(bcen[0][0]-xerrs[0][0], 1.1*yPos, r'$\boldsymbol{{\sigma_{{\mathrm{{fid}}}}\left({{{0}}}<{{{1}}}\,\textnormal{{{2}}}\right)/{{{3}}}\,\textnormal{{{2}}}}}$'.format(xlb[0], binns[0][1], xlb[1], binns[0][1] - binns[0][0]), fontsize=10, transform=tops[0].transData) #, horizontalalignment='center')
        else:
            fig.text(bcen[0][0]-xerrs[0][0], 1.1*yPos, r'$\boldsymbol{{\sigma_{{\mathrm{{fid}}}}\left({{{0}}}<{{{1}}}\right)/{{{2}}}}}$'.format(xlb, binns[0][1], binns[0][1] - binns[0][0]), fontsize=10, transform=tops[0].transData) #, horizontalalignment='center')
    if options.cmsText is not None:
        drawIntLumi(tops[-1], intL=137)
    drawCMSLogo(tops[0], opt=options.cmsText)
    
    suff = '_cap' if options.cap is not None else ''
    if options.xlog:
        suff += '_xlog'
    if options.ylog:
        suff += '_ylog'
    if options.showSyst:
        suff += '_wSyst'
    if options.showStat:
        suff += '_wStat'
    if options.suffix is not None:
        suff += options.suffix
    if any(['powhegVBFdipoleRecoil' in predF for predF in options.extraPredFiles]):
        suff += '_wVBFdipR'
    if options.commonPredComp:
        suff += '_commonPred{}'.format(commonPredLDict[commonPredKey])
    fig.savefig('{}/sigStr_{}{}.png'.format(options.outDir, options.extension, suff), bbox_inches='tight')
    fig.savefig('{}/sigStr_{}{}.pdf'.format(options.outDir, options.extension, suff), bbox_inches='tight')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    requiredArgs = parser.add_argument_group()
    requiredArgs.add_argument(
        '--directory', '-d', action='store', type=str, required=True)
    requiredArgs.add_argument(
        '--extension', '-e', action='store', type=str)
    # requiredArgs.add_argument(
    #     '--variable', '-v', action='store', type=str, required=True)
    requiredArgs.add_argument(
        '--outDir', '-o', action='store', type=str, required=True)
    requiredArgs.add_argument(
        '--predFile', '-p', action='store', type=str, required=True)
    optionalArgs = parser.add_argument_group()
    optionalArgs.add_argument(
        '--filename', '-f', action='store', type=str)
    optionalArgs.add_argument(
        '--filenameStat', action='store', type=str)
    optionalArgs.add_argument(
        '--extraPredFiles', nargs='+', action='store', type=str)
    optionalArgs.add_argument(
        '--commonPredComp', action='store', type=str)
    optionalArgs.add_argument(
        '--cap', '-c', action='store', type=float)
    optionalArgs.add_argument(
        '--ylog', action='store_true', default=False)
    optionalArgs.add_argument(
        '--xlog', action='store_true', default=False)
    optionalArgs.add_argument(
        '--showSyst', action='store_true', default=False)
    optionalArgs.add_argument(
        '--showTheoryUnc', action='store_true', default=False)
    optionalArgs.add_argument(
        '--thErrToSkip', nargs='+', action='store', default=[], type=str)
    optionalArgs.add_argument(
        '--suffix', '-s', action='store', type=str)
    optionalArgs.add_argument(
        '--predToSkipRatio', nargs='+', action='store', default=[], type=str)
    optionalArgs.add_argument(
        '--thErrFiles', nargs='+', action='store', type=str)
    optionalArgs.add_argument(
        '--cmsText', action='store', type=str, default=False)
    optionalArgs.add_argument(
        '--showStat', action='store_true', default=False)
    optionalArgs.add_argument(
        '--showUnderflow', action='store_true', default=False)
    optionalArgs.add_argument(
        '--showPvalue', action='store_true', default=False)
    optionalArgs.add_argument(
        '--fnameSMCompat', action='store', type=str)
    optionalArgs.add_argument('--yLim', nargs='+', required=False, type=float)
    optionalArgs.add_argument('--mergeBins', nargs='+', required=False, type=int)
    optionalArgs.add_argument('--dumpResults', action='store_true', default=False, required=False)

    options = parser.parse_args()
    main(options)
