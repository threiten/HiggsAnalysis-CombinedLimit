import root_pandas
import scipy.interpolate as itr
import scipy.optimize as opt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from matplotlib.ticker import MultipleLocator, NullLocator
import pickle as pkl
import oyaml as yaml
import argparse
import copy
from scipy.stats import chi2
from plotBinnedSigStr import xlabels, varDic, lastBins, find1Sig, getNLL, getBins, getPred, ePredLDict
from plotInclusiveScan import findUpDown

def main(options):

    rcP = {'text.usetex': True,
           'font.family': 'sans-serif',
           'font.sans-serif': ['Helvetica'],
           'pdf.fonttype': 42,
           'axes.labelsize': 16,
           'font.size': 16,
           'pgf.rcfonts': True,
           'text.latex.preamble': r"\usepackage{bm, xspace, amsmath}"}

    plt.rcParams.update(rcP)
    predColors = list(cm.Dark2.colors)
    systColor = 'dodgerblue'
    with open(options.configFile, 'r') as f:
        config = yaml.load(f)
        f.close()

    print(config)

    with open(config['global']['totalXSFile'], 'rb') as f:
        totalXSDic = pkl.load(f)
        f.close()
    with open(config['global']['fidXSFile'], 'rb') as g:
        fidXSDic = pkl.load(g)
        g.close()

    accNom = fidXSDic['nominal']/totalXSDic['nominal']
    accUncDownPdf, accUncUpPdf = findUpDown('pdfWeights', range(1,60), accNom, fidXSDic, totalXSDic)
    accUncDownScale, accUncUpScale = findUpDown('scaleWeights', [4,8], accNom, fidXSDic, totalXSDic)
    accUncDownAlpS, accUncUpAlpS = findUpDown('alphaSWeights', range(2), accNom, fidXSDic, totalXSDic)
    accUncDown = np.sqrt(((accNom - accUncDownPdf)/accNom)**2 + ((accNom - accUncDownScale)/accNom)**2 + ((accNom - accUncDownAlpS)/accNom)**2)
    accUncUp = np.sqrt(((accUncUpPdf - accNom)/accNom)**2 + ((accUncUpScale - accNom)/accNom)**2 + ((accUncUpAlpS - accNom)/accNom)**2)
    accUnc = np.array([accUncDown, accUncUp])
    
    cvals = {}
    errs = {}
    preds = {}
    errsPred = {}
    errsTh = {}
    errsStat = {}
    errsSyst = {}
    for key, itm in config.items():
        if key == 'global':
            continue
        print(itm)
        binns, topLvl, topVar, _ = getBins(varDic[itm['extension']], itm['extension'])
        for k in range(len(binns)):
            print('Bins before:', binns[k])
            if binns[k][-1] > 2 * binns[k][-2]:
                binns[k][-1] = 2 * binns[k][-2] - binns[k][-3] if len(binns[k])>2 else binns[k][-1]
            if varDic[itm['extension']] in lastBins.keys():
                binns[k][-1] = lastBins[varDic[itm['extension']]] 
            print('Bins after:', binns[k])
        df = root_pandas.read_root(itm['fNameFit'], 'limit')
        sStr, deltaNLL = getNLL('r{}'.format(itm['bin']), df)
        dfStat = root_pandas.read_root(itm['fNameFitStat'], 'limit')
        sStrStat, deltaNLLStat = getNLL('r{}'.format(itm['bin']), dfStat)
        cvals[key] = sStr[0]
        errs[key] = np.array(find1Sig(sStr, deltaNLL))
        errsStat[key] = np.array(find1Sig(sStrStat, deltaNLLStat, cvals[key]))
        errsSyst[key] = np.sqrt(np.power(errs[key], 2) - np.power(errsStat[key], 2))
        with open(itm['predFile'], 'rb') as f:
            predF = pkl.load(f)
            f.close()
        pred = np.array([predF[itm['bin']]])
        mass = df['mh'].values[0]
        preds[key] = (binns[0][itm['bin']+1] - binns[0][itm['bin']]) * getPred(pred, mass, 1, weights=[1., 2.3, 1.])[0]
        errsPred[key] = preds[key] * errs[key]
        errsSyst[key] = preds[key] * errsSyst[key]
        with open(itm['thErrFile'], 'rb') as f:
            thErr = pkl.load(f)[:, itm['bin']]
            f.close()
        print('thErr abs: ', thErr)
        thErr = np.divide((binns[0][itm['bin']+1] - binns[0][itm['bin']]) * thErr, preds[key])
        print('thErr rel: ', thErr)
        print('pred :', preds[key])
        errsTh[key] = np.sqrt(np.power(accUnc, 2)+np.power(thErr, 2)+np.power(config['global']['BRUnc']*np.ones_like(thErr), 2)) * preds[key]
        print('total thErr: ', errsTh[key])
    print('cvals: ', cvals)
    print('errs: ', errs)
    print('preds: ', preds)
    print('errsPred: ', errsPred)
    print('errsTh: ', errsTh)

    fig, axes = plt.subplots(1, 1, figsize=(7+len(config.keys()), 7), squeeze=False)
    fig.tight_layout()
    pt = fig.axes[0]
    hatchstr = 11*'/'

    cvPlot = []
    tckLbls = []
    errsPlot = []
    predsPlot = []
    errsThPlot = []
    errsSystPlot = []
    yPoss = []
    offs = 0
    for i, key in enumerate(config.keys()):
        if key == 'global':
            offs = -1
            continue
        cvPlot.append(cvals[key] * preds[key])
        errsPlot.append(errsPred[key])
        predsPlot.append(preds[key])
        tckLbls.append(config[key]['label'])
        errsThPlot.append(errsTh[key])
        errsSystPlot.append(errsSyst[key])
        yPoss.append(i+offs)
    cvPlot = np.array(cvPlot)
    errsPlot = np.abs(np.array(errsPlot).T)
    errsSystPlot = np.array(errsSystPlot)
    errsThPlot = np.array(errsThPlot)
    yPoss = np.array(yPoss)
    print(cvPlot)
    print(yPoss)
    print(errsSystPlot)
    (_,caps,_) = pt.errorbar(yPoss, cvPlot, yerr=errsPlot, marker='.', linestyle='none', markersize=8, capsize=3, linewidth=1, color='black', label=r'Data, stat$\oplus$syst unc.')
    for cap in caps:
        cap.set_markeredgewidth(0)
    for i, prd in enumerate(predsPlot):
        pt.plot([yPoss[i]-0.4, yPoss[i]+0.4], [prd, prd], '-', color=predColors[0], label=ePredLDict['nominal'])
    for i, cv in enumerate(cvPlot):
        pt.text(yPoss[i], (cv+errsPlot[1, i])*1.1, r'${{{0:.2f}}}_{{\,-{1:.2f}}}^{{\,+{2:.2f}}}$'.format(cv, errsPlot[0, i], errsPlot[1, i]), ha='center', va='bottom')
    errorboxesTh = [matplotlib.patches.Rectangle((x-xe[0], y-ye[0]), xe.sum(), ye.sum()) for x, y, xe, ye in zip(yPoss, predsPlot, 0.4*np.ones_like(errsThPlot), errsThPlot)]
    print(errorboxesTh[0])
    pcBottom = matplotlib.collections.PatchCollection(errorboxesTh, facecolor='None', alpha=0.85, hatch=hatchstr, edgecolor=predColors[0], linewidth=0, zorder=3)
    pt.add_collection(pcBottom)
    systBoxes = [matplotlib.patches.Rectangle((x-xe[0], y-ye[0]), xe.sum(), ye.sum()) for x, y, xe, ye in zip(yPoss, cvPlot, 0.4*np.ones_like(errsSystPlot), errsSystPlot)]
    print(systBoxes)
    pcSystBoxes = matplotlib.collections.PatchCollection(systBoxes, facecolor=systColor, alpha=0.5, edgecolor='none', zorder=3)
    pt.add_collection(pcSystBoxes)
    pt.legend()
    LegHandles, LegLabels = pt.get_legend_handles_labels()
    lblCache = []
    LegLabelsUpd = []
    LegHandlesUpd = []
    for i in range(len(LegLabels)):
        if LegLabels[i] not in lblCache:
            LegLabelsUpd.append(LegLabels[i])
            LegHandlesUpd.append(LegHandles[i])
            lblCache.append(LegLabels[i])
    legPatch = matplotlib.patches.Patch(edgecolor=predColors[0], facecolor='None', hatch=hatchstr[:int(len(hatchstr)/2)], linewidth=0, zorder=3)
    legLine = mlines.Line2D([], [], color=predColors[0], marker='|', linestyle='None', zorder=3, markeredgewidth=1, markersize=10)
    legPatchSyst = matplotlib.patches.Patch(edgecolor='None', facecolor=systColor, alpha=0.5, zorder=3)
    LegHandlesUpd.append(legPatchSyst)
    LegLabelsUpd.append(r'syst unc.')
    LegHandlesUpd[LegLabelsUpd.index(ePredLDict['nominal'])] = (legPatch, legLine)
    pt.legend(labels=LegLabelsUpd, handles=LegHandlesUpd, fontsize=14, framealpha=0)
    pt.set_xlim((-1,yPoss.shape[0]))
    pt.xaxis.set_major_locator(NullLocator())
    pt.xaxis.set_minor_locator(MultipleLocator(2))
    pt.set_xticks(yPoss)
    pt.set_xticklabels(tckLbls, rotation=50)
    # pt.invert_xaxis()
    pt.set_yscale('log')
    pt.set_ylabel(r'$\boldsymbol{\sigma_{fid}}\left(\textnormal{fb}\right)$', loc='top')
    pt.set_ylim((pt.get_ylim()[0], 10**1.3*np.max(cvPlot+errsPlot[1,:])))
    fig.text(0.04, 0.965, r'\textbf{CMS} \textit{Preliminary}', fontsize=18)
    fig.text(0.87, 0.965, r'138\mbox{\ensuremath{\,\text{fb}^{-1}}}\xspace (13\ensuremath{\,\text{Te\hspace{-.08em}V}}\xspace)', fontsize=18)
    fig.savefig('/eos/home-t/threiten/www/plots/Hgg/plotFidXS_flipped.png', bbox_inches='tight')
    fig.savefig('/eos/home-t/threiten/www/plots/Hgg/plotFidXS_flipped.pdf', bbox_inches='tight')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    requiredArgs = parser.add_argument_group()
    requiredArgs.add_argument(
        '--outDir', '-o', action='store', type=str, required=True)
    requiredArgs.add_argument(
        '--configFile', '-c', action='store', type=str, required=True)
    # optionalArgs = parser.add_argument_group()
    # optionalArgs.add_argument(
    #     '--filename', '-f', action='store', type=str)
    # optionalArgs.add_argument(
    #     '--filenameStat', action='store', type=str)
    # optionalArgs.add_argument(
    #     '--extraPredFiles', nargs='+', action='store', type=str)
    # optionalArgs.add_argument(
    #     '--commonPredComp', action='store', type=str)
    # optionalArgs.add_argument(
    #     '--cap', '-c', action='store', type=float)
    # optionalArgs.add_argument(
    #     '--ylog', action='store_true', default=False)
    # optionalArgs.add_argument(
    #     '--xlog', action='store_true', default=False)
    # optionalArgs.add_argument(
    #     '--showSyst', action='store_true', default=False)
    # optionalArgs.add_argument(
    #     '--showTheoryUnc', action='store_true', default=False)
    # optionalArgs.add_argument(
    #     '--thErrToSkip', nargs='+', action='store', default=[], type=str)
    # optionalArgs.add_argument(
    #     '--suffix', '-s', action='store', type=str)
    # optionalArgs.add_argument(
    #     '--predToSkipRatio', nargs='+', action='store', default=[], type=str)
    # optionalArgs.add_argument(
    #     '--thErrFiles', nargs='+', action='store', type=str)
    # optionalArgs.add_argument(
    #     '--prelim', action='store_true', default=False)
    # optionalArgs.add_argument(
    #     '--wip', action='store_true', default=False)
    # optionalArgs.add_argument(
    #     '--showStat', action='store_true', default=False)
    # optionalArgs.add_argument(
    #     '--showUnderflow', action='store_true', default=False)
    # optionalArgs.add_argument(
    #     '--showPvalue', action='store_true', default=False)
    # optionalArgs.add_argument(
    #     '--fnameSMCompat', action='store', type=str)
    # optionalArgs.add_argument('--yLim', nargs='+', required=False, type=float)
    # optionalArgs.add_argument('--mergeBins', nargs='+', required=False, type=int)

    options = parser.parse_args()
    main(options)
