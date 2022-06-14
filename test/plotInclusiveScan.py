import root_pandas
import scipy.interpolate as itr
import scipy.optimize as opt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
import pickle as pkl
import oyaml as yaml
import argparse
import copy
from plotBinnedSigStr import getNLL, getBins, getPred, find1Sig, xlabels, ePredLDict, varDic, cutStr, drawCMSLogo, drawIntLumi


legLabel = {
    'nominal' : r'\textsc{MG}\xspace{}5\_a\textsc{mc@nlo}, \textsc{nnlops}'
}

def findUpDown(name, rge, nomVal, fidDic, totDic):
    accDown = nomVal
    accUp = nomVal
    for i in rge:
        accVar = fidDic['{}{}'.format(name, i)]/totDic['{}{}'.format(name, i)]
        if accVar < accDown:
            accDown = accVar
        if accVar > accUp:
            accUp = accVar

    return accDown, accUp


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
    predColors = list(cm.Dark2.colors)

    variable = varDic[options.extension]
    binns, topLvl, topVar, _ = getBins(variable, options.extension)

    nBins = len(binns)

    fname = 'higgsCombineAsimovPostFitScanFit' if options.filename is None else options.filename
    errs = np.zeros((2,nBins))
    cval = np.zeros((nBins))
    for n in range(nBins):
        par = 'r{}'.format(n)
        df = root_pandas.read_root('{}/{}_{}.MultiDimFit.mH125.38.root'.format(options.directory, fname, par), 'limit')
        sStr, deltaNLL = getNLL(par, df)
        cval[n] = sStr[0]
        neg, pos, NLLSpl = find1Sig(sStr, deltaNLL, retSpl=True)
        errs[:, n] = np.array((neg, pos))

    errs = np.abs(errs)
    with open(options.predFile, 'rb') as f:
        preds = pkl.load(f)
        f.close()

    mass = df['mh'].values[0]
    pred = getPred(preds, mass, nBins)
    pred *= 10000.

    if options.showTheoryUnc:
        with open(options.totalXSFile, 'rb') as f:
            totalXSDic = pkl.load(f)
            f.close()
        with open(options.fidXSFile, 'rb') as g:
            fidXSDic = pkl.load(g)
            g.close()

        accNom = fidXSDic['nominal']/totalXSDic['nominal']
        accUncDownPdf, accUncUpPdf = findUpDown('pdfWeights', range(1,60), accNom, fidXSDic, totalXSDic)
        accUncDownScale, accUncUpScale = findUpDown('scaleWeights', [4,8], accNom, fidXSDic, totalXSDic)
        accUncDownAlpS, accUncUpAlpS = findUpDown('alphaSWeights', range(2), accNom, fidXSDic, totalXSDic)

        uncDown = np.sqrt(((accNom - accUncDownPdf)/accNom)**2 + ((accNom - accUncDownScale)/accNom)**2 + ((accNom - accUncDownAlpS)/accNom)**2)
        uncUp = np.sqrt(((accUncUpPdf - accNom)/accNom)**2 + ((accUncUpScale - accNom)/accNom)**2 + ((accUncUpAlpS - accNom)/accNom)**2)
        if options.BRunc is not None:
            uncDown = np.sqrt(uncDown**2 + options.BRunc**2)
            uncUp = np.sqrt(uncUp**2 + options.BRunc**2)
        print(uncDown, uncUp)

    fig = plt.figure(figsize=(8,8))
    fig.tight_layout()
    pt = fig.add_subplot(111)

    xrge = options.xrge if options.xrge is not None else 0.2
    xPred = np.linspace((cval[0]-xrge)*pred[0], (cval[0]+xrge)*pred[0], 1000)
    yNLL  = np.array([NLLSpl(x/pred[0]) for x in xPred])

    print('Predicted xs: {} +{} -{}'.format(pred[0], uncUp*pred[0], uncDown*pred[0]))
    pt.plot(xPred, yNLL, '-', color='black', linewidth=3, zorder=6)
    ylim = (0., 1.2*pt.get_ylim()[1])
    pt.plot([pred[0], pred[0]], [ylim[0], 0.9*ylim[1]], '-', color='red', label=legLabel['nominal'], linewidth=3)
    pt.set_ylim(ylim)
    pt.plot([(cval[0]-xrge)*pred[0], (cval[0]+xrge)*pred[0]], [1.,1.], '--', alpha=0.8, color='lightslategray', linewidth=3, zorder=5)
    pt.plot([pred[0]*(cval[0]-errs[0,0])]*2, [0., 1.], '--', alpha=0.8, color='lightslategray', linewidth=3, zorder=5)
    pt.plot([pred[0]*(cval[0]+errs[1,0])]*2, [0., 1.], '--', alpha=0.8, color='lightslategray', linewidth=3, zorder=5)
    pt.set_xlim([(cval[0]-xrge)*pred[0], (cval[0]+xrge)*pred[0]])
    pt.set_xlabel(r'\ensuremath{\boldsymbol{\sigma}_{\text{fid}}(\text{fb})}') #, loc='right')
    pt.set_ylabel(r'\ensuremath{-2\Delta\ln \text{L}}') #, loc='top')
    # pt.grid()
    try:
        cutString = cutStr[options.extension]
    except KeyError:
        cutString = cutStr[variable]
        
    topLegHandles, topLegLabels = pt.get_legend_handles_labels()
    if options.showTheoryUnc:
        hatchstr = '/'*5
        errPatch = matplotlib.patches.Rectangle((pred[0]*(1.-uncDown), ylim[0]), pred[0]*(uncUp+uncDown), ylim[0]+0.9*ylim[1], edgecolor='red', facecolor='None', hatch=hatchstr, linewidth=0, zorder=3)
        pt.add_patch(errPatch)

        legPatch = matplotlib.patches.Patch(edgecolor=errPatch.get_edgecolor(), facecolor='None', hatch=hatchstr, linewidth=0, zorder=3)
        legLine = mlines.Line2D([], [], color=errPatch.get_edgecolor(), marker='None', zorder=3, linewidth=1)
        topLegHandles[topLegLabels.index(legLabel['nominal'])] = (legPatch, legLine)
    pt.legend(labels=topLegLabels, handles=topLegHandles, fontsize=14, framealpha=0)
    fig.canvas.draw()
    
    legPos = pt.get_legend().get_window_extent()
    figsize = fig.get_size_inches()*fig.dpi
    # cutText = fig.text(legPos.x0/figsize[0], legPos.y0/figsize[1], cutString, fontsize=16, verticalalignment='top')
    cutText = pt.text(0.02, 0.97, cutString, fontsize=16, verticalalignment='top', transform=pt.transAxes)
    txtPos = cutText.get_window_extent().transformed(pt.transAxes.inverted())
    if legPos.x0/figsize[0] < 0.5:
        pt.legend(labels=topLegLabels, handles=topLegHandles, fontsize=14, framealpha=0, bbox_to_anchor=(txtPos.x0, txtPos.y0), loc='upper left', borderpad=0., handletextpad=0.1, borderaxespad=0.2, columnspacing=0.)
    if options.cmsText is not None:
        drawCMSLogo(pt, opt=options.cmsText)
    drawIntLumi(pt, intL=137)
    fig.savefig('{}/NLL_{}.png'.format(options.outDir, options.extension), bbox_inches='tight')
    fig.savefig('{}/NLL_{}.pdf'.format(options.outDir, options.extension), bbox_inches='tight')
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    requiredArgs = parser.add_argument_group()
    requiredArgs.add_argument(
        '--directory', '-d', action='store', type=str, required=True)
    requiredArgs.add_argument(
        '--outDir', '-o', action='store', type=str, required=True)
    requiredArgs.add_argument(
        '--extension', '-e', action='store', type=str)
    requiredArgs.add_argument(
        '--predFile', '-p', action='store', type=str, required=True)
    optionalArgs = parser.add_argument_group()
    optionalArgs.add_argument(
        '--filename', '-f', action='store', type=str)
    optionalArgs.add_argument(
        '--filenameStat', action='store', type=str)
    optionalArgs.add_argument(
        '--showSyst', action='store_true', default=False)
    optionalArgs.add_argument(
        '--showTheoryUnc', action='store_true', default=False)
    optionalArgs.add_argument(
        '--fidXSFile', action='store', type=str)
    optionalArgs.add_argument(
        '--totalXSFile', action='store', type=str)
    optionalArgs.add_argument(
        '--xrge', action='store', type=float)
    optionalArgs.add_argument(
        '--BRunc', action='store', type=float)
    optionalArgs.add_argument(
        '--showStat', action='store_true', default=False)
    optionalArgs.add_argument('--yLim', nargs='+', required=False, type=float)
    optionalArgs.add_argument('--cmsText', action='store', type=str)

    options = parser.parse_args()
    main(options)
