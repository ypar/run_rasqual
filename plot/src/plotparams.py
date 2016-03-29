#!/usr/bin/env python

###
# YoSon
# 08/20/2015
# @ run_rasqual/plot/plotparams.py
###

from plotpanels import *
from plotparams import *
from plotsave import *

from plotimports import *

if __name__ == '__main__':
    '''check input data and output fig name'''


def read_rqout(d):
    # read rqout with pvalues
    data = pd.read_csv(d, sep='\t', header=0, compression='gzip')
    # pos = data[data.columns[3]]
    # prev release of rasqual had phi and delta columns switched
    # in case earlier results contain such info, call in by col location
    phi = data[data.columns[13]]
    delta = data[data.columns[12]]
    pos = data.pos
    #phi = data.ref_bias_phi
    #delta = data.seq_err_rate_delta
    r2_rsnps = data.r2_rsnps
    pval = data.chisq_pval
    return (data, pos, phi, delta, pval, r2_rsnps)


def mbps(x, pos):
    #The two args are the value and tick position
    return '%1.1fMB' % (x*1e-6)


def simple(ax):
    # simplify plots by eliminating reduncant axis, etc
    # also add formatting for chrom bp into megabase positions
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(direction='out', pad=5)
    ax.xaxis.set_major_formatter(FuncFormatter(mbps))
    return ax


def setax(ax, maxpos, ty):
    ax.set_xlim(left=-5)
    ax.set_xlim(right=maxpos)
    ax.set_ylim(bottom=0)
    ax.set_ylim(top=int(ty))
    return ax


def draw_rect(ax, spos, rectwidth):
    # find the most significant gene and draw a box around it
    trans = transforms.blended_transform_factory(
    ax.transData, ax.transAxes)
    rect = mpatches.Rectangle((int(spos), 0), width=rectwidth, height=1,
                         transform=trans, color='grey',
                         alpha=0.2)
    return rect


def setlims(data):
    # get some set values for plotting purposes
    logp = -1 * log10(data.chisq_pval)
    maxpos = int(data.pos.max() + 2e6)
    rectwidth = int(int(maxpos) / 100)
    maxgene = data[data.chisq_pval == data.chisq_pval.min()].max().feature_id
    maxgpos = int(data[data.chisq_pval == data.chisq_pval.min()].max().pos)
    maxx = int(maxgpos)
    maxy = -log10(data.chisq_pval.min())
    annotx = int(maxx) + 5e5
    annoty = int(maxy) + 10
    return (logp, maxpos, rectwidth, maxgene, maxgpos, maxx, maxy, annotx, annoty)



