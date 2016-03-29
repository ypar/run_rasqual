#!/usr/bin/env python

###
# YoSon
# 08/20/2015
# @ run_rasqual/plot/plotpanels.py
###

from plotpanels import *
from plotparams import *
from plotsave import *

from plotimports import *


if __name__ == '__main__':
    '''check the wrapper'''



def plotpanels(data):
    phi = data[data.columns[13]]
    delta = data[data.columns[12]]
    pos = data.pos
    r2_rsnps = data.r2_rsnps
    pval = data.chisq_pval
    logp = -1 * log10(data.chisq_pval)
    
    #scatter plot of -log10 pvals
    fig, ax1 = plt.subplots()
    fig.set_size_inches(14.5, 8.5, forward=True)
    gs = gridspec.GridSpec(4, 4)
    
    fig.suptitle('ASE analysis results and estimated parameters', fontsize=20, color='0.3')
    logp, maxpos, rectwidth, maxgene, maxgpos, maxx, maxy, annotx, annoty = setlims(data)
    
    ax1 = plt.subplot(gs[0, :])
    plt.scatter(pos, logp, color='#1e90ff', alpha=0.5, s=15)
    ax1 = simple(ax1)
    #annotate the most significant gene on the pval plot
    ax1.annotate(maxgene, xy=(maxx, maxy), xytext=(annotx, annoty), arrowprops=dict(arrowstyle="-", facecolor='grey'))
    plt.gca().set_xlim(left=-5)
    plt.gca().set_xlim(right=maxpos)
    plt.ylabel(r'$\mathit{-log(p)}$', fontsize=13, color='0.4')
    plt.gca().set_ylim(bottom=0)
    plt.gcf().subplots_adjust(bottom=0.15)
    # check fdr values for empirical threshold
    pline = -1 * log10(10e-8)
    plt.axhline(y=pline, color = 'r', alpha=0.5)
    ax1.add_patch(draw_rect(ax1, int(maxx-1e5), rectwidth))
    
    #scatter plot of delta
    ax2 = plt.subplot(gs[1, :], sharex=ax1)
    plt.scatter(pos, delta, color='#4b0082', alpha=0.3, s=15)
    ax2 = simple(ax2)
    ax2 = setax(ax2, maxpos, 1.0)
    plt.ylabel(r'$\hat \delta$', fontsize=13, color='0.4')
    plt.axhline(y=0.01, color = 'r', alpha=0.5)
    ax2.add_patch(draw_rect(ax2, int(maxx-1e5), rectwidth))
    
    #scatter plot of phi
    ax3 = plt.subplot(gs[2, :], sharex=ax1)
    plt.scatter(pos, phi, color='#f4a460', alpha=0.3, s=15)
    ax3 = simple(ax3)
    ax3 = setax(ax3, maxpos, 1.0)
    plt.ylabel(r'$\hat \phi$', fontsize=13, color='0.4')
    plt.axhline(y=0.25, color = 'r', alpha=0.5)
    ax3.add_patch(draw_rect(ax3, int(maxx-1e5), rectwidth))
    
    #r2_rsnps
    ax4 = plt.subplot(gs[3, :], sharex=ax1)
    plt.scatter(pos, r2_rsnps, color='#228b22', alpha=0.3, s=15)
    ax4 = simple(ax4)
    ax4 = setax(ax4, maxpos, 1.0)
    plt.axhline(y=0.5, color = 'r', alpha=0.5)
    ax4.add_patch(draw_rect(ax4, int(maxx-1e5), rectwidth))
    plt.ylabel(r'$R^2 \mathit{rsnps}$', fontsize=13, color='0.4')
    
    return (fig, ax1, ax2, ax3, ax4)



