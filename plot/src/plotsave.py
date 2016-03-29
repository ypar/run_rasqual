#!/usr/bin/env python

###
# YoSon
# 08/20/2015
# @ run_rasqual/plot/plotsave.py
###

from plotpanels import *
from plotparams import *
from plotsave import *

from plotimports import *


if __name__ == '__main__':
    '''check the wrapper'''


def saveplot(chrom, fig, ax1, ax2, ax3, ax4):
    # adjust panels and save figure
    chromlabel = chrom + ' position (MB)'
    plt.xlabel(chromlabel, fontsize=16, color='0.4')
    plt.subplots_adjust(left=0.125, bottom=0.125, hspace=0.4)
    figfile = chrom + '_rqout_pval_delta_phi_r2rsnp.png'
    plt.savefig(figfile, dpi=300)



