#!/usr/bin/env python

###
# YoSon
# 08/20/2015
# @ run_rasqual/plot/plotrasqual.py
###

import sys, os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as transforms
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib.ticker import FuncFormatter
import pylab
from scipy import log10

from plotpanels import *
from plotparams import *

if __name__ == '__main__':
    '''check input data and output fig name'''


chrom = str(sys.argv[1])
d = str(sys.argv[2])

data, pos, phi, delta, pval, r2_rsnps = read_rqout(d)
logp, maxpos, rectwidth, maxgene, maxgpos, maxx, maxy, annotx, annoty = setlims(data)

fig, ax1, ax2, ax3, ax4 = plotpanels(data)
saveplot(fig, ax1, ax2, ax3, ax4)



