#!/usr/bin/evn python

###
# YoSon
# 08/20/2015
# @ run_rasqual/plot/run_test.py
###

from plotimports import *
from plotpanels import *
from plotparams import *
from plotsave import *

#chrom = 'chr1'

chrom = sys.argv[1]
d = sys.argv[2]

data, pos, phi, delta, pval, r2_rsnps = read_rqout(d)

logp, maxpos, rectwidth, maxgene, maxgpos, maxx, maxy, annotx, annoty = setlims(data)

fig, ax1, ax2, ax3, ax4 = plotpanels(data)
saveplot(chrom, fig, ax1, ax2, ax3, ax4)





