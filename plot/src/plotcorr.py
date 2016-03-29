#!/usr/bin/env python

###
# YoSon
# 08/20/2015
# @ run_rasqual/plot/plotcorr.py
###

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from plotpanels import *
from plotparams import *
from plotsave import *

from plotimports import *


def read_allchr(d):
    # read all chr rqout with pvalues
    data = pd.read_csv(d, sep='\t', header=None, compression='gzip', names=['feature_id', 'rs_id', 'chrom', 'pos', 'ref', 'alt', 'allele_freq', 'hwe_chisq', 'impute_qual', 'impute_qua_rsq', 'chisq', 'effect_size', 'seq_err_rate_delta', 'ref_bias_phi', 'overdisp', 'snp_id_region', 'no_feature_snps', 'no_tested_snps', 'no_iter_h0', 'no_iter_h1', 'rand_loc_ties', 'likelihood_h0', 'convergence', 'r2_fsnps', 'r2_rsnps', 'chisq_pval'])
    return (data)


def allchr(wkdir, suffix):
    
    heads = ['feature_id', 'rs_id', 'chrom', 'pos', 'ref', 'alt', 'allele_freq', 'hwe_chisq', 'impute_qual', 'impute_qua_rsq', 'chisq', 'effect_size', 'seq_err_rate_delta', 'ref_bias_phi', 'overdisp', 'snp_id_region', 'no_feature_snps', 'no_tested_snps', 'no_iter_h0', 'no_iter_h1', 'rand_loc_ties', 'likelihood_h0', 'convergence', 'r2_fsnps', 'r2_rsnps', 'chisq_pval']
    
    chr = str(1)
    chrom = "chr" + str(chr)
    o = pd.DataFrame(read_rqout(filename))
    o.info()
    for chr in range(2,23):
        chr = str(chr)
        chrom = "chr" + str(chr)
        temp = pd.DataFrame(read_rqout(ofile))
        final = pd.concat([o, temp])
        o.info()
        temp.info()
        final.info()
        o = final
        o.columns = heads
        del final
        
    outfile = wkdir + "temp_allchr" + suffix
    o.to_csv(outfile, sep="\t", index=False)
    
    return o
    



if __name__ == '__main__':
    '''checks pval, permutation pval, and null model pvals from rasqual runs to plot correlations'''


wkdir = sys.argv[1]
ofile = wkdir + sys.argv[2]
xfile = wkdir + sys.argv[3]
nfile = wkdir + sys.argv[4]

odata = read_allchr(ofile)
odata = odata[odata.ref != "X"]
xdata = read_allchr(xfile)
xdata = xdata[xdata.ref != "X"]
ndata = read_allchr(nfile)
ndata = ndata[ndata.ref != "X"]

df12 = odata.merge(xdata, left_on=['chrom', 'feature_id', 'rs_id'], right_on=['chrom', 'feature_id', 'rs_id'])
temp = df12[['chisq_pval_x', 'chisq_pval_y']]
df12 = temp
temp.to_csv('temp_allchr_tx_a01_pval_xpval.out', sep='\t')
df12['chisq_pval_x'].corr(df12['chisq_pval_y'])

df13 = odata.merge(ndata, left_on=['chrom', 'feature_id', 'rs_id'], right_on=['chrom', 'feature_id', 'rs_id'])
temp = df13[['chisq_pval_x', 'chisq_pval_y']]
df13 = temp
temp.to_csv('temp_allchr_tx_a01_pval_npval.out', sep='\t')
df13['chisq_pval_x'].corr(df13['chisq_pval_y'])

df23 = xdata.merge(ndata, left_on=['chrom', 'feature_id', 'rs_id'], right_on=['chrom', 'feature_id', 'rs_id'])
temp = df23[['chisq_pval_x', 'chisq_pval_y']]
df23 = temp
temp.to_csv('temp_allchr_tx_a01_xpval_npval.out', sep='\t')
df23['chisq_pval_x'].corr(df23['chisq_pval_y'])


p = pd.DataFrame()
p['x'] = -log10(odata['chisq_pval_x'])
p['y'] = -log10(odata['chisq_pval_y'])

def simple(ax):
    # simplify plots by eliminating reduncant axis, etc
    # also add formatting for chrom bp into megabase positions
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(direction='out', pad=5)
    return ax


plt.figure(1, figsize=(8,8))
fig, ax = plt.subplots()
fig.set_size_inches(8.5, 8.5, forward=True)
gs = gridspec.GridSpec(4, 4)
plt.scatter(p.x, p.y, color='#FF4500', alpha=0.3, s=15)

pcorr = p.corr().y[0]
plt.text(330, 330,'$R^2 = %0.5f$'% pcorr, fontsize=30)



plt.scatter(p.x,p.y)

ax = simple(ax)
ax.set_title('ASE analysis results \n with and without feature count svd covariates', fontsize=14, color='0.3')
ax.set_xlabel('unadjusted pvalues', fontsize=12)
ax.set_ylabel('covariate adjusted pvalues', fontsize=12)
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
figfile = ofile + ".png"
plt.savefig(figfile, dpi=300)



###
par = np.polyfit(p.x, p.y, 1, full=True)
par = np.polyfit(p.x, p.y, 1)
slope=par[0][0]
intercept=par[0][1]
xl = [min(p.x), max(p.x)]
yl = [slope*xx + intercept  for xx in xl]
variance = np.var(p.y)
residuals = np.var([(slope*xx + intercept - yy)  for xx,yy in zip(p.x, p.y)])
Rsqr = np.round(1-residuals/variance, decimals=2)
plt.text(.9*max(p.x)+.1*min(p.x),.9*max(p.y)+.1*min(p.y),'$R^2 = %0.5f$'% Rsqr, fontsize=30)
###



nullfmt   = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left+width+0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(1, figsize=(8,8))

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
axScatter.scatter(x, y)

# now determine nice limits by hand:
binwidth = 0.25
xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
lim = ( int(xymax/binwidth) + 1) * binwidth

axScatter.set_xlim( (-lim, lim) )
axScatter.set_ylim( (-lim, lim) )

bins = np.arange(-lim, lim + binwidth, binwidth)
axHistx.hist(x, bins=bins)
axHisty.hist(y, bins=bins, orientation='horizontal')

axHistx.set_xlim( axScatter.get_xlim() )
axHisty.set_ylim( axScatter.get_ylim() )

plt.show()



