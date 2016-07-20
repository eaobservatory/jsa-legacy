import datetime
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from astropy.table import Table, join, Column


matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{apjfonts},\usepackage{sfmath}'
matplotlib.rcParams['figure.autolayout'] = True
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['legend.fontsize'] = 'medium'
matplotlib.rcParams['axes.linewidth'] = '0.5'



# Calstats info, from legacy release calibrators and also from CAL DB
# regular reductions.
lrlogfile = 'combined_calstats.log'
caldblogfile = 'caldb-arcsecfcfs.txt'

lrvalues = Table.read(lrlogfile, format='ascii', header_start=1)
regularcals = Table.read(caldblogfile, format='ascii')

# Join the two tables.
# First fix up the headers: join on HST and obsnum
regular_hsts = [datetime.datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in regularcals['hstdate']]
regularcals.add_column(Column(regular_hsts, 'HST_FIX'))
lr_hsts = [datetime.datetime.strptime(i, '%Y-%m-%dT%H:%M:%S') for i in lrvalues['HST']]
lrvalues.add_column(Column(lr_hsts, 'HST_FIX'))
regularcals.rename_column('obsnum', 'Obs')

# Only find observations which are in both calibration DB and in legacy release calibrations
join_cals = join(lrvalues, regularcals, keys=('HST_FIX', 'Obs'), join_type='inner')



def show_results(names, Sources, FCFs, FCF_errors):
    clippedresults = {}
    for i in sources:
        if i == 'All':
            mask = Sources!= 'not a source'
        else:
            mask = Sources == i
        fcfs  = FCFs[mask]
        errs = FCF_errors[mask]
        median = np.median(fcfs)
        mean = np.mean(fcfs)
        weighted_mean = np.average(fcfs, weights=1./errs)
        std = np.std(fcfs)
        clipped, lowval, highval = scipy.stats.sigmaclip(fcfs, low=5.0, high=5.0)
        clippedmean = np.mean(clipped)
        clippedstd = np.std(clipped)
        print('{:>7}: median={:.2f}, mean={:.2f}, weighted_mean={:.2f}, std={:.4f}, clipped 5sig mean={:.2f} and std={:.2f} ({:.2f} to {:.2f} )'.format(
                i, median, mean, weighted_mean, std, clippedmean, clippedstd, lowval, highval))
        clippedresults[i] = clippedmean, clippedstd, lowval, highval
    return clippedresults

print('Results from all legacy calibrators')
sources = list(set(lrvalues['Source']))
sources= ['All'] + sources
lrclippedresults = show_results(sources, lrvalues['Source'], lrvalues['FCFasec'], lrvalues['FCFasec_err'])

print('\n')
print('Results from cal DB values with counterpart in legacy calibrators')

show_results(sources, join_cals['Source'], join_cals['fcfasec'], join_cals['fcfasecerr'])

# Source masks
crl2688 = lrvalues['Source'] == 'CRL2688'
crl618 = lrvalues['Source'] == 'CRL618'
uranus = lrvalues['Source'] == 'URANUS'


# Filter out observations with OMP state='Bad', 'Questionable', 'Junk' or 'Rejected' # Actually this was done already.
#plt.ion()
uranuscolor = '#12eda8' #'#66c2a5'
crl2688color ='#ff7f4d'# '#fc8d62'
crl618color = '#7495e2'#'#8da0cb'

# Histogram of legacy FCFS (including by source)
fig = plt.figure(figsize=(3.394,2.5))
ax = fig.add_subplot(111)
ax.hist(lrvalues['FCFasec'], bins=50, range=(2,3), histtype='step', align='mid',
        color='0.7', facecolor='0.7', fill=True,
        label='All: {:.2f}$\pm${:.2f}'.format(lrclippedresults['All'][0],
                                                   lrclippedresults['All'][1]))

ax.hist(lrvalues['FCFasec'][uranus], bins=50, range=(2,3), histtype='step', align='mid',
        color=uranuscolor, edgecolor='black', hatch='', fill=True, facecolor=uranuscolor, alpha=0.5,
        label='Uranus: {:.2f}$\pm${:.2f}'.format(lrclippedresults['URANUS'][0],
                                                   lrclippedresults['URANUS'][1]))

ax.hist(lrvalues['FCFasec'][crl618], bins=50, range=(2,3), histtype='step', align='mid',
        color=crl618color, edgecolor='black',
        hatch='', fill=True, facecolor=crl618color, alpha=0.5,
        label='CRL 618: {:.2f}$\pm${:.2f}'.format(lrclippedresults['CRL618'][0],
                                                   lrclippedresults['CRL618'][1]))

ax.hist(lrvalues['FCFasec'][crl2688], bins=50, range=(2,3), histtype='step',
        align='mid', edgecolor='black',
        color=crl2688color, hatch='', fill=True, facecolor=crl2688color, alpha=0.5,
        label='CRL 2688: {:.2f}$\pm${:.2f}'.format(lrclippedresults['CRL2688'][0],
                                                   lrclippedresults['CRL2688'][1]))
ax.legend(frameon=False, loc=2)
ax.set_xlabel('arcsecond FCF')
ax.set_ylabel('Count')
fig.tight_layout()
xlim = ax.get_xlim()
ylim = ax.get_ylim()
ax.vlines(lrclippedresults['All'][0], ylim[0], ylim[1], linestyle='dashed', color='black')
#ax.fill_betweenx(ylim, lrclippedresults['All'][0] - lrclippedresults['All'][1],
#                 lrclippedresults['All'][0] + lrclippedresults['All'][1],
#                 color='0.5', alpha=0.3, zorder=-5)
ax.vlines(lrclippedresults['All'][0] + lrclippedresults['All'][1], ylim[0], ylim[1], linestyle='solid', color='0.5', alpha=0.3)
ax.vlines(lrclippedresults['All'][0] - lrclippedresults['All'][1], ylim[0], ylim[1], linestyle='solid', color='0.5', alpha=0.3)
fig.savefig('lrvalues-histo.pdf', bbox_inches='tight', pad_inches=0.0)

# Plot of Legacy FCFs versions transmission.
fig = plt.figure(figsize=(3.394,2.5))
ax = fig.add_subplot(111)

#cont1 = ax.errorbar(lrvalues['Trans'][crl2688], lrvalues['FCFasec'][crl2688], yerr=lrvalues['FCFasec_err'][crl2688], ecolor=crl2688color, mec=crl2688color,elinewidth=0.5, capsize=3, fmt='o', mfc='none',label='CRL2688', alpha=1.0)
#cont2 = ax.errorbar(lrvalues['Trans'][crl618], lrvalues['FCFasec'][crl618], yerr=lrvalues['FCFasec_err'][crl618], ecolor=crl618color, mec=crl618color, elinewidth=0.5, capsize=3, fmt='x', mfc='none', label='CRL618')
#cont3 = ax.errorbar(lrvalues['Trans'][uranus], lrvalues['FCFasec'][uranus], yerr=lrvalues['FCFasec_err'][uranus], ecolor=uranuscolor, mec=uranuscolor,elinewidth=0.5, capsize=3, fmt='^', mfc='none', label='Uranus')
ax.scatter(lrvalues['Trans'][uranus], lrvalues['FCFasec'][uranus], s=20, marker='^', edgecolors=uranuscolor, label='Uranus', facecolor='none')
ax.scatter(lrvalues['Trans'][crl2688], lrvalues['FCFasec'][crl2688], s=50, marker='+', facecolors=crl2688color, label='CRL 2688')
ax.scatter(lrvalues['Trans'][crl618], lrvalues['FCFasec'][crl618], s=20, marker='x', facecolors=crl618color, label='CRL 618')

#ax.set_title('Legacy Release arcsecond FCFs vs Transmission')
ax.legend(frameon=False, fontsize='small', title='LR-reduced Calibrators')
ax.set_xlabel('850$\mu$m Transmission')
ax.set_ylabel('arcsecond FCF')
ax.set_xlim(0, 1.06)
xlim= ax.get_xlim()
ax.set_ylim(1.56, 4.0)
ylim = ax.get_ylim()
clipped, low, high = scipy.stats.sigmaclip(lrvalues['FCFasec'], low=5.0, high=5.0)
clippedmean = np.mean(clipped)
clippedstd = np.std(clipped)
ax.fill_between(xlim, clippedmean-clippedstd, clippedmean+clippedstd, alpha=0.2, color='0.7')
ax.hlines(clippedmean, ax.get_xlim()[0], ax.get_xlim()[1], color='black', linestyle='dashed', zorder=1000)
ax.hlines(clippedmean+clippedstd, ax.get_xlim()[0], ax.get_xlim()[1], color='0.7', linestyle='dashed')
ax.hlines(clippedmean-clippedstd, ax.get_xlim()[0], ax.get_xlim()[1], color='0.7', linestyle='dashed')
ax.vlines(0.87, ylim[0]+0.25, 3.3, color='black', linestyle='dotted')
fig.tight_layout()
fig.savefig('legacyFCF-vs-transmission.pdf', bbox_inches='tight', pad_inches=0.01)




# Plot of Legacy FCFs vs time.
hstdates = np.array([datetime.datetime.strptime(i, '%Y-%m-%dT%H:%M:%S') for i in lrvalues['HST']])
hsttimes = [i.time() for i in hstdates]
hsttimes = np.array([i.hour + i.minute/60 for i in hsttimes])
fig = plt.figure(figsize=(7.1,2.0))
ax = fig.add_subplot(111)
ax.scatter( hstdates[crl2688], lrvalues['FCFasec'][crl2688], marker='s', color=crl2688color,facecolor='none', label='CRL 2688')
ax.scatter( hstdates[crl618], lrvalues['FCFasec'][crl618], marker='o', color=crl618color, facecolor='none',label='CRL 618')
ax.scatter( hstdates[uranus], lrvalues['FCFasec'][uranus],  marker='x', color=uranuscolor, label='Uranus ')
ax.set_ylim(1.5, 4.5)
ax.legend(frameon=False, title='LR FCFs by date', loc=2)
ax.set_ylabel('arcsecond FCF')
ax.set_xlabel('Date (HST)')
fig.tight_layout()
fig.savefig('legacyFCF-vs-date.pdf', bbox_inches='tight', pad_inches=0.01)

# fig = plt.figure()
# ax= fig.add_subplot(111)
# bins = None
# print('Uranus arcsecond FCFs from the caldb')
# for i in range(5):
#     if  bins is None:
#         bins = 20
#     fcfs = regularcals['fcfasec'][ruranus & (years==(2011+i))]
#     res = ax.hist(fcfs, bins=20, range=(2.0,3.0),histtype='step', align='mid', label='%i' % (2011+i), fill='solid', alpha=0.2)
#     print('{}: mean={:.2f}, median={:.2f}, std={:.2f}, count={}'.format(2011+i, np.mean(fcfs), np.median(fcfs), np.std(fcfs), len(fcfs)))
#     bins = res[1]

# ax.legend(frameon=False)
# ax.set_title('Histogram of Uranus arcsecond FCFs from Calibrator DB by year')


# Comparisons of CAL DB and Legacy release FCFS: histograms.
crl2688 = join_cals['Source'] == 'CRL2688'
crl618 = join_cals['Source'] == 'CRL618'
uranus = join_cals['Source'] == 'URANUS'

fig = plt.figure(figsize=(7.11,4))
ax = fig.add_subplot(121)
res = ax.hist(join_cals['FCFasec'], bins=50, range=(2,3), histtype='step', align='mid', color='0.7', facecolor='0.7', label='Legacy Release', fill= True)
res2 = ax.hist(join_cals['fcfasec'], bins=50, range=(2,3), histtype='step', align='mid', color='black', label='Std. Cal. Reduction', hatch='////')
ax.legend(frameon=False, title='All Calibrator sources')
ax.set_xlabel('arcsecond FCF')
ax.set_ylabel('Count')

ax = fig.add_subplot(322)
ax.hist(join_cals['FCFasec'][uranus], bins=50, range=(2,3), histtype='step', align='mid', color=uranuscolor, fill=True, facecolor=uranuscolor, label='LR reduction')
ax.hist(join_cals['fcfasec'][uranus], bins=50, range=(2,3), histtype='step', align='mid', color='black', hatch='////', label='Std. Cal. reduction')
ax.legend((),(), frameon=False, title='Uranus')
ax.tick_params(labelbottom='off')

ax = fig.add_subplot(324)
ax.hist(join_cals['FCFasec'][crl618], bins=50, range=(2,3), histtype='step', align='mid', color=crl618color, fill=True, facecolor=crl618color, label='LR reduction')
ax.hist(join_cals['fcfasec'][crl618], bins=50, range=(2,3), histtype='step', align='mid', color='black', hatch='////', label='Std. Cal. reduction')
ax.legend((),(), frameon=False, title='CRL 618')
ax.tick_params(labelbottom='off')

ax = fig.add_subplot(326)
ax.hist(join_cals['FCFasec'][crl2688], bins=50, range=(2,3), histtype='step', align='mid', color=crl2688color, fill=True, facecolor=crl2688color, label='LR reduction')
ax.hist(join_cals['fcfasec'][crl2688], bins=50, range=(2,3), histtype='step', align='mid', color='black', hatch='////', label='Std. Cal. reduction')
ax.legend((),(), frameon=False, title='CRL 2688')
ax.set_xlabel('arcsecond FCF')
fig.tight_layout()
fig.savefig('legacyFCF-caldbFCF-histograms.pdf', bbox_inches='tight', pad_inches=0.01)

# Comparisons of CAL DB and Legacy release FCFS: scatter plot.
fig = plt.figure(figsize=(3.394, 2.5))
ax = fig.add_subplot(111)
ax.scatter(join_cals['FCFasec'][uranus], join_cals['fcfasec'][uranus], marker='^', color=uranuscolor, label='Uranus', facecolor='none')
ax.scatter(join_cals['FCFasec'][crl2688], join_cals['fcfasec'][crl2688], s=50, marker='+', color=crl2688color, label='CRL 2688')
ax.scatter(join_cals['FCFasec'][crl618], join_cals['fcfasec'][crl618], s=20, marker='x', facecolors=crl618color, label='CRL 618')
ax.set_xlim(1.5,4.0)
ax.set_ylim(1.5, 4.0)
ax.set_xlabel('Legacy Reduction arcsecond FCF')
ax.set_ylabel('Std. Cal. reduction arcsecond FCF')
ax.plot([0,11], np.array([0,11]) - (2.48-2.34), linestyle='dashed', color='black', label='1:1 line')
ax.legend(frameon=False, title='Comparison of LR and Std. Cal. FCFs', loc=2)
fig.tight_layout()
fig.savefig('legacyFCF-caldbFCF-scatter.pdf', bbox_inches='tight', pad_inches=0.01)

