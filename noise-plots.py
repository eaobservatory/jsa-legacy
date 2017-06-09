
from astropy.io import fits
# from jsa_proc.config import get_database
# from jsa_proc.state import JSAProcState
# from jsa_proc.admin.directories import get_log_dir, get_output_dir, make_temp_scratch_dir, get_scratch_dir
import os
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
plt.ion()

#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{apjfonts},\usepackage{sfmath}'
matplotlib.rcParams['figure.autolayout'] = True
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['legend.fontsize'] = 'medium'
matplotlib.rcParams['axes.linewidth'] = '0.5'
matplotlib.rcParams['axes.formatter.limits'] = [-3,3]

db = get_database()
coadds = db.find_jobs(task='hpx-s2-850-r2-coadd', state=JSAProcState.PROCESSED, outputs='%.fits')



bins1 = np.linspace(1e-6, 500, num=1e5)
bins2 = np.logspace(np.log10(1e-6), np.log10(500), num=1e5)
summedcounts1 = np.zeros(len(bins1) -1)
summedcounts2 = np.zeros(len(bins2) -1)


meanvars = {}
for j in coadds:
    print(j.id, coadds.index(j))
    outputdir = get_output_dir(j.id)
    fitsfile = os.path.join(outputdir, j.outputs[0])
    var = fits.getdata(fitsfile, ext=1)[0,:,:]
    varflat = var[np.isfinite(var)]

    pixelcount = np.isfinite(var).sum()
    mean = np.mean(varflat)
    meanvars[j.id] = mean, pixelcount

    counts1, edges = np.histogram(varflat, bins=bins1)
    summedcounts1 += counts1
    counts2, edges = np.histogram(varflat, bins=bins2)
    summedcounts2 += counts2
    del(var)
    del(varflat)

rmsbins2 = np.sqrt(bins2)
rmsbins1 = np.sqrt(bins1)

rmsbins2centres = rmsbins2[:-1] + 0.5*(rmsbins2[1:]  - rmsbins2[:-1])
rmsbins1centres = rmsbins1[:-1] + 0.5*(rmsbins1[1:] - rmsbins1[:-1])
histtext = open('histogram_values_logbins', 'w')
for i,j in zip(rmsbins2centres, summedcounts2):
    histtext.write('{},{}\n'.format(i,j))
histtext.close()

histtext = open('histogram_values_linbins', 'w')
for i,j in zip(rmsbins1centres, summedcounts1):
    histtext.write('{},{}\n'.format(i,j))
histtext.close()

histtext = open('histogram_values_logbins','r')
vals = [i.strip() for i in histtext.readlines()]
rmsbins2centres = [float(i.split(',')[0]) for i in vals]
summedcounts2 = [float(i.split(',')[1]) for i in vals]

histtext = open('histogram_values_linbins','r')
vals1 = [i.strip() for i in histtext.readlines()]
rmsbins1centres = [float(i.split(',')[0]) for i in vals1]
summedcounts1 = [float(i.split(',')[1]) for i in vals1]



fig = plt.figure(figsize=(3.394,2.5))
ax = fig.add_subplot(111)
ax.plot( rmsbins2centres, summedcounts2, color='black')
ax.set_ylabel('Pixel count')
ax.set_yscale('symlog')
ax.set_xlabel(r'RMS noise on a pixel (mJy\,arcsec\textsuperscript{$-$2})')
ax.set_xscale('log')

ax2 = fig.add_axes([0.72, 0.67, 0.2, 0.2])
ax2.tick_params(labelsize=4, pad=2)


peakrms = rmsbins2centres[np.argmax(summedcounts2)]
print('Peak value found at {:.2F} RMS'.format(peakrms))

lin = ax.vlines(0.005, ax.get_ylim()[0], ax.get_ylim()[1], linestyle='dotted', alpha=0.5)
lin = ax.vlines(peakrms, ax.get_ylim()[0], ax.get_ylim()[1], linestyle='dashed', alpha=0.5)


ax2.plot(rmsbins2centres, summedcounts2, color='black')
ax2.set_xscale('log')
ax2.minorticks_off()
fig.tight_layout()
fig.savefig('coadds-noise-histogram.pdf', bbox_inches='tight', pad_inches=0.01, dpi=200)

