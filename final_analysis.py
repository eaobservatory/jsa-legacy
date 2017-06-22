from astropy.table import Table, Column, join
import astropy.units as u
import matplotlib.pyplot as plt
from scipy import stats

import datetime
import glob
import numpy as np
import matplotlib
import os

from starlink import kappa
import astropy.coordinates as coord
#plt.ion()
FCF_arcsec = 2.46
FCF_beam = 571

matplotlib.rcParams['font.size']=7
matplotlib.rcParams['lines.linewidth']=0.5
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['xtick.major.width']=0.5
matplotlib.rcParams['xtick.minor.width']=0.5
matplotlib.rcParams['ytick.major.width']=0.5
matplotlib.rcParams['ytick.minor.width']=0.5
matplotlib.rcParams['font.family'] = 'serif'
#oldserif = matplotlib.rcParams['font.serif']
#matplotlib.rcParams['font.serif'] = ['Times New Roman, Times']
matplotlib.rcParams['font.weight'] = 'medium'

cogtable = 'deep_coadds_flux_by_radius.csv'


fcftable = 'calanalysis-2017-customphotom-105-120-multiplesourceradii.csv'
fcfbkgrnd = (105,120) * u.arcsec
sourceradius = 60.0 * u.arcsec
asttable_filename = 'astmask_area_vs_fluxscale_simulation.csv'

simulatedfcfs = 'simfcfs.csv'
astmaskareatable = 'astmaskarea.csv'
fakemapfile = '/net/kapili/export/data/sgraves/LegacyCalibrator_2017/simulation/simulation2/fakesource850.sdf'

sourcecolordict = {'CRL618': 'C0',
                   'CRL2688': 'C1',
                   'URANUS': 'C2',
                   'Arp220': 'C3'}

pointingtable = '/home/sgraves/jsa-legacy/pointing-offsets.csv'
allskynoise = '/home/sgraves/jsa-legacy/average-rms-per-tile-lessthan200pixelsexcluded.csv'
# Calibration Plots:
# 1. curve of growth for deep coadds. (all three sources)
# 2. histograms of 850um FCFs (beam and arcsec.)
# 3. relative FCF for simulated bright source at various apertures &brightness
# 4. Area  included in AST mask at different brightnesses
# 5. Relative FCF vs area in AST mask (colored by brightness.)

# Maybe? Flux in background for various backgrounds at various brightnesses (simulations)

#1. cog cumulative and non cumulative (Total flux included at a given
#radius, relative to flux included in radius of 60"). Derived from
#deep coadds of cal sources.

cogres = Table.read(cogtable)


fig = plt.figure(figsize=(3,4))
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax)

# remove the bkground values we don't want to use.
cogres = cogres[cogres['type']=='source']

# Go through each source, create
coggroup = cogres.group_by('source')
for g in coggroup.groups:

    # cumulative cog
    normalisation = g['total'][g['radius']==sourceradius.value][0]
    ax.plot(g['radius'], g['total']/normalisation,
            label=g['source'][0], drawstyle='steps', color=sourcecolordict[g['source'][0]])

    # non cumulative cog.
    noncumulative = g['total'][1:] - g['total'][:-1]
    normindex = np.where(g['radius']==sourceradius.value)[0].item()
    # normalize relative to values in annulus with 60.0 as outer radii
    noncum_normalize = noncumulative[0]
    ax2.plot(g['radius'][1:], noncumulative /noncum_normalize,
             label=g['source'][0], drawstyle='steps', color=sourcecolordict[g['source'][0]])


# Indicate on the plots the 60" aperture and the 105-120 background
# aperture used in the FCF analysis (note no background was removed
# for this cog plot.)
ax.axvspan(xmin=fcfbkgrnd[0].value, xmax=fcfbkgrnd[1].value, facecolor='black', alpha=0.1,
             label='AP background')
ax.axvline(x=sourceradius.value, color='black', linestyle='dashed', zorder=0,
           linewidth=0.5, label='AP source')
ax.legend(frameon=False)

ax.tick_params(labelbottom=False)
ax.set_ylabel('Relative cumulative flux')

ax2.axvspan(xmin=fcfbkgrnd[0].value, xmax=fcfbkgrnd[1].value, facecolor='black', alpha=0.1,
             label='AP background')
ax2.axvline(x=sourceradius.value, color='black', linestyle='dashed', zorder=0,
           linewidth=0.5, label='AP Source')


ax2.set_xlabel('Radius (")')
ax2.set_ylabel('Relative average flux')
ax2.axhline(0.0, color='black', linewidth=0.5)
ax.set_xlim(0,300)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
fig.savefig('cog-relative60arcs.pdf',
            bbox_inches='tight', pad_inches=0.01)



# 2. histograms of 850um FCFs (beam and arcsec.)
fcfres = Table.read(fcftable)

# Extract only the source and bkground levels we are using here.
mask = (fcfres['bkgrnd_inner_radius_asec'] == fcfbkgrnd[0].value) & \
       (fcfres['bkgrnd_outer_radius_asec'] == fcfbkgrnd[1].value) & \
       (fcfres['source_radius_asec'] == sourceradius.value)

fcfres = fcfres[mask]

noturanusmask = fcfres['source']!='URANUS'
fcf_nou = fcfres[noturanusmask]

# Get the non uranus values.
clippeda, lowa, higha = stats.sigmaclip(fcf_nou['FCF_asec'], low=5.0, high=5.0)
print('FCF_asec (CRL sources): {:.2f} +/- {:.2f} (from {:.2F} to {:.2f})'.format(clippeda.mean(),
                                                                                 clippeda.std(), lowa, higha))
fcf_asec_final = clippeda.mean()
clippedb, lowb, highb = stats.sigmaclip(fcf_nou['FCF_beam'], low=5.0, high=5.0)
print('FCF_beam (CRL sources): {:.2f} +/- {:.2f} (from {:.2F} to {:.2f})'.format(clippedb.mean(),
                                                                                 clippedb.std(), lowb, highb))

# Go through each source separately.
fcfgroups = fcfres.group_by('source')


fig = plt.figure(figsize=(7,4))
fig2 = plt.figure(figsize=(3.5,3.5))
axweather = fig2.add_subplot(211)
axtime = fig2.add_subplot(212)

binsb = 50
binsa = 50
count=1
markeriter = iter(['x','o','^'])
for g in fcfgroups.groups:

    axbeam=fig.add_subplot(3,2,count)
    axasec=fig.add_subplot(3,2,count+1)

    fcfbeam = g['FCF_beam']
    fcfasec = g['FCF_asec']
    source = g['source'][0]
    trans850 = g['850trans']
    dateobs = g['date-obs']
    dateobs = np.array([datetime.datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%f') for i in dateobs])
    color = sourcecolordict[source]
    marker = markeriter.__next__()

    # Clipped FCF values.
    clippedb, lowb, highb = stats.sigmaclip(fcfbeam, low=5.0, high=5.0)
    meanb = clippedb.mean()
    medianb = np.median(clippedb)
    medianb_all = np.median(fcfbeam)
    stdb = clippedb.std()

    clippeda, lowa, higha = stats.sigmaclip(fcfasec, low=5.0, high=5.0)
    meana = clippeda.mean()
    mediana = np.median(clippeda)
    mediana_all = np.median(fcfasec)
    stda = clippeda.std()

    print('{}: {} observations.'.format(source, len(g)))
    print('{}: Aperture median FCF: {:.2F}'.format(source, mediana_all))
    print('{}: Aperture clipped mean fCF: {:.2f} +/- {:.2f} (from {:.2F} to {:.2f})'.format(source,
                                                                                            meana, stda,
                                                                                            lowa, higha))
    print('{}: Beam median FCF: {:.2F}'.format(source, medianb_all))
    print('{}: Beam clipped mean fCF: {:.2f} +/- {:.2f} (from {:.2F} to {:.2f})'.format(source,meanb,
                                                                                            stdb, lowb, highb))



    # Histograms of all 3.

    valsb, binsb, patchesb = axbeam.hist(fcfbeam, bins=binsb, label=source, range=(450.0,700), histtype='step',
                                         color=color)
    axbeam.hist(fcf_nou['FCF_beam'], bins=binsb, label='CRL618+CRL2688', histtype='step',
                facecolor='black', edgecolor='none', fill=True, alpha=0.2, zorder=0)
    valsa, binsa, patchesa = axasec.hist(fcfasec, bins=binsa, label=source, histtype='step',range=(1.9, 3.2),
                                         color=color)

    axasec.hist(fcf_nou['FCF_asec'], bins=binsa, label='CRL618+CRL2688', histtype='step',
                facecolor='black',edgecolor='none', fill=True, alpha=0.2, zorder=0)


    axbeam.axvline(meanb, color=color, alpha=1.0)
    axbeam.axvline(medianb_all, color=color, linestyle='dashed')

    axbeam.hlines(valsb.max(), xmin=meanb - stdb, xmax=meanb + stdb, color=color, alpha=1.0)
    axbeam.axvspan(xmin=meanb - stdb, xmax=meanb + stdb, facecolor=color, alpha=0.1)

    axasec.axvline(meana, color=color, alpha=1.0)

    axasec.axvline(mediana_all, color=color, linestyle='dashed')
    axasec.hlines(valsa.max(), xmin=meana - stda, xmax=meana + stda, color=color, alpha=1.0)
    axasec.axvspan(xmin=meana - stda, xmax=meana + stda, facecolor=color, alpha=0.1)
    axbeam.text(450, 100.0, source, ha='left')

    axasec.yaxis.tick_right()
    if count !=5:
        axasec.tick_params(labelbottom=False)
        axbeam.tick_params(labelbottom=False)
    count += 2

    # Plot the weather and time histograms for arcsec.
    coll=axweather.scatter(trans850, fcfasec,  marker='.', color=color, edgecolor='none',
                           label=source, alpha=0.5, s=12)
    #coll.set_edgecolor(sourcecolordict[source])
    coll=axtime.scatter(dateobs, fcfasec,  marker='.', edgecolor='none',
                        color=color, label=source, alpha=0.5, s=12)
    #coll.set_edgecolor(sourcecolordict[source])

[i.set_ha('right') for i in axtime.get_xticklabels()]
[i.set_rotation(30.0) for i in axtime.get_xticklabels()]
[i.set_rotation_mode('anchor') for i in axtime.get_xticklabels()]
axweather.set_ylim(1.8,4.0)
axtime.set_ylim(1.8,4.0)
axweather.legend(frameon=False)
# make legend entries not have alpha level.
[i.set_alpha(1.0) for i in axweather.legend_.legendHandles]
axweather.set_xlabel('850 $\mu$m transmission')
axtime.set_ylabel('FCF aperture')
axtime.set_xlabel('Date Obs')
fig2.subplots_adjust(hspace=0.5)
#Draw  a horizontal line at the FCF value.
axweather.axhline(fcf_asec_final, color='black', label=r'FCF$_{\mathrm{asec}}$' + '={:.2F}'.format(fcf_asec_final))
axweather.axhspan(ymin = fcf_asec_final - clippeda.std(), ymax = fcf_asec_final + clippeda.std(),
                  facecolor='black', alpha=0.2, zorder=0, edgecolor='none')
axtime.axhline(fcf_asec_final, color='black', label=r'FCF$_{\mathrm{asec}}$' + '={:.2F}'.format(fcf_asec_final))
axtime.axhspan(ymin = fcf_asec_final - clippeda.std(), ymax = fcf_asec_final + clippeda.std(),
               facecolor='black', alpha=0.2, zorder=0, edgecolor='none')
axtime.spines['top'].set_visible(False)
axweather.spines['top'].set_visible(False)
axtime.spines['right'].set_visible(False)
axweather.spines['right'].set_visible(False)
fig2.savefig('fcf-extra.pdf', bbox_inches='tight', pad_inches=0.01)
fig.axes[0].set_title('Beam FCF')
fig.axes[1].set_title('Aperture FCF')
fig.subplots_adjust(wspace=0.2)
fig.axes[4].set_xlabel(r'FCF$_{\mathrm{beam}}$ Jy pW$^{-1}$ beam$^{-1}$')
fig.axes[5].set_xlabel(r'FCF$_{\mathrm{aperture}}$ Jy pW$^{-1}$ arcsec$^{-2}$')
for i in fig.axes:
    i.spines['top'].set_visible(False)

fig.axes[0].spines['right'].set_visible(False)
fig.axes[2].spines['right'].set_visible(False)
fig.axes[4].spines['right'].set_visible(False)

fig.axes[1].spines['left'].set_visible(False)
fig.axes[3].spines['left'].set_visible(False)
fig.axes[5].spines['left'].set_visible(False)

fig.savefig('fcf-histogram.pdf', bbox_inches='tight', pad_inches=0.01)



# Simulation figure
# 3. relative FCF for simulated bright source at various apertures &brightness
simtable = Table.read(simulatedfcfs)

# Get the peak of the fakemap used as input.
fakepeak = kappa.stats(fakemapfile).maximum


# Only look at our chosen background radius (no real difference in
# qualitative effect, but easire to examine for only one.)
simtable = simtable[(simtable['bkgrnd_inner_radius_asec']==fcfbkgrnd[0].value) & \
                    (simtable['bkgrnd_outer_radius_asec']==fcfbkgrnd[1].value)]

# Get the values for the input fakesource and the simulations seperately.
simmeasures = simtable[simtable['source']!='URANUS']
inputmeasures = simtable[simtable['source']=='URANUS']

#Now plot for each souraceradius.
fig = plt.figure(figsize=(2.5,3))
ax = fig.add_subplot(111)


sourceradii = list(set(simmeasures['source_radius_asec']))
sourceradii.sort()
sourceradii = [i for i in sourceradii if i < 100.0]

# Remove the colors for CRLs & uranus; also add line style cycler.

coliter =  iter(['black', '#d62728', '#9467bd', '#8c564b', '#e377c2'])
lsiter  = iter(['-', '--', '-.', ':',(0, (6,2,2,4)) ])

# Get the ast mask size at various flux radii.
basepath = '/net/kapili/export/data/sgraves/LegacyCalibrator_2017/simulation/simulation2/'
pathtypes = ['simulatedbrightsource_peak0_1.sdf',
             'simulatedbrightsource_peak0_0[0-9].sdf',
             'simulatedbrightsource_peak0_0[0-9][0-9].sdf',
             ]

files = []
for p in pathtypes:
    files+=glob.glob(os.path.join(basepath, p))

#Get ast mask in all files.
# output = []
# for filename in files:
#     peakvalue = float(os.path.splitext(filename)[0].split('peak')[1].replace('_', '.'))

#     # Get the quality count: have to do this from stdout, as showqual
#     # doesn't write the count into the ADAM file.
#     res, stdout = kappa.showqual(filename, count=True, returnStdOut=True)
#     not_in_astmask_count =  int([i.strip() for i in stdout.split('\n')
#                      if i.strip().startswith('AST')][0].split()[-1].strip('(').strip(')'))
#     not_in_fltmask_count = int([i.strip() for i in stdout.split('\n')
#                      if i.strip().startswith('FLT')][0].split()[-1].strip('(').strip(')'))
#     total = kappa.stats(filename).numpix
#     output += [[peakvalue, total - not_in_astmask_count, total - not_in_fltmask_count]]



# asttable = Table(np.array(output), names=('fluxscale', 'astcount', 'fltcount'))
# asttable.sort('fluxscale')
# asttable.write(asttable_filename)

asttable = Table.read(asttable_filename)
fig3 = plt.figure(figsize=(3,2.5))
ax3 = fig3.add_subplot(111)
ax3.plot(asttable['fluxscale'], 3.22**2 * asttable['astcount'], marker='x', color='black')
ax3.set_xlabel('Simulated peak brightness (pW)')
ax3.set_ylabel('Area in AST mask (arcsec$^{2}$)')
fig3.savefig('astmask-peakbrightness.pdf', bbox_inches='tight', pad_inches=0.01)


for sr in sourceradii:
    # Only look at that source radius.
    fakemap = inputmeasures[inputmeasures['source_radius_asec']==sr]
    simulation = simmeasures[simmeasures['source_radius_asec']==sr]

    # Get the expected scale of each simulation from the filename,
    # then the expected multiplicative factor from the fakemap.
    scales = np.array([float(i.split('peak')[1].replace('_', '.')) for i in simulation['filename']])
    factors = scales/fakepeak

    #Calculate the 'expected' FCF, based on the simulation multiplied by the factors.
    expected_FCF = factors * (fakemap['source_fluxd_sum'] -
                              (fakemap['bkgrnd_mean'] * fakemap['source_numgood']))
    simulated_FCF = (simulation['source_fluxd_sum'] -
                                   (simulation['bkgrnd_mean'] * simulation['source_numgood']))
    relative_FCF = expected_FCF/simulated_FCF

    # Create a table and sort by the fluxscale (for plotting)
    fcftemptable = Table([Column(name='fluxscale', data=scales),
                          Column(name='relative_FCF_asec', data=relative_FCF.data),
                          Column(name='simulated_FCF_asec', data=simulated_FCF.data),
                          Column(name='simulated_flux', data=simulation['source_fluxd_sum']),
                          Column(name='factors', data=factors)])
    fcftemptable.sort('fluxscale')

    # Now plot these values.
    color = coliter.__next__()
    ls = lsiter.__next__()
    ax.plot(fcftemptable['fluxscale'].data, fcftemptable['relative_FCF_asec'].data,label=str(sr),
            color='black',ls=ls)

    # Now add a plot of area in ast mask vs realtive FCF.
    #mytab = join(fcftemptable, asttable, keys=['fluxscale'])
    #mytab.sort('fluxscale')
    #ax2.plot(mytab['relative_FCF_asec'], 1/(np.sqrt((3.22**2) * mytab['astcount']/np.pi)),
    #         label=str(sr),
    #         color=color, ls=ls)
    #ax3.plot(mytab['relative_FCF_asec'].data, 1/(np.sqrt((3.22**2) * mytab['fltcount']/np.pi)), label=str(sr),
    #         color=color, ls=ls)
#ax3.set_title('expected_FCF_asec')
#ax2.set_title('simualted_FCF_asec')


leg1=ax.legend(title='Aperture radius (")', frameon=False, loc=9)
ax.set_ylabel('Relative arcsec FCF')
ax.set_xlabel('Simulated source brightness (pW)')

# Label mean peak pW of CRL618, CRL2688 and Uranus
lines=[]
for src in ['CRL618', 'CRL2688', 'URANUS']:
    peak = np.mean(fcfres['amplitude'][fcfres['source']==src])
    line = ax.axvline(peak, color=sourcecolordict[src], linestyle='solid', label=src)
    lines.append(line)

leg2 = ax.legend(handles=lines, loc=10, title='Average source peaks', frameon=False)
ax.add_artist(leg1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig('simulated-fcfs.pdf', bbox_inches='tight', pad_inches=0.01)



# Pointing plot
pointtab = Table.read(pointingtable)
pointtab.add_column(Column(name='radial_offset', data=np.sqrt(pointtab['offset1']**2 + pointtab['offset2']**2)))
dateobs = np.array([datetime.datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in pointtab['date-obs']])
dateend = np.array([datetime.datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in pointtab['date-end']])
pointtab.replace_column('date-end',dateend)
pointtab.replace_column('date-obs',dateobs)


# Do all pointings and calibrators.

fig = plt.figure(figsize=(2.5, 2.5))
ax = fig.add_subplot(111)
ax.hist(pointtab['radial_offset'], bins=30, range=(0,10), color='black', label=None, histtype='step')
median = np.median(pointtab['radial_offset'])
print('median radial pointing offset is: {}'.format(median))
ax.axvline(median, color='black', linestyle='solid', label='Median: {:.2F}"'.format(median))

# for s in ['CRL2688', 'CRL618', 'Arp220']:

#     values = calibs[calibs['source']==s]
#     offsets = np.sqrt(values['offset1']**2 + values['offset2']**2)
#     color = sourcecolordict[s]
#     facecolor = matplotlib.colors.ColorConverter.to_rgba(color, alpha=0.1)
#     n, bins, patches = ax.hist(offsets, bins=30, range=(0,7), label=s,
#                                edgecolor=color, facecolor=facecolor,
#                                fill=True,
#                                histtype='step')

#ax.legend()
ax.set_xlabel('Radial offset (")')
ax.set_ylabel('Count')
ax.set_xlim(0,10)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.savefig('pointing-offset.pdf', bbox_inches='tight', pad_inches=0.01)


#all sky noise map.
import pymoc
import healpy

averagerms = Table.read(allskynoise)

# Get a dictionary with HPX tilenumber as key.
hpxdict = dict(zip(averagerms['Tilenum'], averagerms['RMS']))

# Turn it into a MOC
mymoc = pymoc.MOC(order=6, cells=hpxdict.keys())

# Turn it into a np array of the right setup for healpy.
order = mymoc.order
skymap = np.zeros(12 * 4 ** order)
for cell in mymoc.flattened(order):
    skymap[cell] = hpxdict[cell]


maskedmap = healpy.ma(skymap, badval=0)
fig = plt.figure(figsize=(7,5))

cmap = matplotlib.cm.gist_heat
cmap.set_bad('0.7')
cmap.set_under('white')
healpy.visufunc.mollview(maskedmap, unit=r'mJy arcsec$^{-2}$',
                         nest=True, norm='log', cmap=cmap,
                         min=0.024, max=5, coord=('C'),
                         notext=True, cbar=True, title='', fig=fig.number)
ax = fig.axes[1]
text = ax.texts[0]
text.set_size('medium')
ax=fig.axes[0]
[i.set_linewidth(0.5) for i in ax.lines]

healpy.visufunc.graticule(coord='G')
fig.savefig('mollweide-average-noise-galacticaxes.pdf', bbox_inches='tight', pad_inches=0.01, dpi=200)
healpy.visufunc.delgraticules()
healpy.visufunc.graticule(coord='C')
fig.savefig('mollweide-average-noise-equataxes.pdf', bbox_inches='tight', pad_inches=0.01, dpi=200)
import subprocess
subprocess.call(['pdfcrop', 'mollweide-average-noise-galacticaxes.pdf', 'mollweide-average-noise-galacticaxes-crop.pdf'])
subprocess.call(['pdfcrop', 'mollweide-average-noise-equataxes.pdf', 'mollweide-average-noise-equataxes-crop.pdf'])
# Note figure will require cropping.



# Plots of noise across map.
noisehistogramlinear = '/home/sgraves/jsa-legacy/histogram_values_linbins'
noisehistogramlog = '/home/sgraves/jsa-legacy/histogram_values_logbins'

histtext = open(noisehistogramlog, 'r')
vals = [i.strip() for i in histtext.readlines()]
rmsbins2centres = [float(i.split(',')[0]) for i in vals]
summedcounts2 = [float(i.split(',')[1]) for i in vals]

histtext = open(noisehistogramlinear,'r')
vals1 = [i.strip() for i in histtext.readlines()]
rmsbins1centres = [float(i.split(',')[0]) for i in vals1]
summedcounts1 = [float(i.split(',')[1]) for i in vals1]

fig = plt.figure(figsize=(3.394,2.5))
ax = fig.add_subplot(111)
ax.plot( rmsbins2centres, summedcounts2, color='black')
ax.set_ylabel('Pixel count')
ax.set_yscale('symlog')
ax.set_xlabel(r'RMS noise on a pixel (mJy$\,$arcsec$^{-2}$)')
ax.set_xscale('log')

ax2 = fig.add_axes([0.72, 0.7, 0.2, 0.2])
ax2.tick_params(labelsize=6, pad=2)


peakrms = rmsbins2centres[np.argmax(summedcounts2)]
print('Peak RMS value found at {:.2F}'.format(peakrms))

lin = ax.axvline(0.005, linestyle='dotted', alpha=0.5, linewidth=0.7, color='black')
lin = ax.axvline(peakrms, linestyle='dashed', alpha=0.5, linewidth=0.7, color='black')


ax2.plot(rmsbins2centres, summedcounts2, color='black')
ax2.set_xscale('log')
ax2.minorticks_off()
ax2.set_xticks([1e-2, 1e0])
fig.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
fig.savefig('coadds-noise-histogram.pdf', bbox_inches='tight', pad_inches=0.01, dpi=200)

# Images of maps.
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.wcs import WCS
from wcsaxes import WCSAxes, WCSAxesSubplot
from astropy.visualization.mpl_normalize import ImageNormalize

def formataxes(ax):
    ax.grid(color='white', alpha=0.1, linestyle='solid', linewidth=1.0)
    ax.coords.frame.set_linewidth(0.5)
    racoord=ax.coords[0]
    racoord.set_axislabel('RA -- HPX')
    racoord.set_major_formatter('hh:mm:ss')
    racoord.set_ticks_position('bt')
    racoord.set_ticks(size=8)
    racoord.set_ticks(spacing=10*(1/60.)*u.deg)
    racoord.set_minor_frequency(4)
    #racoord.set_separator((r'\textsuperscript{h}',r'\textsuperscript{m}',r'\textsuperscript{s}'))
    racoord.set_separator((r'$^{h}$', '$^{m}$', '$^{s}$'))
    racoord.display_minor_ticks(True)

    deccoord=ax.coords[1]
    deccoord.set_axislabel('Dec -- HPX', minpad=0)
    deccoord.set_major_formatter('dd:mm')
    deccoord.set_ticks_position('lr')
    deccoord.set_ticks(size=8)
    deccoord.set_ticks(spacing=10*(1/60.)*u.deg)
    deccoord.set_minor_frequency(4)
    deccoord.display_minor_ticks(True)


def add_colorbar(ax):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0, axes_class=matplotlib.axes.Axes)
    cbar = plt.colorbar(mappable=ax.images[0], cax=cax, orientation='horizontal')
    cbar.set_label(r'mJy arcsec$^{-2}$')
    cax.xaxis.tick_top()
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation=90)
    cbar.ax.xaxis.set_label_position('top')

s=matplotlib._cm.cubehelix(s=0, r=-0.5, gamma=0.75)
matplotlib.cm.register_cmap(name='test', data=s, lut=128)

cmap = matplotlib.cm.get_cmap('test')
cmap.set_bad('0.7')

cmap_variance = matplotlib.cm.get_cmap('gist_heat_r')
cmap_variance.set_bad('0.7')
# g34plot
#jcmts850um_healpix030318_pub_000.fits
# Plot tile 30318 (G34.3)
from astropy.io import fits
import astropy

g34filename = 'jcmts850um_healpix030318_pub_000.fits'
hdu = fits.open(g34filename)[0]
wcs = WCS(hdu.header)
wcs = wcs.dropaxis(2)
data = hdu.data[0,:,:]

fig = plt.figure(figsize=(6.6,11))
ax = WCSAxesSubplot(fig, 1,2,1, wcs=wcs)
fig.add_axes(ax)

# color image.
norm = ImageNormalize(vmin=-0.1, vmax=20, stretch=astropy.visualization.SqrtStretch())

ax.imshow(data, cmap='test', origin='lower', norm=norm, interpolation='none')

formataxes(ax)

# Label
#ax.text(0.95,0.9, r'Tile 30318 (G34.27+0.15)' + '\n' + r'850$\mu$m emission', transform=ax.transAxes, size='medium', ha='right',va='top')
# Colorbar
add_colorbar(ax)

# Second image: noise array.
ax2 = WCSAxesSubplot(fig, 1,2,2, wcs=wcs)
fig.add_axes(ax2)


# Get variance data.
hdu_var = fits.open(g34filename)[1]
wcs = WCS(hdu_var.header).dropaxis(2)


norm = ImageNormalize(vmin=np.sqrt(0.0004), vmax=0.2, stretch=astropy.visualization.PowerStretch(a=0.46*2))
cmap_variance = matplotlib.cm.get_cmap('gist_heat_r')
cmap_variance.set_bad('0.7')
ax2.imshow(np.sqrt(hdu_var.data[0,:,:]), cmap=cmap_variance, origin='lower', norm=norm, interpolation='none')
formataxes(ax2)
add_colorbar(ax2)


#ax2.text(0.95,0.9, r'Tile 30318 (G34.27+0.15)' + '\n' + r'850$\mu$m RMS noise',
#         transform=ax2.transAxes, size='medium', ha='right',va='top')



ax2.coords[1].set_ticklabel_position('r')
ax2.coords[1].set_axislabel_position('')

fig.tight_layout()
#fig.subplots_adjust(wspace=0.035)

plt.draw()
fig.savefig('tile30318-g34-coadd-noise.pdf', bbox_inches='tight', pad_inches=0.01)

# Second figure: extents and peaks.
from starlink import wrapper, kappa, convert
maskedfile = '30318_extent_contourmask.fits'
if os.path.isfile(maskedfile):
    os.remove(maskedfile)
kappa.ndfcopy('jcmts850um_extent-mask030318_pub_000.fits(,)', out='temp.sdf')
kappa.nomagic('temp.sdf', repval=0, out=maskedfile)


fig = plt.figure(figsize=(6.6,11))
ax =  WCSAxesSubplot(fig, 1,2,1, wcs=wcs)
fig.add_axes(ax)
norm = ImageNormalize(vmin=-0.1, vmax=20, stretch=astropy.visualization.SqrtStretch())


hdu_extents = fits.open(maskedfile)[0]
wcs_contour = WCS(hdu_extents.header)

ax.imshow(data, cmap='gist_yarg', origin='lower', norm=norm)
formataxes(ax)
#racoord=ax.coords[0]
#racoord.set_axislabel('RA -- HPX')
#racoord.set_major_formatter('hh:mm:ss')
#racoord.set_ticks_position('bt')
#racoord.set_ticklabel_position('bt')
#racoord.set_ticks(size=8)
#racoord.set_ticks(spacing=5*(1/60.)*u.deg)
#racoord.set_minor_frequency(6)
#racoord.set_separator((r'\textsuperscript{h}',r'\textsuperscript{m}',r'\textsuperscript{s}'))
#racoord.display_minor_ticks(True)

#deccoord=ax.coords[1]
#deccoord.set_axislabel('Dec -- HPX', minpad=0)
#deccoord.set_major_formatter('dd:mm')
#deccoord.set_ticks_position('lr')
#deccoord.set_axislabel_position('l')
#deccoord.set_ticklabel_position('l')
#deccoord.set_ticks(size=8)
#deccoord.set_ticks(spacing=5*(1/60.)*u.deg)
#deccoord.set_minor_frequency(6)
#deccoord.display_minor_ticks(True)
ax.grid(color='black', alpha=0.1, linestyle='solid', linewidth=1.0)

ax.contour(hdu_extents.data, transform=ax.get_transform(wcs_contour),
           colors='red', levels=[0.5], linewidths=0.7, label='Extents')
ax.set_xlim(37.0, 618.468)
ax.set_ylim(106.41, 634.41)


# Peaks.
ax2 =  WCSAxesSubplot(fig, 1,2,2, wcs=wcs)
fig.add_axes(ax2)
norm = ImageNormalize(vmin=-0.1, vmax=20, stretch=astropy.visualization.SqrtStretch())


hdu_extents = fits.open(maskedfile)[0]
wcs_contour = WCS(hdu_extents.header)
ax2.imshow(data, cmap='gist_yarg', origin='lower', norm=norm)
formataxes(ax2)
# racoord=ax2.coords[0]
# racoord.set_axislabel('RA -- HPX')
# racoord.set_major_formatter('hh:mm:ss')
# racoord.set_ticks_position('bt')
# racoord.set_ticklabel_position('bt')
# racoord.set_ticks(size=8)
# racoord.set_ticks(spacing=5*(1/60.)*u.deg)
# racoord.set_minor_frequency(6)
# racoord.set_separator((r'\textsuperscript{h}',r'\textsuperscript{m}',r'\textsuperscript{s}'))
# racoord.display_minor_ticks(True)

# deccoord=ax2.coords[1]
# deccoord.set_axislabel('Dec -- HPX', minpad=0)
# deccoord.set_major_formatter('dd:mm')
# deccoord.set_ticks_position('lr')
# deccoord.set_axislabel_position('r')
# deccoord.set_ticklabel_position('r')
# deccoord.set_ticks(size=8)
# deccoord.set_ticks(spacing=5*(1/60.)*u.deg)
# deccoord.set_minor_frequency(6)
# deccoord.display_minor_ticks(True)
# ax2.grid(color='black', alpha=0.1, linestyle='solid', linewidth=1.0)

ax2.contour(hdu_extents.data, transform=ax2.get_transform(wcs_contour),alpha=0.5,
            colors='black', linstyle='dashed',levels=[0.5], linewidths=0.5, label='Extents')
ax2.set_xlim(37.0, 618.468)
ax2.set_ylim(106.41, 634.41)
peakcatfile = 'jcmts850um_peak-cat030318_pub_000.fits'

# # New figure, showing extents and peaks.
# matplotlib.rcParams['figure.autolayout'] = False
# # contours
# contourfile = '30318_extent_contourmask_2d.fits'
# hdu = fits.open(contourfile)[0]
# wcs_contour = WCS(hdu.header)
# ax2.set_autoscale_on(False)
# ax2.contour(hdu.data, transform=ax2.get_transform(wcs_contour), colors='red', levels=[0.5], linewidths=1.0, label='Extents')

# plt.draw()
fig.subplots_adjust(wspace=0.025)

# #peaks
# peakcat = 'jcmts850um_peak-cat030318_pub_000.fits'


peakcat = Table.read(peakcatfile)
ra = peakcat['RA']
dec = peakcat['DEC']

ax2.scatter(ra, dec, s=45, transform=ax2.get_transform('fk5'), color='red', marker='x',
            label='Peaks', linewidths=0.7)
ax2.collections[0].set_label('Extents')
ax2.collections[1].set_linewidths(0.7)
#ax2.legend(frameon=False, fontsize='medium')
ax2.grid(color='black', alpha=0.2, linestyle='solid', linewidth=1.0)
ax2.coords[1].set_axislabel_position('r')
ax2.coords[1].set_ticklabel_position('r')
fig.tight_layout()
fig.savefig('tile30318-g34-extent-peak.pdf', bbox_inches='tight', pad_inches=0.05)


#CRL618: job 194622, 194624; extents 198838, 198840
t1238 = 'jcmts850um_healpix001238_pub_000.fits'
t1244 = 'jcmts850um_healpix001244_pub_000.fits'

t1238_extent = 'jcmts850um_extent-mask001238_pub_000.fits'
t1244_extent = 'jcmts850um_extent-mask001244_pub_000.fits'

t1238_peakcat = 'jcmts850um_peak-cat001238_pub_000.fits'
t1244_peakcat = 'jcmts850um_peak-cat001244_pub_000.fits'

pasted = 'crl618_pasted.fits'
pastedextents = 'crl618_extentmask_pasted.fits'
kappa.paste(in_=t1238, p1=t1244, out=pasted)
kappa.paste(in_=t1238_extent, p1=t1244_extent, out=pastedextents)



hdu = fits.open(pasted)[0]
wcs = WCS(hdu.header)
wcs2d = wcs.dropaxis(2)
fig = plt.figure(figsize=(6.5,11))
ax = WCSAxesSubplot(fig, 1,2,1, wcs=wcs2d)
fig.add_axes(ax)
data = hdu.data[0,:,:]
norm = ImageNormalize(vmin=-0.1, vmax=4, stretch=astropy.visualization.SqrtStretch())
ax.imshow(data, cmap=cmap, origin='lower', interpolation='none', norm=norm)
formataxes(ax)
add_colorbar(ax)

hduvar = fits.open(pasted)[1]
wcsvar = WCS(hduvar.header)
wcsvar2d = wcsvar.dropaxis(2)
ax2 = WCSAxesSubplot(fig, 1,2,2, wcs=wcs2d)
fig.add_axes(ax2)
datavar = np.sqrt(hduvar.data[0,:,:])
norm = ImageNormalize(vmin=3.5e-3, vmax=1, stretch=astropy.visualization.PowerStretch(0.3))
ax2.imshow(datavar, cmap=cmap_variance, origin='lower', interpolation='none', norm=norm)
conts = ax2.contour(datavar, colors='black', levels=[5e-3,6e-3,7e-3], linewidths=1.0)
formataxes(ax2)
add_colorbar(ax2)

ax2.coords[1].set_axislabel_position('r')
ax2.coords[1].set_ticklabel_position('r')
# #fig.subplots_adjust(wspace=0.09)
fig.subplots_adjust(wspace=0.035)
# plt.draw()
# #fig.set_tight_layout(True)
fig.tight_layout()
fig.savefig('crl618-whole-map.pdf', bbox_inches='tight', pad_inches=0.05)

#Zoomed in figure
fig = plt.figure(figsize=(6.6,2.3))
ax = WCSAxesSubplot(fig, 1,1,1, wcs=wcs2d)
fig.add_axes(ax)
norm = ImageNormalize(vmin=-0.03, vmax=0.3, stretch=astropy.visualization.PowerStretch(0.7))
ax.imshow(data, cmap=cmap, origin='lower', interpolation='none', norm=norm)
ax.set_xlim(212.53703020598158, 356.00391290220477)
ax.set_ylim(166.26655865188837, 287.7506125478838)

# Extent masks
hduextent = fits.open(pastedextents)[0]
wcs_extent = WCS(hduextent.header)
wcs_extent2d = wcs_extent.dropaxis(2)
dataextent = hduextent.data[0,:,:]
dataextent[np.isfinite(dataextent)] = 1.0
dataextent[np.isnan(dataextent)] = 0.0
ax.contour(dataextent, transform=ax.get_transform(wcs_extent2d), levels=[0.5], colors='white')
ax.coords[0].set_major_formatter('hh:mm:ss')
ax.coords[0].set_separator((r'$^{h}$',r'$^{m}$',r'$^{s}$'))
ax.coords[1].set_axislabel('Dec -- HPX', minpad=-0.5)
ax.coords[0].set_axislabel('RA -- HPX')
ax.coords[0].set_ticks(size=6, color='white')
ax.coords[1].set_ticks(size=6, color='white')
ax.coords[0].set_ticks_position('bt')
ax.coords[1].set_ticks_position('lr')


divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="180%", pad=0.1, axes_class=matplotlib.axes.Axes)
xvals = np.arange(201) -100
cax.plot(xvals, data[229-100:229+101, 286] - 0.02, color='black', drawstyle='default')
cax.plot(xvals, data[229-100:229+101, 286] - 0.02, color='cyan', drawstyle='default', linestyle='dashed')
cax.plot(xvals, data[229,286-100:286+101]+0.02, color='black', drawstyle='default', linestyle='solid')
cax.plot(xvals, data[229,286-100:286+101]+0.02, color='red', drawstyle='default', linestyle='dotted')
cax.set_ylim(-0.06, 0.06)
ax.vlines(286, 229-100, 229+100, color='cyan', linestyle='dashed', linewidth=1.0)
ax.hlines(229, 286-100, 286+100, color='Tomato',linestyle='dotted', linewidth=1.5)
cax.tick_params(labelright=True, labelleft=False)
cax.axhline(-0.02, color='0.7')
cax.axhline(+0.02, color='0.7')

cax.set_ylabel(r'mJy arcsec$^{-2}$')
cax.yaxis.set_label_position('right')
cax.tick_params(right=True, left=True)
cax.set_xlabel('Pixels (offset from peak)')
fig.tight_layout()
fig.savefig('crl618-sourceonly.pdf', bbox_inches='tight', pad_inches=0.05)



# Deep cosmos region.
# job 197594: CATALOGUE 201810
mapfile = 'jcmts850um_healpix027258_pub_000.fits'

# Catalogs
# Casey et al paper: UH observations. J/MNRAS/436/1919/table2
uhcat = 'vizier_votable.vot'

# CLS DR1 catalogs:
clscat = 'S2CLS_CATALOGUE_DR1.FITS'

# Our peak catalog
lrcat = 'jcmts850um_peak-cat027258_pub_000.fits'

hdu = fits.open(mapfile)[0]
data = hdu.data[0,:,:]
hduvar = fits.open(mapfile)[1]
wcs = WCS(hdu.header).dropaxis(2)
wcsvar = WCS(hduvar.header).dropaxis(2)
var = np.sqrt(hduvar.data[0,:,:])


fig = plt.figure(figsize=(6.5,7))
ax = WCSAxesSubplot(fig, 1,2,1, wcs=wcs)
fig.add_axes(ax)


norm = ImageNormalize(vmin=-0.01, vmax=0.0596, stretch=astropy.visualization.LinearStretch())
ax.imshow(data, cmap=cmap, origin='lower', interpolation='none', norm=norm)
ax.contour(var, levels=[1*FCF_arcsec/FCF_beam], colors='white')
add_colorbar(ax)
formataxes(ax)
ra=ax.coords[0]
ra.set_ticks(color='white', spacing=0.25*u.deg)
ra.set_minor_frequency(5)
dec = ax.coords[1]
dec.set_ticks(color='white', spacing=0.25*u.deg)
dec.set_minor_frequency(5)


ax2 = WCSAxesSubplot(fig, 1,2,2, wcs=wcsvar)
fig.add_axes(ax2)
norm = ImageNormalize(vmin=0.007, vmax=0.017, stretch=astropy.visualization.LinearStretch())
ax2.imshow(var, cmap=cmap_variance, origin='lower', interpolation='none', norm=norm)
ax2.contour(var, levels=np.array([0.25, 0.3, 0.5, 1, 1.5, 2])*FCF_arcsec/FCF_beam, colors='0.7')

add_colorbar(ax2)
formataxes(ax2)
ra=ax2.coords[0]
ra.set_ticks(color='0.7', spacing=0.25*u.deg)
ra.set_minor_frequency(5)
dec = ax2.coords[1]
dec.set_ticks(color='0.7', spacing=0.25*u.deg)
dec.set_minor_frequency(5)
add_colorbar(ax2)
dec.set_axislabel_position('r')
dec.set_ticklabel_position('r')
fig.tight_layout()
fig.savefig('27258-whole-map.pdf', bbox_inches='tight', pad_inches=0.05)


fig = plt.figure(figsize=(6.5,7))
ax = WCSAxesSubplot(fig, 1,2,1, wcs=wcs)
fig.add_axes(ax)
norm = ImageNormalize(vmin=-0.0097, vmax=0.029, stretch=astropy.visualization.LinearStretch())
ax.imshow(data, cmap=cmap, origin='lower', interpolation='none', norm=norm)
conts=ax.contour(var, levels=[1 * FCF_arcsec/FCF_beam], colors='white')


ax.set_xlim(117.38899509591401, 433.24228632982658)
ax.set_ylim(70.910809319165651, 384.07598743619394)

ax2 = WCSAxesSubplot(fig, 1,2,2, wcs=wcs)
fig.add_axes(ax2)
norm = ImageNormalize(vmin=-0.0097, vmax=0.029, stretch=astropy.visualization.LinearStretch())
ax2.imshow(data, cmap='gist_gray_r', origin='lower', interpolation='none', norm=norm)
conts=ax2.contour(var, levels=[1*FCF_arcsec/FCF_beam], colors='black')
ax2.set_xlim(117.38899509591401, 433.24228632982658)
ax2.set_ylim(70.910809319165651, 384.07598743619394)




uhtab = Table.read(uhcat)

# Get decimal degree ra and dec
dec=uhtab['DEJ2000']
dec=[i.decode() for i in dec]
dec=coord.Angle(dec, unit=u.degree)

ra=uhtab['RAJ2000']
ra=[i.decode() for i in ra]
ra=coord.Angle(ra, unit=u.hour)
ra=ra.to(unit=u.degree)
coll = ax2.scatter(ra.value, dec.value, s=45, transform=ax2.get_transform('fk5'),
                   label='UH',marker='o', color='black')
coll.set_edgecolor(coll.get_facecolor())
coll.set_facecolor('none')


clstab=Table.read(clscat)
ra=clstab['RA_DEG']
dec=clstab['Dec_DEG']
coll=ax2.scatter(ra, dec, s=45, transform=ax2.get_transform('fk5'), label='LCS',marker='^', zorder=6, color='purple')
coll.set_edgecolor(coll.get_facecolor())
coll.set_facecolor('none')


lrtab = Table.read(lrcat)
coll=ax2.scatter(lrtab['RA'], lrtab['DEC'], s=45, transform=ax2.get_transform('fk5'), marker='x', zorder=7, label='LR', color='red')
coll.set_linewidths(1.0)
fig.subplots_adjust(wspace=0.1)

formataxes(ax)
formataxes(ax2)

add_colorbar(ax)
add_colorbar(ax2)
fig.axes[-1].remove()
ax2.coords[1].set_axislabel('')
ax2.coords[1].set_ticklabel_position('r')
ax2.coords[1].set_axislabel_position('r')

fig.subplots_adjust(wspace=0.05)
fig.tight_layout()
fig.savefig('27258-zoomin.pdf', bbox_inches='tight', pad_inches=0.05)


# Beam analysis.
# Fit the beam to uranus, CRL618 and CRL2688. For uranus only, fit two gaussians (primary and secondary).

# Uranus map.
uranusdeep = 'uranus_allfiles_nearest.sdf'
crl618deep =  'crl618_all_rough.sdf'
crl2688deep = 'crl2688_all_rough.sdf'


beamfit_uranus = kappa.beamfit(ndf=uranusdeep, beams=2, mode='interface', polar=False, circular=True,
                               pos='"0,0"', pos2='"0,0"', resid = 'uranus_beamfit_residuals', fixback=0)

primaryamp = beamfit_uranus.amp[0]
secondaryamp = beamfit_uranus.amp[2]

crl618pos = '"4:42:53.5,36:06:54"'
crl2688pos = '"21:02:18.3, 36:41:37"'

beamfit_crl618 = kappa.beamfit(ndf=crl618deep, beams=1, mode='interface', polar=False, circular=True,
                               pos=crl618pos, resid = 'crl618_beamfit_residuals', fixback=0.)
beamfit_crl2688 = kappa.beamfit(ndf=crl2688deep, beams=1, mode='interface', polar=False, circular=True,
                                pos=crl2688pos, resid='crl2688_beamfit_residuals', fixback=0)
# Relative amplitudes
print('Uranus primary beam:')
print('    Primary  FWHM: {:.3f}"'.format((beamfit_uranus.majfwhm[0]*u.radian).to(u.arcsec)))
print('    Secon.   FWHM: {:.3f}"'.format((beamfit_uranus.majfwhm[2]*u.radian).to(u.arcsec)))
print('    Pri Amp (rel): {:.3f}'.format(primaryamp/ (primaryamp + secondaryamp)))
print('    Sec Amp (rel): {:.3f}'.format(secondaryamp/ (primaryamp + secondaryamp)))


# Plot of the Uranus beam map.
# Create fits file
from starlink import convert
if os.path.isfile('uranus_deep.fits'):
    os.remove('uranus_deep.fits')
convert.ndf2fits(uranusdeep, 'uranus_deep.fits')
if os.path.isfile('uranus_beamfit_residuals.fits'):
    os.remove('uranus_beamfit_residuals.fits')
convert.ndf2fits('uranus_beamfit_residuals.sdf', 'uranus_beamfit_residuals.fits')

hdu_deep = fits.open('uranus_deep.fits')[0]
datauranus = hdu_deep.data[0,:,:]

fig = plt.figure(figsize=(6.5,2.5))
ax = fig.add_subplot(121)
smax = 2.32894e-05
smin=-4.86899e-07
norm = ImageNormalize(vmin=smin, vmax=smax, stretch=astropy.visualization.SqrtStretch()) 
ax.imshow(datauranus, cmap='test', origin='lower', interpolation='none', norm=norm)

# zoom in on 100x100 central region, and set up appropriate ticks
xmax = 248
ymax = 219
ax.set_xlim(xmax - 75, xmax + 75)
ax.set_ylim(ymax - 75, ymax + 75)

ticks = np.arange(-350,351, 50)
ax.set_xticks(xmax + (ticks/3.22))
ax.set_xticklabels(ticks)

ax.set_yticks(ymax + (ticks/3.22))
ax.set_yticklabels(ticks)

ax.set_xlabel(u'$\delta$ RA (")')
ax.set_ylabel(u'$\delta$ Dec (")')
ax.set_xlim(xmax - 75, xmax + 75)
ax.set_ylim(ymax - 75, ymax + 75)

ax2 = fig.add_subplot(122)

# Uranus cut through.
cutthrough1 = datauranus[ymax,:]
cutthrough2 = datauranus[:,xmax]

posvals1 = 3.22*(np.arange(len(cutthrough1)) - xmax)
posvals2 = 3.22*(np.arange(len(cutthrough2)) - ymax)
ax2.plot(posvals1, cutthrough1/np.nanmax(cutthrough1), color='black', drawstyle='steps', label='RA')
ax2.plot(posvals2, cutthrough2/np.nanmax(cutthrough2), drawstyle='steps', label='Dec')

def Gauss(x, amp, fwhm):
    return amp * np.exp( (-4 * np.log(2) * x**2)/(fwhm**2))

resid_hdu=fits.open('uranus_beamfit_residuals.fits')[0]
residuals = resid_hdu.data[0,:,:]
cutresiduals1 = residuals[ymax,:]
cutresiduals2 = residuals[:, xmax]

fit = Gauss(posvals1, primaryamp, (beamfit_uranus.majfwhm[0]*u.radian).to(u.arcsec).value) +\
            Gauss(posvals1, secondaryamp, (beamfit_uranus.majfwhm[2]*u.radian).to(u.arcsec).value)
fit_primarybeam = Gauss(posvals1, primaryamp, (beamfit_uranus.majfwhm[0]*u.radian).to(u.arcsec).value)

ax2.plot(posvals1, fit/np.nanmax(fit), label='Beam fit')

ax2.legend(frameon=False)
ax2.set_xlabel('$\delta$ (")')
ax2.set_ylim(-0.001, 0.03)
ax2.set_xlim(-150,150)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
fig.subplots_adjust(wspace=0.6)
fig.savefig('beamfit-uranus.pdf', bbox_inches='tight', pad_inches=0.01)


uranuspeak = 0.11
uranus_center = 50
crl618peak = 0.009
crl618_center = 150
x = np.arange(0, 200)
noise = np.random.normal(np.zeros(x.shape), scale=0.1 * crl618peak)
fwhm = 14
rms =np.sqrt(np.mean(noise**2))



def Gaussian(x, fwhm, peak, center):
    exponent = (-1.0 * 4 * np.log(2) * (x - center)**2)/ (fwhm**2)
    return peak * np.e **(exponent)


gauss_uranus = Gaussian(x, fwhm, uranuspeak, uranus_center)
gauss_crl618 = Gaussian(x, fwhm, crl618peak, crl618_center)

fig = plt.figure(figsize=(2.5,2))
ax = fig.add_subplot(111)

ax.plot(x, gauss_uranus + gauss_crl618, color='black')
ax.plot(x, noise + gauss_uranus + gauss_crl618, color='red')

colluranus = matplotlib.collections.BrokenBarHCollection.span_where(x,
                                                                    ymin=3*rms,
                                                                    ymax=max(gauss_uranus + noise),
                                                                    where=gauss_uranus  > 3*rms,
                                                                    label='AST masked',
                                                                    facecolor='black', edgecolor='None',
                                                                    alpha=0.2)

collcrl618 = matplotlib.collections.BrokenBarHCollection.span_where(x,
                                                                    ymin=3*rms,
                                                                    ymax=max(gauss_crl618 + noise),
                                                                    where=gauss_crl618  > 3*rms,
                                                                    label='AST masked',
                                                                    facecolor='black', edgecolor='None',
                                                                    alpha=0.2)
ax.add_collection(colluranus)
ax.add_collection(collcrl618)


crl618_percentinside = 100* (gauss_crl618[gauss_crl618 > 3*rms].sum())/(gauss_crl618.sum())
uranus_percentinside = 100* (gauss_uranus[gauss_uranus > 3*rms].sum())/(gauss_uranus.sum())

t1 = ax.text(uranus_center, 0.04, '{:.1F}%\n inside AST mask'.format(uranus_percentinside), va='bottom', ha='center')
t2 = ax.text(crl618_center, crl618peak+0.001, '{:.1F}%\n inside AST mask'.format(crl618_percentinside), va='bottom', ha='center')

ax.spines['bottom'].set_position(('data', 0.0))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('Offset (arcseconds)')
ax.set_ylabel('Intensity (pW)')
ax.set_ylim(top=0.055)
fig.set_tight_layout(True)
fig.savefig('toymodel-astmasking.pdf', bbox_inches='tight', pad_inches=0.01)
