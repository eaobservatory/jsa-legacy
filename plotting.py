import os

from astropy.table import Table

import matplotlib
#matplotlib.use('svg')

import matplotlib.pyplot as plt

# from starlink import wrapper, kappa, convert

# if not wrapper.starpath:
#     wrapper.change_starpath('/stardev')

# matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{apjfonts},\usepackage{sfmath}'


# g34map = 'jcmts850um_healpix030318_pub_000.fits'
# s=matplotlib._cm.cubehelix(s=0, r=-0.5, gamma=0.75)
# matplotlib.cm.register_cmap(name='test', data=s, lut=128)

# im = aplpy.FITSFigure(g34map)
# im.tick_labels.set_xformat('hh:mm:ss')
# im.tick_labels.set_yformat('dd:mm')

# im.show_colorscale(cmap='test', vmin=-0.1, vmax=20, stretch='sqrt')
# im.set_nan_color('0.7')
# im.add_colorbar(pad=0)
# im.colorbar.set_axis_label_text(r'mJy/arcsec\textsuperscript{2}')


# im.add_grid()
# im.grid.set_color('white')
# im.grid.set_alpha(0.2)
# im.ticks.set_color('black')

# im.add_label(0.95,0.9, r'{Tile 30318\\G34.27+0.15}', relative=True, layer='title', size='medium', ha='right')

# im.frame.set_linewidth(0.5)
# im.colorbar.set_frame_linewidth(0.5)
# im.axis_labels.set_ypad(-10)
# im._figure.set_tight_layout(tight=True)

# im2 = aplpy.FITSFigure(g34map)
# im2.tick_labels.set_xformat('hh:mm:ss')
# im2.tick_labels.set_yformat('dd:mm')
# im2.ticks.set_color('black')
# im2.show_colorscale(cmap='gist_yarg', vmin=-0.1, vmax=20, stretch='sqrt')


# # Turn mask into correct file.
# extentmask = 'jcmts850um_extent-mask030318_pub_000.fits'
# peakcat = 'jcmts850um_peak-cat030318_pub_000.fits'
# contourmask = '30318_extent_contourmask.fits'

# # Create file for contouring (change bad values to 0, set all other values to 1)
# res = kappa.nomagic(extentmask, 'tempmask', 0)
# res = kappa.thresh('tempmask', contourmask, thrlo=0, thrhi=0.5, newlo=0, newhi=1)
# os.remove('tempmask.sdf')

# # Find positions of peaks.
# peakcat = Table.read(peakcat)
# ra = peakcat['RA']
# dec = peakcat['DEC']
# im2.add_colorbar(pad=0)
# im2.colorbar.set_frame_linewidth(0.5)
# im2.colorbar.set_axis_label_text(r'mJy/arcsec\textsuperscript{2}')


# im2.show_contour(contourmask, levels=[0.5], colors='red', linewidths=1.0, layer='contours')
# contours = im2.get_layer('contours')
# contours.collections[0].set_label('Extents')

# im2.show_markers(ra, dec, layer='peaks', marker='x', zorder=100, c='red', label='Peaks')

# im2._ax1.legend(frameon=False, fontsize='medium')
# im2._figure.set_tight_layout(tight=True)

import astropy
import astropy.units as u
from astropy.table import Table
import wcsaxes
from astropy.io import fits
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from astropy.wcs import WCS
from wcsaxes import WCSAxes, WCSAxesSubplot

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{apjfonts},\usepackage{sfmath}'
matplotlib.rcParams['figure.autolayout'] = True
matplotlib.rcParams['font.size'] = 7

#plt.ion()
hdu = fits.open('30318_2d.fits')[0]
wcs = WCS(hdu.header)

fig = plt.figure(figsize=(6.6,11))

ax = WCSAxesSubplot(fig, 1,2,1, wcs=WCS(hdu.header))

fig.add_axes(ax)

# color image.
norm = ImageNormalize(vmin=-0.1, vmax=20, stretch=astropy.visualization.SqrtStretch())
s=matplotlib._cm.cubehelix(s=0, r=-0.5, gamma=0.75)
matplotlib.cm.register_cmap(name='test', data=s, lut=128)
cmap = matplotlib.cm.get_cmap('test')
cmap.set_bad('0.7')

ax.imshow(hdu.data, cmap='test', origin='lower', norm=norm, interpolation='none')

# Grid
ax.grid(color='white', alpha=0.2, linestyle='solid', linewidth=1.0)

# Ticks.

racoord=ax.coords[0]
racoord.set_axislabel('RA -- HPX')
racoord.set_major_formatter('hh:mm:ss')
racoord.set_ticks_position('bt')
racoord.set_ticks(size=8)
racoord.set_ticks(spacing=5*(1/60.)*u.deg)
racoord.set_minor_frequency(6)
racoord.set_separator((r'\textsuperscript{h}',r'\textsuperscript{m}',r'\textsuperscript{s}'))
racoord.display_minor_ticks(True)

deccoord=ax.coords[1]
deccoord.set_axislabel('Dec -- HPX', minpad=0)
deccoord.set_major_formatter('dd:mm')
deccoord.set_ticks_position('lr')
deccoord.set_ticks(size=8)
deccoord.set_ticks(spacing=5*(1/60.)*u.deg)
deccoord.set_minor_frequency(6)
deccoord.display_minor_ticks(True)

# Label
ax.text(0.95,0.9, r'Tile 30318 (G34.27+0.15)' + '\n' + r'850$\mu$m emission', transform=ax.transAxes, size='medium', ha='right',va='top')
# Colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="5%", pad=0, axes_class=matplotlib.axes.Axes)
cbar = plt.colorbar(mappable=ax.images[0], cax=cax, orientation='horizontal')
cbar.set_label(r'mJy/arcsec\textsuperscript{2}')
cax.xaxis.tick_top()
cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation=90)
cbar.ax.xaxis.set_label_position('top')

# Second image: noise array.
ax2 = WCSAxesSubplot(fig, 1,2,2, wcs=WCS(hdu.header))
fig.add_axes(ax2)


# Get variance data.
hdu_var = fits.open('30318_2d.fits')[1]
wcs = WCS(hdu_var.header)


norm = ImageNormalize(vmin=np.sqrt(0.0004), vmax=0.2, stretch=astropy.visualization.PowerStretch(a=0.46*2))
cmap_variance = matplotlib.cm.get_cmap('gist_heat_r')
cmap_variance.set_bad('0.7')
ax2.imshow(np.sqrt(hdu_var.data), cmap=cmap_variance, origin='lower', norm=norm, interpolation='none')

racoord=ax2.coords[0]
racoord.set_axislabel('RA -- HPX')
racoord.set_major_formatter('hh:mm:ss')
racoord.set_ticks_position('bt')
racoord.set_ticks(size=8)
racoord.set_ticks(spacing=5*(1/60.)*u.deg)
racoord.set_minor_frequency(6)
racoord.set_separator((r'\textsuperscript{h}',r'\textsuperscript{m}',r'\textsuperscript{s}'))
racoord.display_minor_ticks(True)

deccoord=ax2.coords[1]
deccoord.set_axislabel('Dec -- HPX', minpad=0)
deccoord.set_major_formatter('dd:mm')
deccoord.set_ticks_position('lr')
deccoord.set_axislabel_position('r')
deccoord.set_ticklabel_position('r')
deccoord.set_ticks(size=8)
deccoord.set_ticks(spacing=5*(1/60.)*u.deg)
deccoord.set_minor_frequency(6)
deccoord.display_minor_ticks(True)


divider2 = make_axes_locatable(ax2)
cax = divider2.append_axes("top", size="5%", pad=0, axes_class=matplotlib.axes.Axes)
cbar = plt.colorbar(mappable=ax2.images[0], cax=cax, orientation='horizontal')
cbar.set_label(r'mJy/arcsec\textsuperscript{2}')
cax.xaxis.tick_top()
cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation=90)
cbar.ax.xaxis.set_label_position('top')
ax2.grid(color='white', alpha=0.2, linestyle='solid', linewidth=1.0)
ax2.text(0.95,0.9, r'Tile 30318 (G34.27+0.15)' + '\n' + r'850$\mu$m RMS noise',
         transform=ax2.transAxes, size='medium', ha='right',va='top')




fig.subplots_adjust(wspace=0.035)
fig.set_tight_layout(True)
fig.tight_layout()
plt.draw()

fig.savefig('tile30318-g34-coadd-noise.pdf', bbox_inches='tight', pad_inches=0.01)


# Second figure: extents and peaks.
#from starlink import wrapper, kappa, convert
#maskedfile = '30318_extent_contourmask_wlevels.fits'
#if os.path.isfile:
#    os.remove(maskedfile)
#kappa.ndfcopy('jcmts850um_extent-mask030318_pub_000.fits(,)', out='temp.sdf')
#kappa.nomagic('temp.sdf', repval=0, out=maskedfile)


matplotlib.rcParams['figure.autolayout'] = False
fig = plt.figure(figsize=(6.6,11))
wcs = WCS(hdu.header)
ax =  WCSAxesSubplot(fig, 1,2,1, wcs=WCS(hdu.header))
fig.add_axes(ax)
norm = ImageNormalize(vmin=-0.1, vmax=20, stretch=astropy.visualization.SqrtStretch())

maskedfile = '30318_extent_contourmask_2d.fits'
hdu_extents = fits.open(maskedfile)[0]
wcs_contour = WCS(hdu_extents.header)
ax.imshow(hdu.data, cmap='gist_yarg', origin='lower', norm=norm)
racoord=ax.coords[0]
racoord.set_axislabel('RA -- HPX')
racoord.set_major_formatter('hh:mm:ss')
racoord.set_ticks_position('bt')
racoord.set_ticklabel_position('bt')
racoord.set_ticks(size=8)
racoord.set_ticks(spacing=5*(1/60.)*u.deg)
racoord.set_minor_frequency(6)
racoord.set_separator((r'\textsuperscript{h}',r'\textsuperscript{m}',r'\textsuperscript{s}'))
racoord.display_minor_ticks(True)

deccoord=ax.coords[1]
deccoord.set_axislabel('Dec -- HPX', minpad=0)
deccoord.set_major_formatter('dd:mm')
deccoord.set_ticks_position('lr')
deccoord.set_axislabel_position('l')
deccoord.set_ticklabel_position('l')
deccoord.set_ticks(size=8)
deccoord.set_ticks(spacing=5*(1/60.)*u.deg)
deccoord.set_minor_frequency(6)
deccoord.display_minor_ticks(True)
ax.grid(color='black', alpha=0.1, linestyle='solid', linewidth=1.0)

ax.contour(hdu_extents.data, transform=ax.get_transform(wcs_contour), colors='red', levels=[0.5], linewidths=1.0, label='Extents')
ax.set_xlim(37.0, 618.468)
ax.set_ylim(106.41, 634.41)


# Peaks.
ax2 =  WCSAxesSubplot(fig, 1,2,2, wcs=WCS(hdu.header))
fig.add_axes(ax2)
norm = ImageNormalize(vmin=-0.1, vmax=20, stretch=astropy.visualization.SqrtStretch())

maskedfile = '30318_extent_contourmask_2d.fits'
hdu_extents = fits.open(maskedfile)[0]
wcs_contour = WCS(hdu_extents.header)
ax2.imshow(hdu.data, cmap='gist_yarg', origin='lower', norm=norm)
racoord=ax2.coords[0]
racoord.set_axislabel('RA -- HPX')
racoord.set_major_formatter('hh:mm:ss')
racoord.set_ticks_position('bt')
racoord.set_ticklabel_position('bt')
racoord.set_ticks(size=8)
racoord.set_ticks(spacing=5*(1/60.)*u.deg)
racoord.set_minor_frequency(6)
racoord.set_separator((r'\textsuperscript{h}',r'\textsuperscript{m}',r'\textsuperscript{s}'))
racoord.display_minor_ticks(True)

deccoord=ax2.coords[1]
deccoord.set_axislabel('Dec -- HPX', minpad=0)
deccoord.set_major_formatter('dd:mm')
deccoord.set_ticks_position('lr')
deccoord.set_axislabel_position('r')
deccoord.set_ticklabel_position('r')
deccoord.set_ticks(size=8)
deccoord.set_ticks(spacing=5*(1/60.)*u.deg)
deccoord.set_minor_frequency(6)
deccoord.display_minor_ticks(True)
ax2.grid(color='black', alpha=0.1, linestyle='solid', linewidth=1.0)

ax2.contour(hdu_extents.data, transform=ax2.get_transform(wcs_contour),alpha=0.5, colors='black', linstyle='dashed',levels=[0.5], linewidths=0.5, label='Extents')
ax2.set_xlim(37.0, 618.468)
ax2.set_ylim(106.41, 634.41)
peakcat = 'jcmts850um_peak-cat030318_pub_000.fits'

# # New figure, showing extents and peaks.
# matplotlib.rcParams['figure.autolayout'] = False
# # contours
# contourfile = '30318_extent_contourmask_2d.fits'
# hdu = fits.open(contourfile)[0]
# wcs_contour = WCS(hdu.header)
# ax2.set_autoscale_on(False)
# ax2.contour(hdu.data, transform=ax2.get_transform(wcs_contour), colors='red', levels=[0.5], linewidths=1.0, label='Extents')
# fig.set_tight_layout(tight=False)
# plt.draw()
fig.subplots_adjust(wspace=0.025)

# #peaks
# peakcat = 'jcmts850um_peak-cat030318_pub_000.fits'


peakcat = Table.read(peakcat)
ra = peakcat['RA']
dec = peakcat['DEC']

ax2.scatter(ra, dec, s=45, transform=ax2.get_transform('fk5'), color='red', marker='x', label='Peaks')
ax2.collections[0].set_label('Extents')
#ax2.legend(frameon=False, fontsize='medium')
ax2.grid(color='black', alpha=0.2, linestyle='solid', linewidth=1.0)
fig.tight_layout()
fig.savefig('tile30318-g34-extent-peak.pdf', bbox_inches='tight', pad_inches=0.05)
