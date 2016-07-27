import os
from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
import astropy.coordinates as coord 
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

g34plot = None
crl618plot = True
cosmosmap = None

s=matplotlib._cm.cubehelix(s=0, r=-0.5, gamma=0.75)
matplotlib.cm.register_cmap(name='test', data=s, lut=128)

cmap = matplotlib.cm.get_cmap('test')
cmap.set_bad('0.7')

cmap_variance = matplotlib.cm.get_cmap('gist_heat_r')
cmap_variance.set_bad('0.7')

def formataxes(ax):
    ax.grid(color='white', alpha=0.1, linestyle='solid', linewidth=1.0)
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

def add_colorbar(ax):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0, axes_class=matplotlib.axes.Axes)
    cbar = plt.colorbar(mappable=ax.images[0], cax=cax, orientation='horizontal')
    cbar.set_label(r'mJy/arcsec\textsuperscript{2}')
    cax.xaxis.tick_top()
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation=90)
    cbar.ax.xaxis.set_label_position('top')

if g34plot:
    # Plot tile 30318 (G34.3)
    hdu = fits.open('30318_2d.fits')[0]
    wcs = WCS(hdu.header)
    fig = plt.figure(figsize=(6.6,11))
    ax = WCSAxesSubplot(fig, 1,2,1, wcs=WCS(hdu.header))
    fig.add_axes(ax)

    # color image.
    norm = ImageNormalize(vmin=-0.1, vmax=20, stretch=astropy.visualization.SqrtStretch())
 
    ax.imshow(hdu.data, cmap='test', origin='lower', norm=norm, interpolation='none')

    formataxes(ax)

    # Label
    ax.text(0.95,0.9, r'Tile 30318 (G34.27+0.15)' + '\n' + r'850$\mu$m emission', transform=ax.transAxes, size='medium', ha='right',va='top')
    # Colorbar
    add_colorbar(ax)

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
    formataxes(ax2)
    add_colorbar(ax2)


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


# Plot tiles 1238 and 1244 (CRL618 -- tile boundary right through middle)
#cp /net/kamaka/export/data/jsa_proc/output/000/000194/000194622/jcmts850um_healpix001238_pub_000.fits
#cp /net/kamaka/export/data/jsa_proc/output/000/000194/000194624/jcmts850um_healpix001244_pub_000.fits .

if crl618plot:

    t1238 = 'jcmts850um_healpix001238_pub_000.fits'
    t1244 = 'jcmts850um_healpix001244_pub_000.fits'
    pasted = 'crl618_pasted.fits'
    pastedextents = 'crl618_extentmask_pasted.fits'




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
    #fig.subplots_adjust(wspace=0.09)
    fig.subplots_adjust(wspace=0.035)
    plt.draw()
    #fig.set_tight_layout(True)
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
    ax.coords[0].set_separator((r'\textsuperscript{h}',r'\textsuperscript{m}',r'\textsuperscript{s}'))
    ax.coords[1].set_axislabel('Dec -- HPX', minpad=-0.5)
    ax.coords[0].set_axislabel('RA -- HPX')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="180%", pad=0.1, axes_class=matplotlib.axes.Axes)
    xvals = np.arange(201) -100
    cax.plot(xvals, data[229-100:229+101, 286] - 0.02, color='black', drawstyle='default')
    cax.plot(xvals, data[229-100:229+101, 286] - 0.02, color='cyan', drawstyle='default', linestyle='dotted')
    cax.plot(xvals, data[229,286-100:286+101]+0.02, color='black', drawstyle='default', linestyle='solid')
    cax.plot(xvals, data[229,286-100:286+101]+0.02, color='red', drawstyle='default', linestyle='dotted')
    cax.set_ylim(-0.06, 0.06)
    ax.vlines(286, 229-100, 229+100, color='cyan', linestyle='dotted')
    ax.hlines(229, 286-100, 286+100, color='red',linestyle='dotted')
    cax.tick_params(labelright=True, labelleft=False)
    cax.hlines(-0.02, cax.get_xlim()[0], cax.get_xlim()[1], color='0.7')
    cax.hlines(+0.02, cax.get_xlim()[0], cax.get_xlim()[1], color='0.7')

    cax.set_ylabel(r'Jy\,arcsec\textsuperscript{-2}')
    cax.yaxis.set_label_position('right')
    cax.set_xlabel('Pixels (offset from peak)')
    fig.tight_layout()
    fig.savefig('crl618-sourceonly.pdf', bbox_inches='tight', pad_inches=0.05)


if cosmosmap:
    # Deep cosmos region.
    mapfile = 'jcmts850um_healpix027258_pub_000.fits'

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
    ax2.contour(var, levels=[1.75e-3, 2e-3, 2.5e-3, 5e-3,7.5e-3,10e-3], colors='0.7')

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
    conts=ax.contour(var, levels=[5e-3], colors='white')


    ax.set_xlim(117.38899509591401, 433.24228632982658)
    ax.set_ylim(70.910809319165651, 384.07598743619394)

    ax2 = WCSAxesSubplot(fig, 1,2,2, wcs=wcs)
    fig.add_axes(ax2)
    norm = ImageNormalize(vmin=-0.0097, vmax=0.029, stretch=astropy.visualization.LinearStretch())
    ax2.imshow(data, cmap='gist_gray_r', origin='lower', interpolation='none', norm=norm)
    conts=ax2.contour(var, levels=[5e-3], colors='black')
    ax2.set_xlim(117.38899509591401, 433.24228632982658)
    ax2.set_ylim(70.910809319165651, 384.07598743619394)



    # Catalogs
    # Casey et al paper: UH observations. J/MNRAS/436/1919/table2
    uhcat = 'vizier_votable.vot'

    # CLS DR1 catalogs:
    clscat = 'S2CLS_CATALOGUE_DR1.FITS'

    # Our peak catalog
    lrcat = 'jcmts850um_peak-cat027258_pub_000.fits'

    uhtab = Table.read(uhcat)

    # Get decimal degree ra and dec
    dec=uhtab['DEJ2000']
    dec=[i.decode() for i in dec]
    dec=coord.Angle(dec, unit=u.degree)

    ra=uhtab['RAJ2000']
    ra=[i.decode() for i in ra]
    ra=coord.Angle(ra, unit=u.hour)
    ra=ra.to(unit=u.degree)
    ax2.scatter(ra, dec, s=45, transform=ax2.get_transform('fk5'), label='UH', color='red', marker='o', facecolor='none',)


    clstab=Table.read(clscat)
    ra=clstab['RA_DEG']
    dec=clstab['Dec_DEG']
    ax2.scatter(ra, dec, s=45, transform=ax2.get_transform('fk5'), label='LCS', color='blue', marker='^', facecolor='none', zorder=6)


    lrtab = Table.read(lrcat)
    ax2.scatter(lrtab['RA'], lrtab['DEC'], s=45, transform=ax2.get_transform('fk5'), color='cyan', marker='x', zorder=7, label='LR')
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

    fig.savefig('27258-zoomin.pdf', bbox_inches='tight', pad_inches=0.05)
