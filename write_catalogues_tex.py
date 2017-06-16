from astropy.table import Table



extentcat = Table.read('jcmts850um_extent-cat030318_pub_000.fits')

extentcat['RA'].format = '%.4F'
extentcat['DEC'].format = '%.4F'
extentcat['TOTAL_FLUX'].format = '%.3E'
extentcat['PEAK_FLUX'].format = '%.3E'
extentcat['AREA'].format = '%.3E'
# Fix up units so they aren't impossible to compile in latex.

units={'RA':'deg', 'DEC': 'deg',
       'PEAK_FLUX': r'mJy\,arcsec$^{-2}$',
       'AREA':r'arcsec$^{-2}$'}

caption = r'An excerpt from the catalog of extents for Tile 30318.\label{tab:extents}'

extentcat[0:5].write('extentcat-t30318.tex',
                     format='ascii.aastex',
                     caption=caption,
                     col_align='c c c c c c p{7cm}',
                     latexdict={'units':units})

peakcat = Table.read('jcmts850um_peak-cat030318_pub_000.fits')
# Fix up units so they aren't impossible to compile in latex.
peakcat['RA'].format = '%.4F'
peakcat['DEC'].format = '%.4F'
peakcat['PEAK_FLUX'].format = '%.3E'


caption = r'An excerpt from the catalog of peaks for Tile 30318.\label{tab:peaks}'

peakcat[0:40].write('peakcat-t30318.tex',
                     format='ascii.aastex',
                     caption=caption,
                     latexdict={'units':units})
