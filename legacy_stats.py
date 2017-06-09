import logging
logger = logging.getLogger(__name__)
import numpy as np
from jsa_proc.config import get_database
from jsa_proc.state import JSAProcState
from jsa_proc.qastate import JSAQAState
from astropy.table import Table

db = get_database()

# Find all completed jobs
allcompleted = db.find_jobs(task=['hpx-s2-850-r2', 'hpx-s2-850-r1'], state=JSAProcState.COMPLETE)

# Get observation info about them all
obsinfo = []
for i in allcompleted:
    obs = db.get_obs_info(i.id)
    if len(obs) > 1:
        logger.warning('job {}:Multiple observations found!'.format(i.id))
        obsinfo.append(obs)
    else:
        obsinfo.append(obs[0])

tiles = []
for i in allcompleted:
    tile = db.get_tilelist(i.id)
    tiles.append(list(tile))

tiles = np.array(tiles)

# ob

# What do we want to know?
# 1. Total time and number obs, all, pointing and science, then by coadded and not coadded.
# 2. Time and numberin JCMTCAL, Legacy projects, PI projects



# Masks
times = np.array([o.date_end - o.date_obs for o in obsinfo])

science = np.array([True if o.obstype=='science' else False for o in obsinfo])
pointing = np.array([True if o.obstype=='pointing' else False for o in obsinfo])


jcmtcal = np.array([True if ('CAL' in o.project or 'EC30' in o.project) else False for o in obsinfo])
jls = np.array([True if o.project.startswith('MJLS') else False for o in obsinfo])

#engineering = np.array([True if 'EC' in o.project else False for o in obsinfo])
other = np.invert((jcmtcal) | (jls) )

qastatus_good = np.array([True if j.qa_state in [JSAQAState.GOOD, JSAQAState.QUESTIONABLE] else False for j in allcompleted])

coaddmask = (qastatus_good) & (science)

# for each type of selection, show number obs, number coadd obs, total time, total coadd time, 


#selections are: all, all science, all pointing, all CAL, , all JLS, engineering, other

def get_stats(times, tiles, coaddmask, categorymask):

    alltiles = set([j for i in tiles[categorymask] for j in i])
    allcoaddtiles = set([j for i in tiles[categorymask & coaddmask] for j in i])
    # Get d
    numberobs = len(times[categorymask])
    totaltime = times[categorymask].sum()
    if totaltime:
        totaltime=totaltime.total_seconds()/(60.0*60.0)
    else:
        totaltime=0.0
    numbercoaddobs = (coaddmask & categorymask).sum()
    timecoaddobs = times[coaddmask & categorymask].sum()
    if timecoaddobs:
        timecoaddobs = timecoaddobs.total_seconds()/(60.0*60.0)
    else:
        timecoaddobs = 0.0
    return numberobs, totaltime, len(alltiles), numbercoaddobs, timecoaddobs, len(allcoaddtiles)


allmask = np.array([True] * len(allcompleted))


res = Table(names = ['Type', 'Num. obs.', 'Tot. obs. time (hrs)', 'Tiles', 'Num. coadded obs', 'Tot. coadd time (hrs)', 'Coadded tiles'],
            dtype = ['S10', int, float, int,  int, float, int],)
#            format = ['%s', '%i', '%.1f', '%i', '%i', '%.1f', '%i'])
res['Tot. obs. time (hrs)'].format = '{:.1f}'
res['Tot. coadd time (hrs)'].format = '{:.1f}'


res.add_row(('All',) +  get_stats(times, tiles, coaddmask, allmask))
res.add_row(('Science',) +  get_stats(times, tiles, coaddmask, science))
res.add_row(('Pointing',) +  get_stats(times, tiles, coaddmask, pointing))
res.add_row(('Calib.',) +  get_stats(times, tiles, coaddmask, (jcmtcal & science)))
res.add_row(('JLS',) +  get_stats(times, tiles, coaddmask, jls))
#res.add_row(('EC30',) +  get_stats(times, tiles, coaddmask, engineering))
res.add_row(('PI',) +  get_stats(times, tiles, coaddmask, other))

