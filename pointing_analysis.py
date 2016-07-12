from jsa_proc.config import get_database
from jsa_proc.admin.directories import get_log_dir

import os
import glob
import matplotlib
import matplotlib.pyplot as plt
plt.ion()
import numpy as np

def get_pointing_offsets(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    matchline = "Determined&nbsp;pointing&nbsp;offsets&nbsp;of&nbsp;"
    matchline = '<span CLASS="green">&nbsp;Shifting&nbsp;'
    offset_outputs = []
    matches = [i for i in lines if matchline in i]
    for l in matches:
        offsets = [float(i.strip(',')) for i in l.split('&nbsp;')[-4:-2]]

        if len(offsets) != 2:
            raise StandardError("Couldn't find the offsets!")
        offset_outputs.append(offsets)
        #print('{} {}'.format(offsets[0], offsets[1]))
    return offset_outputs

def plot_offset(offset_dict, ax, label, color, marker):
    offsets = np.array(offset_dict.values())
    ax.scatter(np.pi + np.arctan2(offsets[:,0], offsets[:,1]), np.sqrt(offsets[:,0]**2 + offsets[:,1]**2), marker=marker, label=label, color=color)




def find_offsets(jobs):
    offset_dict = {}
    for j in jobs:
        logdir = get_log_dir(j.id)
        oraclogs = glob.glob(os.path.join(logdir,'oracdr_*.html'))
        oraclogs.sort(key=os.path.getmtime)
        oraclog = oraclogs[-1]
        try:
            pointingoffsets =  get_pointing_offsets(oraclog)
            for p in pointingoffsets:
                offset_dict[j.id] = p
        except:
            sourcename = db.get_obs_info(j.id)[0].sourcename
            print('No pointing for %s, job=%i' % (sourcename,j.id))
    return offset_dict

def plot_offset(offset_dict, ax, label, color, marker):
    offsets = np.array(offset_dict.values())
    ax.scatter(np.pi + np.arctan2(offsets[:,0], offsets[:,1]), np.sqrt(offsets[:,0]**2 + offsets[:,1]**2), marker=marker, label=label, color=color)




db = get_database()
crl618_jobs = db.find_jobs(task=['hpx-s2-850-r1', 'hpx-s2-850-r2'], state='Y', obsquery={'project':'JCMTCAL', 'obstype':'science', 'sourcename':'CRL618', 'instrument':'SCUBA-2','subsys':'850' })
crl2688_jobs = db.find_jobs(task=['hpx-s2-850-r1', 'hpx-s2-850-r2'], state='Y', obsquery={'project':'JCMTCAL', 'obstype':'science', 'sourcename':'CRL2688', 'instrument':'SCUBA-2','subsys':'850' })
arp220_jobs = db.find_jobs(task=['hpx-s2-850-r1', 'hpx-s2-850-r2'], state='Y', obsquery={'project':'JCMTCAL', 'obstype':'science', 'sourcename':'Arp220', 'instrument':'SCUBA-2','subsys':'850' })
allcal_jobs = db.find_jobs(task=['hpx-s2-850-r1', 'hpx-s2-850-r2'], state='Y', obsquery={'project':'JCMTCAL', 'obstype':'science', 'instrument':'SCUBA-2','subsys':'850' })

print('Finding CRL 618 offsets')
crl618_offsets = find_offsets(crl618_jobs)
print('Finding CRL 2688 offsets')
crl2688_offsets = find_offsets(crl2688_jobs)
print('Finding Arp 220 offsets')
arp220_offsets = find_offsets(arp220_jobs)
print('Finding All Calibration offsets')
allcal_offsets = find_offsets(allcal_jobs)



crl618_mag = np.array(crl618_offsets.values())
crl618_mag = np.sqrt(crl618_mag[:,0]**2 + crl618_mag[:,1]**2)
crl2688_mag = np.array(crl2688_offsets.values())
crl2688_mag = np.sqrt(crl2688_mag[:,0]**2 + crl2688_mag[:,1]**2)
arp220_mag = np.array(arp220_offsets.values())
arp220_mag = np.sqrt(arp220_mag[:,0]**2 + arp220_mag[:,1]**2)
allcal_mag = np.array(allcal_offsets.values())
allcal_mag = np.sqrt(allcal_mag[:,0]**2 + allcal_mag[:,1]**2)

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{apjfonts},\usepackage{sfmath}'
matplotlib.rcParams['figure.autolayout'] = True
matplotlib.rcParams['font.size'] = 7
fig = plt.figure(figsize=(3.394, 3.5))
ax = fig.add_subplot(111)
values, bins, patches = ax.hist(allcal_mag, bins=25, range=(0,10), histtype='step', color='black', label='All Calibrations')
ax.hist(crl618_mag, bins=bins, histtype='step', color='red', hatch='\\\\\\\\', label='CRL 618')
ax.hist(crl2688_mag, bins=bins, histtype='step', color='blue', hatch='////', label='CRL 2688')
ax.hist(arp220_mag, bins=bins, histtype='step', color='purple', hatch='++', label='Arp 220')

ax.set_xlabel('Radial offset from known position (arcseconds)')
ax.set_ylabel('Number observations')
ylim = ax.get_ylim()
ax.vlines(np.median(allcal_mag), ylim[0], ylim[1], linestyle='dashed', color='black', label='Median (%.2f)'%(np.median(allcal_mag)))
ax.set_ylim(ylim)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig('pointing-offsets-by-source.pdf', bbox_inches='tight', pad_inches=0.01)

print('Mean: %.2f' % np.mean(allcal_mag))
print('Std: %.2f' % np.std(allcal_mag))
print('Median: %.2f' % np.median(allcal_mag))
print('CRL618 Mean: %.2f' % np.mean(crl618_mag))
print('CRL618 Std: %.2f' % np.std(crl618_mag))
print('CRL618 Median: %.2f' % np.median(crl618_mag))
print('CRL2688 Mean: %.2f' % np.mean(crl2688_mag))
print('CRL2688 Std: %.2f' % np.std(crl2688_mag))
print('CRL2688 Median: %.2f' % np.median(crl2688_mag))
print('ARP220 Mean: %.2f' % np.mean(arp220_mag))
print('ARP220 Std: %.2f' % np.std(arp220_mag))
print('ARP220 Median: %.2f' % np.median(arp220_mag))

#ax1 = fig.add_subplot(111, polar=True)
#ax2 = fig.add_subplot(222, polar=True)
#ax3 = fig.add_subplot(223, polar=True)
#ax4 = fig.add_subplot(224, polar=True)
#plot_offset(allcal_offsets, ax1, 'All Cal', 'black', 'x')
#plot_offset(crl618_offsets, ax1, 'CRL 618', 'orange', 'x')
#plot_offset(crl2688_offsets, ax1, 'CRL 2688', 'red', 'x')
#plot_offset(arp220_offsets, ax1, 'Arp 220', 'green', 'x')
# ax1.legend()
# ax1.set_title('CRL618 offsets (")')

# plot_offset(crl2688_point_offsets, ax2, 'CRL 2688 (pointings)', 'black', '+')
# plot_offset(crl2688_offsets, ax2, 'CRL 2688', 'orange', 'x')
# ax2.legend()
# ax2.set_title('CRL2688 offsets (")')
# fig.suptitle('Offset in measured source position from expected position, in reductions of pointings and JCMTCAL scans, as reported in the DR logs for jcmt-nightly and jcmt-reproc')

#fig.savefig('pointing_offsets_polarplot.pdf', bbox_inches='tight', pad_inches=0.1

# from astropy.table import Table

# vals = np.asarray(pointing_dict.values())

# total_offset = np.sqrt((vals[:,1]**2 + vals[:,0]**2))

# fields = ('id',
#  'job_id',
#  'obsid',
#  'obsidss',
#  'utdate',
#  'obsnum',
#  'instrument',
#  'backend',
#  'subsys',
#  'project',
#  'survey',
#  'scanmode',
#  'sourcename',
#  'obstype',
#  'association',
#  'date_obs',
#  'omp_status',
#  'tau',
#  'seeing',
#  'date_end')

# import datetime
# tablecols = ['job_id', 'obsid', 'utdate', 'obsnum', 'project', 'scanmode',
#              'sourcename', 'obstype', 'tau', 'date_obs', 'seeing', 'date_end',
#              'pointoffx', 'pointoffy', 'pointoff', 'omp_status']
# dtypes = ['i4', 'S16', datetime.date, 'i4', 'S20', 'S20', 'S20', 'S20', float, datetime.datetime, float, datetime.datetime, float, float, float, int]
# obstable = Table(names=tablecols, dtype=dtypes)
# for i in pointing_dict:
#     obsinfo = db.get_obs_info(i)[0]
#     obstable.add_row((i, obsinfo.obsid, obsinfo.utdate, obsinfo.obsnum,
#                      obsinfo.project, obsinfo.scanmode,
#                      obsinfo.sourcename, obsinfo.obstype,
#                      obsinfo.tau, obsinfo.date_obs,
#                      obsinfo.seeing,
#                      obsinfo.date_end,
#                      pointing_dict[i][0],
#                      pointing_dict[i][1],
#                      np.sqrt(pointing_dict[i][0]**2 + pointing_dict[i][1]**2),
#                      obsinfo.omp_status))
