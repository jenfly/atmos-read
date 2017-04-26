import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import xarray as xray
import numpy as np
import atmos as atm

filenm1 = 'troubleshooting/data-check2.log'
filenm2 = 'troubleshooting/data-probs2.log'

version = 'merra2'
years = np.arange(1980, 2016)

datadir = atm.homedir() + 'datastore/' + version + '/daily/'
months = np.arange(1, 13)

#varnms = ['U', 'DUDP', 'OMEGA', 'DOMEGADP', 'V', 'H', 'DUDTANA']
#plevs = [1000,925,850,775,700,600,500,400,300,250,200,150,100,70,50,30,20]
varnms = ['TLML', 'QLML']
plevs = [None]

latlonstr = '40E-120E_90S-90N'

def get_filename(year, varnm, plev, datadir, version, latlonstr):
    if plev is None:
        nm = varnm
    else:
        nm = '%s%d' % (varnm, plev)
    filenm = ('%s%d/%s_%s_%s_%d.nc' %
              (datadir, year, version, nm, latlonstr, year))
    return filenm


def check_var(var, max_magnitude=1e10):
    problem_days = []
    days = atm.get_coord(var, 'day')
    for i, day in enumerate(days):
        #print(day)
        vals = abs(var[i].values)
        if np.nanmax(vals) > max_magnitude:
            problem_days.append(day)
    return problem_days


def initialize_logs(filenm1, filenm2):
    with open(filenm1, 'w') as f1:
        f1.write('filename, year, varnm, plev\n')

    with open(filenm2, 'w') as f2:
        f2.write('filename, year, varnm, plev, jday, yyyymmdd\n')


def write_logs(filenm1, filenm2, datafile, year, varnm, plev, problem_days):
    s = '%s, %d, %s' % (datafile, year, varnm)
    if plev is not None:
        s = s + ', %d' % plev
    if len(problem_days) == 0:
        with open(filenm1, 'a') as f1:
            f1.write(s + ' OK\n')
    else:
        with open(filenm1, 'a') as f1:
            f1.write(s + ' ***** PROBLEM ****\n')
        for jday in problem_days:
            month, day = atm.jday_to_mmdd(jday, year)
            datestr = '%d%02d%02d' % (year, month, day)
            with open(filenm2, 'a') as f2:
                f2.write(s + ', %d, %s\n' % (jday, datestr))


initialize_logs(filenm1, filenm2)
for year in years:
    for nm in varnms:
        for plev in plevs:
            datafile = get_filename(year, nm, plev, datadir, version, latlonstr)
            print('Loading ' + datafile)
            with xray.open_dataset(datafile) as ds:
                var = ds[nm]
                problem_days = check_var(var)
                write_logs(filenm1, filenm2, datafile, year, nm, plev, problem_days)