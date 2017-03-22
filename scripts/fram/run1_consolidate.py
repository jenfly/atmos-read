import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import xray
import numpy as np
import atmos as atm


# ----------------------------------------------------------------------
# Consolidate daily data files into yearly files

version = 'merra2'
years = [2010, 2011]
months = np.arange(1, 13)
datadir = atm.homedir() + 'eady/datastore/' + version + '/daily/'
filestr = datadir + '%s_%s.nc'

varnms = ['U', 'V', 'OMEGA', 'T', 'QV', 'H', 'DUDTANA', 'PS',
          'UFLXCPT', 'VFLXCPT', 'UFLXPHI', 'VFLXPHI']

def get_filename(var, version, datadir, year, month=None, day=None):
    """Return a filename for a variable."""
    filenm = datadir + version + '_' + var.attrs['filestr'] + '_%d' % year
    if month is not None:
        filenm = filenm + '%02d' % month
    if day is not None:
        filenm = filenm + '%02d' % day
    filenm = filenm + '.nc'
    return filenm

# Read data and concatenate
for year in years:
    dates = []
    dailyfiles = {}
    for month in months:
        days = range(1, atm.days_this_month(year, month) + 1)
        dates = dates + ['%d%02d%02d' % (year, month, d) for d in days]
    for nm in varnms:
        dailyfiles[nm] = [filestr % (nm, date) for date in dates]

    # Consolidate daily files into yearly files and delete daily files
    for nm in dailyfiles:
        data = atm.load_concat(dailyfiles[nm], concat_dim='day')
        for varnm in data.data_vars:
            var = data[varnm]
            filenm = get_filename(var, version, datadir, year)
            var.name = var.attrs.get('varnm', varnm)
            print('Saving to ' + filenm)
            atm.save_nc(filenm, var)
        print('Deleting daily files')
        for filenm in dailyfiles[nm]:
            print(filenm)
            os.remove(filenm)
