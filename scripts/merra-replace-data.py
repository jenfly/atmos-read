"""
Replace corrupted data files with daily data re-downloaded with wget
"""

import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import shutil
import xarray as xray
import numpy as np
import collections
import time
import matplotlib.pyplot as plt
import pandas as pd
import atmos as atm
import precipdat
import merra

# ----------------------------------------------------------------------
datadir = '/net/eady/data1/jwalker/datastore/merra2/wget/'
savedir = '/net/eady/data1/jwalker/datastore/merra2/merged/'
probdata = pd.read_csv('scripts/merra_urls/merge_data.csv', index_col=0)

# For each corrupted data file:
# - load the corrupted data file
# - load the new downloaded file for the problem day
# - calculate d/dp and other stuff
# - merge the data for the affected day
# - save into data file for the year

def latlon_filestr(lat1, lat2, lon1, lon2):
    """Return nicely formatted string for lat-lon range."""
    latstr = atm.latlon_str(lat1, lat2, 'lat')
    lonstr = atm.latlon_str(lon1, lon2, 'lon')
    return lonstr + '_' + latstr

def latlon_data(var, lat1, lat2, lon1, lon2, plev=None):
    """Extract lat-lon subset of data."""
    name = var.name
    varnm = name
    subset_dict = {'lat' : (lat1, lat2), 'lon' : (lon1, lon2)}
    latlonstr = latlon_filestr(lat1, lat2, lon1, lon2)
    if plev is not None:
        name = name + '%d' % plev
        subset_dict['plev'] = (plev, plev)
    var = atm.subset(var, subset_dict, copy=False, squeeze=True)
    var.name = name
    var.attrs['filestr'] = '%s_%s' % (name, latlonstr)
    var.attrs['varnm'] = varnm
    return var

def pgradient(var, lat1, lat2, lon1, lon2, plev):
    """Return d/dp of a lat-lon variable."""
    pwidth = 100
    p1, p2 = plev - pwidth, plev + pwidth
    var = atm.subset(var, {'lat' : (lat1, lat2), 'lon' : (lon1, lon2),
                           'plev' : (p1, p2)}, copy=False, squeeze=True)
    latlonstr = latlon_filestr(lat1, lat2, lon1, lon2)
    attrs = var.attrs
    pname = atm.get_coord(var, 'plev', 'name')
    pdim = atm.get_coord(var, 'plev', 'dim')
    pres = var[pname]
    pres = atm.pres_convert(pres, pres.attrs['units'], 'Pa')
    dvar_dp = atm.gradient(var, pres, axis=pdim)
    dvar_dp = atm.subset(dvar_dp, {pname : (plev, plev)}, copy=False,
                         squeeze=True)
    varnm = 'D%sDP' % var.name
    name = '%s%d' % (varnm, plev)
    dvar_dp.name = name
    attrs['long_name'] = 'd/dp of ' + var.attrs['long_name']
    attrs['standard_name'] = 'd/dp of ' + var.attrs['standard_name']
    attrs['units'] = ('(%s)/Pa' % attrs['units'])
    attrs[pname] = plev
    attrs['filestr'] = '%s_%s' % (name, latlonstr)
    attrs['varnm'] = varnm
    dvar_dp.attrs = attrs
    return dvar_dp


def var_calcs(filenm, varnm, plev, latlon=(-90, 90, 40, 120)):
    """Process a single variable from a single day."""
    lat1, lat2, lon1, lon2 = latlon
    if varnm == 'DUDP':
        nm, dp = 'U', True
    elif varnm == 'DOMEGADP':
        nm, dp = 'OMEGA', True
    else:
        nm, dp = varnm, False
    with xray.open_dataset(filenm) as ds:
        var = ds[nm].load()
    if dp:
        print('Computing d/dp')
        var = pgradient(var, lat1, lat2, lon1, lon2, plev)
    else:
        var = latlon_data(var, lat1, lat2, lon1, lon2, plev)

    return var

def process_row(row, datadir, savedir):
    filenm1 = row['filename']
    year = row['year']
    varnm = row['varnm']
    plev = row['plev']
    jday = row['jday']
    filenm2 = datadir + row['datfile']
    savefile1 = filenm1
    savefile2 = savedir + os.path.split(filenm1)[1]

    print('%d, %s, plev=%d' % (year, varnm, plev))
    print('Reading original data from ' + filenm1)
    with xray.open_dataset(filenm1) as ds:
        var1 = ds[varnm].load()
    print('Processing new data from ' + filenm2)
    var2 = var_calcs(filenm2, varnm, plev)

    print('Merging data for jday %d' % jday)
    var = var1.copy()
    ind = jday - 1
    days = atm.get_coord(var1, 'day')
    if not days[ind] == jday:
        raise ValueError('Days not indexed from 1, need to edit code to handle')
    var[ind] = var2
    print('Saving to ' + savefile1)
    var.to_netcdf(savefile1)
    print('Saving to ' + savefile2)
    var.to_netcdf(savefile2)

    data = {'orig' : var1, 'new' : var2, 'merged' : var}
    return data

# Make a copy of each of the original files -- only run this code once!
# for filenm in probdata['filename']:
#     shutil.copyfile(filenm, filenm.replace('.nc', '_orig.nc'))

for i, row in probdata.iterrows():
    data = process_row(row, datadir, savedir)


