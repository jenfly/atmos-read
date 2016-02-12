import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import numpy as np
import xray
import pandas as pd
import matplotlib.pyplot as plt

import atmos as atm
import merra

# ----------------------------------------------------------------------
# Read daily data and calculate d/dp

varnm = 'U'
lon1, lon2 = 40, 120
plev = 200
pmin, pmax = 100, 300
subset_dict = {'XDim' : (lon1, lon2), 'Height' : (pmin, pmax)}
subset = '_40E-120E_90S-90N'

datadir = atm.homedir() + 'datastore/merra/daily/'
years = np.arange(1979, 2015)
months = np.arange(1, 13)

def savefile(datadir, varnm, plev, subset, year, month):
    filestr = '%smerra_D%sDP%d%s_%d%02d.nc'
    return filestr % (datadir, varnm, plev, subset, year, month)

def daily_calc(varnm, year, month, day, urls, subset_dict, plev):
    """Read daily data from file, compute d/dp and extract pressure level.
    """
    url = urls['%d%02d%02d' % (year, month, day)]
    print('Reading ' + url)
    with xray.open_dataset(url) as ds:
        var = atm.subset(ds[varnm], subset_dict, copy=False)
        name, attrs, _, _  = atm.meta(var)
        var = var.mean(dim='TIME')

    print('Computing d/dp')
    pres = var['Height']
    pres = atm.pres_convert(pres, pres.attrs['units'], 'Pa')
    dvar_dp = atm.gradient(var, pres, axis=0)
    # dp = np.gradient(pres)
    # dims = var.shape
    # dvar_dp = np.nan * var
    # for i in range(dims[1]):
    #     for j in range(dims[2]):
    #         dvar_dp.values[:, i, j] = np.gradient(var[:, i, j], dp)

    dvar_dp = atm.subset(dvar_dp, {'Height' : (plev, plev)}, copy=False)
    attrs['long_name'] = 'd/dp of ' + attrs['long_name']
    attrs['standard_name'] = 'd/dp of ' + attrs['standard_name']
    attrs['units'] = ('(%s)/Pa' % attrs['units'])
    dvar_dp.name = 'D%sDP' % name
    dvar_dp.attrs = attrs
    return dvar_dp

for year in years:
    urls = merra.merra_urls([year])
    for month in months:
        days = np.arange(1, atm.days_this_month(year, month) + 1)
        for d, day in enumerate(days):
            var_in = daily_calc(varnm, year, month, day, urls, subset_dict,
                                plev)
            jday = atm.mmdd_to_jday(month, day,year)
            var_in = atm.expand_dims(var_in, 'day', jday, axis=0)
            if d == 0:
                var = var_in
            else:
                var = xray.concat([var, var_in], dim='day')
        var = atm.squeeze(var)
        filenm = savefile(datadir, varnm, plev, subset, year, month)
        print('Saving to ' + filenm)
        atm.save_nc(filenm, var)

# ----------------------------------------------------------------------
# Consolidate monthly files into yearly files

def yrlyfile(datadir, varnm, year, subset):
    return '%smerra_%s%s_%d.nc' % (datadir, varnm, subset, year)

varname = 'D%sDP%d' % (varnm, plev)
for year in years:
    files = [savefile(datadir, varnm, plev, subset, year, m) for m in months]
    var = atm.load_concat(files, 'D%sDP' % varnm, concat_dim='day')
    filenm = yrlyfile(datadir, varname, year, subset)
    print('Saving to ' + filenm)
    atm.save_nc(filenm, var)
