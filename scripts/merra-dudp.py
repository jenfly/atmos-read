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

#varnm = 'U'
varnm = 'OMEGA'
save_var = True     # True to save var and dvar_dp, False to save only dvar_dp
lon1, lon2 = 40, 120
plev = 200
pmin, pmax = 100, 300
subset_dict = {'XDim' : (lon1, lon2), 'Height' : (pmin, pmax)}
subset = '_40E-120E_90S-90N'

datadir = atm.homedir() + 'datastore/merra/daily/'
years = np.arange(1979, 2015)
months = np.arange(1, 13)

def monthfile(datadir, varnm, plev, subset, year, month):
    filestr = '%smerra_%s%d%s_%d%02d.nc'
    return filestr % (datadir, varnm, plev, subset, year, month)

def savemonth(var, datadir, varnm, plev, subset, year, month):
    filenm = monthfile(datadir, varnm, plev, subset, year, month)
    print('Saving to ' + filenm)
    atm.save_nc(filenm, var)

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
    dvar_dp = atm.subset(dvar_dp, {'Height' : (plev, plev)}, copy=False)
    dvar_dp = atm.squeeze(dvar_dp)
    attrs['long_name'] = 'd/dp of ' + attrs['long_name']
    attrs['standard_name'] = 'd/dp of ' + attrs['standard_name']
    attrs['units'] = ('(%s)/Pa' % attrs['units'])
    dvar_dp.name = 'D%sDP' % name
    dvar_dp.attrs = attrs
    var = atm.subset(var, {'Height' : (plev, plev)}, copy=False)
    var = atm.squeeze(var)

    # Add day dimension
    jday = atm.mmdd_to_jday(month, day, year)
    var = atm.expand_dims(var, 'day', jday, axis=0)
    dvar_dp = atm.expand_dims(dvar_dp, 'day', jday, axis=0)

    return var, dvar_dp

def consolidate_year(datadir, varnm, plev, subset, year, months):
    files = [monthfile(datadir, varnm, plev, subset, year, m) for m in months]
    var = atm.load_concat(files, varnm, concat_dim='day')
    filenm = '%smerra_%s%d%s_%d.nc' % (datadir, varnm, plev, subset, year)
    print('Saving to ' + filenm)
    atm.save_nc(filenm, var)

for year in years:
    urls = merra.merra_urls([year])
    for month in months:
        days = np.arange(1, atm.days_this_month(year, month) + 1)
        for d, day in enumerate(days):
            var_in, dvar_in = daily_calc(varnm, year, month, day, urls,
                                         subset_dict, plev)
            if d == 0:
                var, dvar = var_in, dvar_in
            else:
                var = xray.concat([var, var_in], dim='day')
                dvar = xray.concat([dvar, dvar_in], dim='day')
        if save_var:
            savemonth(var, datadir, varnm, plev, subset, year, month)
        savemonth(dvar, datadir, 'D%sDP' % varnm, plev, subset, year, month)

    # Consolidate monthly files into yearly files
    if save_var:
        consolidate_year(datadir, varnm, plev, subset, year, months)
    consolidate_year(datadir, 'D%sDP' % varnm, plev, subset, year, months)
