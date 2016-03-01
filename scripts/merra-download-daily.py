import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import xray
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import atmos as atm
import precipdat
import merra

# ----------------------------------------------------------------------
# Download daily data

# version, years = 'merra2', np.arange(1980, 2016)
version, years = 'merra', np.arange(1979, 2016)
basedir = atm.homedir() + 'datastore/' + version + '/'
months = np.arange(1, 13)

sector = True
varnms = ['U']
plev = None

#sector = False
#varnms = ['DUDTANA']
#plev = 200

if sector:
    datadir =  basedir + 'daily_latpres/'
    lon1, lon2, lat1, lat2 = None, None, None, None
    sector_lon1, sector_lon2 = 60, 100
else:
    datadir = basedir + 'daily/'
    lon1, lon2, lat1, lat2 = 40, 120, -90, 90
    sector_lon1, sector_lon2 = None, None

def sector_and_zonal_mean(var, lon1=None, lon2=None, incl_global=True):
    name = var.name

    # Sector mean
    varbar = atm.dim_mean(var, 'lon', lon1, lon2)
    varbar.name = name + '_' + 'SEC'
    varbar.attrs['varnm'] = name
    varbar.attrs['lonstr'] = atm.latlon_str(lon1, lon2, 'lon')

    # Global zonal mean
    if incl_global:
        varzon = atm.dim_mean(var, 'lon')
        varzon.name = name + '_' + 'ZON'
        varzon.attrs['varnm'] = name
        varzon.attrs['lonstr'] = atm.latlon_str(0, 360, 'lon')
        data_out = xray.Dataset({varzon.name : varzon, varbar.name: varbar})
    else:
        data_out = varbar

    return data_out


def monthlyfile(datadir, varnm, year, month, subset):
    return '%smerra_%s%s_%d%02d.nc' % (datadir, varnm, subset, year, month)

def yrlyfile(datadir, varnm, year, subset):
    return '%smerra_%s%s_%d.nc' % (datadir, varnm, subset, year)

def get_opts(version, plev, lon1, lon2, lat1, lat2, sector, sector_lon1,
             sector_lon2):
    subset_dict, subset = {}, ''
    if plev is not None:
        subset_dict['plev'] = (plev, plev)
        subset = '%d' % plev + subset
    if not sector:
        func, func_kw = None, None
        subset_dict['lon'] = (lon1, lon2)
        subset_dict['lat'] = (lat1, lat2)
        for d1, d2, nm in zip([lon1, lat1], [lon2, lat2], ['lon', 'lat']):
            if d1 is not None:
                subset = subset + '_' + atm.latlon_str(d1, d2, nm)
    else:
        subset_dict = None
        func = sector_and_zonal_mean
        func_kw = {'lon1' : sector_lon1, 'lon2' : sector_lon2}

    time_dim = {'merra' : 'TIME', 'merra2' : 'time'}[version]
    nc_fmt = {'merra' : None, 'merra2' : 'NETCDF4_classic'}[version]
    nc_eng = {'merra' : None, 'merra2' : 'netcdf4'}[version]

    return subset_dict, subset, time_dim, nc_fmt, nc_eng, func, func_kw


opts = get_opts(version, plev, lon1, lon2, lat1, lat2, sector, sector_lon1,
                sector_lon2)
subset_dict, subset, time_dim, nc_fmt, nc_eng, func, func_kw = opts

# Read data and concatenate
for varnm in varnms:
    for year in years:
        for month in months:
            url_dict = merra.get_urls(year, month, version, varnm)
            days = range(1, atm.days_this_month(year, month) + 1)
            jdays = atm.season_days(atm.month_str(month), atm.isleap(year))
            urls = [url_dict['%d%02d%02d' % (year, month, day)] for day in days]
            data = atm.load_concat(urls, varnm, concat_dim=time_dim,
                                   subset_dict=subset_dict, func=func,
                                   func_kw=func_kw)
            nperday = len(data[time_dim]) / len(days)
            data = atm.daily_from_subdaily(data, nperday, dayname='day',
                                           dayvals=jdays)
            if isinstance(data, xray.DataArray):
                data = data.to_dataset()

            # Save monthly files
            for nm in data.data_vars:
                var = data[nm]
                if sector:
                    nm = nm.split('_')[0]
                    subset = '_' + var.attrs['lonstr']
                    var.name = nm
                filenm = monthlyfile(datadir, nm, year, month, subset)
                print('Saving to ' + filenm)
                ds = var.to_dataset()
                ds.to_netcdf(filenm, format=nc_fmt, engine=nc_eng)

                # Check if output is corrupted
                with xray.open_dataset(filenm) as ds_check:
                    print(ds.dims.keys())
                    print(ds.data_vars.keys())
                    if len(ds.data_vars.keys()) > 1:
                        raise ValueError('Corrupted monthly output file')

        # Consolidate monthly files into yearly files
        if sector:
            subset_list = [atm.latlon_str(sector_lon1, sector_lon2, 'lon'),
                           atm.latlon_str(0, 360, 'lon')]
        else:
            subset_list = [subset]
        for sub in subset_list:
            files = [monthlyfile(datadir, varnm, year, m, sub) for m in months]
            var = atm.load_concat(files, varnm, concat_dim='day')
            filenm = yrlyfile(datadir, varnm, year, sub)
            print('Saving to ' + filenm)
            var.to_dataset().to_netcdf(filenm)
            print('Deleting monthly files:')
            for filenm in files:
                print(filenm)
                os.remove(filenm)
