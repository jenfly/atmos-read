"""
Surface fluxes and vertically integrated variables:
---------------------------------------------------
plev, sector, dp = None, False, False
varnms = ['PRECTOT', 'EVAP', 'EFLUX', 'HFLUX', 'QLML', 'TLML', 'TQV',
          'UFLXQV', 'VFLXQV', 'VFLXCPT', 'VFLXPHI', 'DQVDT_ANA', 'PS', 'SLP']

3-D variables - lat-lon subsets and sector means
------------------------------------------------
plev = 850 | 200
sector, dp = False, False
varnms = ['V', 'T', 'QV', 'H', 'DUDTANA']

plev = 850 | 200
sector, dp = False, True
varnms = ['U', 'OMEGA']

plev = None
sector, dp = True, False
varnms = ['U', 'V', 'OMEGA', 'T', 'QV', 'H', 'DUDTANA']
"""

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

# version = 'merra'
# years = np.arange(1979, 2016)
version = 'merra2'
# years = np.arange(1980, 2016)
years = [1980, 1981]

basedir = atm.homedir() + 'datastore/' + version + '/daily/'
months = np.arange(1, 13)

varnms = ['U', 'V', 'PRECTOT']

latlon=(-90, 90, 40, 120)
plevs=(850, 200)
sector_lons=(60, 100)

nc_fmt = {'merra' : None, 'merra2' : 'NETCDF4_classic'}[version]
nc_eng = {'merra' : None, 'merra2' : 'netcdf4'}[version]

def get_filename(var, version, datadir, year, month=None):
    filenm = datadir + version + '_' + var.attrs['filestr']
    if month is not None:
        filenm = filenm + '_%d%02d.nc' % (year, month)
    else:
        filenm = filenm + '_%d.nc' % year
    return filenm

def latlon_filestr(lat1, lat2, lon1, lon2):
    latstr = atm.latlon_str(lat1, lat2, 'lat')
    lonstr = atm.latlon_str(lon1, lon2, 'lon')
    return lonstr + '_' + latstr

def latlon_data(var, lat1, lat2, lon1, lon2, plev=None):
    name = var.name
    subset_dict = {'lat' : (lat1, lat2), 'lon' : (lon1, lon2)}
    latlonstr = latlon_filestr(lat1, lat2, lon1, lon2)
    if plev is not None:
        name = name + '%d' % plev
        subset_dict['plev'] = (plev, plev)
    var = atm.subset(var, subset_dict, copy=False, squeeze=True)
    var.name = name
    var.attrs['filestr'] = '%s_%s' % (name, latlonstr)
    return var

def pgradient(var, lat1, lat2, lon1, lon2, plev):
    pwidth = 100
    p1, p2 = plev - pwidth, plev + pwidth
    var = atm.subset(var, {'lat' : (lat1, lat2), 'lon' : (lon1, lon2),
                           'plev' : (p1, p2)}, copy=False)
    latlonstr = latlon_filestr(lat1, lat2, lon1, lon2)
    attrs = var.attrs
    pname = atm.get_coord(var, 'plev', 'name')
    pdim = atm.get_coord(var, 'plev', 'dim')
    pres = var[pname]
    pres = atm.pres_convert(pres, pres.attrs['units'], 'Pa')
    dvar_dp = atm.gradient(var, pres, axis=pdim)
    dvar_dp = atm.subset(dvar_dp, {pname : (plev, plev)}, copy=False,
                         squeeze=True)
    name = 'D%sDP%d' % (var.name, plev)
    dvar_dp.name = name
    attrs['long_name'] = 'd/dp of ' + var.attrs['long_name']
    attrs['standard_name'] = 'd/dp of ' + var.attrs['standard_name']
    attrs['units'] = ('(%s)/Pa' % attrs['units'])
    attrs[pname] = plev
    attrs['filestr'] = '%s_%s' % (name, latlonstr)
    dvar_dp.attrs = attrs
    return dvar_dp

def sector_mean(var, lon1, lon2):
    name = var.name
    lonstr = atm.latlon_str(lon1, lon2, 'lon')
    varbar = atm.dim_mean(var, 'lon', lon1, lon2)
    varbar.name = name + '_' + 'SEC'
    varbar.attrs['varnm'] = name
    varbar.attrs['lonstr'] = lonstr
    varbar.attrs['filestr'] = '%s_sector_%s' % (name, lonstr)
    return varbar

def calc_data(var, jday=0, latlon=(-90, 90, 40, 120), plevs=(850, 200),
             sector_lons=(60, 100)):

    lat1, lat2, lon1, lon2 = latlon
    opts = merra.url_opts(var.name)
    vertical = opts['vertical']
    if vertical == 'X':
        plevs = [None]
    if var.name in ['U', 'OMEGA']:
        dp = True
    else:
        dp = False
    data = xray.Dataset()

    # Lat-lon data
    print('Lat-lon data')
    for plev in plevs:
        print('plev', plev)
        var_out = latlon_data(var, lat1, lat2, lon1, lon2, plev)
        data[var_out.name] = var_out
        if dp:
            print('Computing d/dp')
            var_out = pgradient(var, lat1, lat2, lon1, lon2, plev)
        data[var_out.name] = var_out

    # Sector and zonal mean data
    print('Computing zonal mean')
    var_out = sector_mean(var, 0, 360)
    data[var_out.name] = var_out
    if vertical == 'P':
        print('Computing sector mean')
        var_out = sector_mean(var, sector_lons[0], sector_lons[1])
        data[var_out.name] = var_out

    # Compute daily data from subdaily data
    nperday = len(atm.get_coord(data, 'time'))
    data = atm.daily_from_subdaily(data, nperday, dayname='day',
                                   dayvals=[jday])

    return data


def get_kw(jdays, latlon, plevs, sector_lons):
    kw = [{'jday' : jday, 'latlon' : latlon, 'plevs' : plevs,
           'sector_lons' : sector_lons} for jday in jdays]
    return kw


# Read data and concatenate
for year in years:
    for varnm in varnms:
        for month in months:
            url_dict = merra.get_urls(year, month, version, varnm)
            days = range(1, atm.days_this_month(year, month) + 1)
            jdays = atm.season_days(atm.month_str(month), atm.isleap(year))
            urls = [url_dict['%d%02d%02d' % (year, month, day)] for day in days]
            func_kw = get_kw(jdays, latlon, plevs, sector_lons)
            data = atm.load_concat(urls, varnm, concat_dim='day', func=calc_data,
                                   func_kw=func_kw)
            if isinstance(data, xray.DataArray):
                data = data.to_dataset()

            # Save monthly files

            for nm in data.data_vars:
                var = data[nm]
                filenm = get_filename(var, version, datadir, year, month)

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
        for nm in data.data_vars:
            files = [get_filename(data[nm], version, datadir, year, m) for m in months]
            var = atm.load_concat(files, varnm, concat_dim='day')
            filenm = get_filename(data[nm], version, datadir, year)
            print('Saving to ' + filenm)
            var.to_dataset().to_netcdf(filenm)
            print('Deleting monthly files:')
            for filenm in files:
                print(filenm)
                os.remove(filenm)



#
#
# def sector_and_zonal_mean(var, lon1=None, lon2=None, incl_global=True):
#     """Return variable mean over sector and (optional) global zonal mean."""
#     name = var.name
#     # -- Sector mean
#     varbar = atm.dim_mean(var, 'lon', lon1, lon2)
#     varbar.name = name + '_' + 'SEC'
#     varbar.attrs['varnm'] = name
#     varbar.attrs['lonstr'] = atm.latlon_str(lon1, lon2, 'lon')
#     # -- Global zonal mean
#     if incl_global:
#         varzon = atm.dim_mean(var, 'lon')
#         varzon.name = name + '_' + 'ZON'
#         varzon.attrs['varnm'] = name
#         varzon.attrs['lonstr'] = atm.latlon_str(0, 360, 'lon')
#         data_out = xray.Dataset({varzon.name : varzon, varbar.name: varbar})
#     else:
#         data_out = varbar
#     return data_out

#
#
#
# def get_opts(version, plev, lon1, lon2, lat1, lat2, sector, sector_lon1,
#              sector_lon2, dp):
#     subset_dict, subset = {}, ''
#     if plev is not None:
#         subset = '%d' % plev + subset
#         if dp:
#             subset_dict['plev'] = (plev - 100, plev + 100)
#         else:
#             subset_dict['plev'] = (plev, plev)
#     if sector:
#         subset_dict = None
#         func = sector_and_zonal_mean
#         func_kw = {'lon1' : sector_lon1, 'lon2' : sector_lon2}
#     else:
#         subset_dict['lon'] = (lon1, lon2)
#         subset_dict['lat'] = (lat1, lat2)
#         for d1, d2, nm in zip([lon1, lat1], [lon2, lat2], ['lon', 'lat']):
#             if d1 is not None:
#                 subset = subset + '_' + atm.latlon_str(d1, d2, nm)
#         if dp:
#             func = var_and_dp
#             func_kw = {'plev' : plev}
#         else:
#             func, func_kw = None, None
#
#     time_dim = {'merra' : 'TIME', 'merra2' : 'time'}[version]
#     nc_fmt = {'merra' : None, 'merra2' : 'NETCDF4_classic'}[version]
#     nc_eng = {'merra' : None, 'merra2' : 'netcdf4'}[version]
#
#     return subset_dict, subset, time_dim, nc_fmt, nc_eng, func, func_kw
#
#
# opts = get_opts(version, plev, lon1, lon2, lat1, lat2, sector, sector_lon1,
#                 sector_lon2, dp)
# subset_dict, subset, time_dim, nc_fmt, nc_eng, func, func_kw = opts
