"""
3-D variables - lat-lon subsets and sector means
------------------------------------------------
varnms = ['U', 'V', 'OMEGA', 'T', 'QV', 'H', 'DUDTANA']

Surface fluxes and vertically integrated variables:
---------------------------------------------------
varnms = ['PRECTOT', 'EVAP', 'EFLUX', 'HFLUX', 'QLML', 'TLML', 'TQV',
          'UFLXQV', 'VFLXQV', 'VFLXCPT', 'VFLXPHI', 'DQVDT_ANA', 'PS', 'SLP']
"""

import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import xray
import numpy as np
import collections
import time
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

datadir = atm.homedir() + 'datastore/' + version + '/daily/'
months = np.arange(1, 13)

varnms = ['U', 'PRECTOT', 'V', 'TQV', 'EVAP', 'OMEGA']

dp_vars = ['U', 'OMEGA']

latlon=(-90, 90, 40, 120)
plevs=(850, 200)
sector_lons=(60, 100)

nc_fmt = {'merra' : None, 'merra2' : 'NETCDF4_classic'}[version]
nc_eng = {'merra' : None, 'merra2' : 'netcdf4'}[version]


def group_variables(varnms, version):
    """Group variables together according to URL."""
    def get_group(varnm, version):
        opts = merra.url_opts(varnm, version)
        group = '%s%s_%s_%s' % (opts['res'], opts['vertical'], opts['kind'],
                                opts['time_kind'])
        return group

    groups = {nm : get_group(nm, version) for nm in varnms}
    keys = set(groups.values())
    vargroups = collections.defaultdict(list)
    for nm, key in groups.iteritems():
        vargroups[key] += [nm]
    return vargroups

def get_filename(var, version, datadir, year, month=None, day=None):
    """Return a filename for a variable."""
    filenm = datadir + version + '_' + var.attrs['filestr'] + '_%d' % year
    if month is not None:
        filenm = filenm + '%02d' % month
    if day is not None:
        filenm = filenm + '%02d' % day
    filenm = filenm + '.nc'
    return filenm

def latlon_filestr(lat1, lat2, lon1, lon2):
    """Return nicely formatted string for lat-lon range."""
    latstr = atm.latlon_str(lat1, lat2, 'lat')
    lonstr = atm.latlon_str(lon1, lon2, 'lon')
    return lonstr + '_' + latstr

def latlon_data(var, lat1, lat2, lon1, lon2, plev=None):
    """Extract lat-lon subset of data."""
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
    """Return d/dp of a lat-lon variable."""
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
    """Return the sector mean of a variable."""
    name = var.name
    lonstr = atm.latlon_str(lon1, lon2, 'lon')
    varbar = atm.dim_mean(var, 'lon', lon1, lon2)
    varbar.name = name + '_' + 'SEC'
    varbar.attrs['varnm'] = name
    varbar.attrs['lonstr'] = lonstr
    varbar.attrs['filestr'] = '%s_sector_%s' % (name, lonstr)
    return varbar

def var_calcs(var, jday=0, latlon=(-90, 90, 40, 120), plevs=(850, 200),
              dp_vars=['U', 'OMEGA'], sector_lons=(60, 100)):
    """Process a single variable from a single day."""
    lat1, lat2, lon1, lon2 = latlon
    opts = merra.url_opts(var.name)
    vertical = opts['vertical']
    if vertical == 'X':
        plevs = [None]
    if dp_vars is not None and var.name in dp_vars:
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

    # Make sure output is in a Dataset
    if isinstance(data, xray.DataArray):
        data = data.to_dataset()

    return data

def all_data(ds, varnms, datadir, year, month, day, jday, calc_kw, nc_kw):
    """Process selected variables in a dataset and save each to file."""
    files = {}
    for nm in varnms:
        print(nm)
        data = var_calcs(ds[nm], jday, **calc_kw)
        filenm = '%s%s_%d%02d%02d.nc' % (datadir, nm, year, month, day)
        print('Saving to ' + filenm)
        data.to_netcdf(filenm, **nc_kw)
        files[nm] = filenm
    return files

def read_url(url, varnms, datadir, year, month, day, jday, calc_kw, nc_kw):
    """Open url and process selected variables."""
    # Number of times to attempt opening url (in case of server problems)
    NMAX = 3
    # Wait time (seconds) between attempts
    WAIT = 5

    print('Loading ' + url)
    attempt = 0
    while attempt < NMAX:
        try:
            with xray.open_dataset(url) as ds:
                files = all_data(ds, varnms, datadir, year, month, day, jday,
                                 calc_kw, nc_kw)
            attempt = NMAX
        except RuntimeError as err:
            attempt += 1
            if attempt < NMAX:
                print('Error reading file.  Attempting again in %d s' % WAIT)
                time.sleep(WAIT)
            else:
                raise err
    return files

def read_groups(url_dict, vargroups, datadir, year, month, day, jday, calc_kw,
                nc_kw):
    """Process variables for a day, grouped by URL."""
    files = {}
    for key, varids in vargroups.iteritems():
        url = url_dict[key]['%d%02d%02d' % (year, month, day)]
        datafiles = read_url(url, varids, datadir, year, month, day, jday,
                             calc_kw, nc_kw)
        files.update(datafiles)
    return files

def get_url_dict(year, month, version, vargroups):
    """Return dict of urls for the variable groups."""
    url_dict = {}
    for key in vargroups:
        nm = vargroups[key][0]
        url_dict[key] = merra.get_urls(year, month, version, nm)
    return url_dict

# Initial setup
vargroups = group_variables(varnms, version)
calc_kw = {'latlon' : latlon, 'plevs' : plevs, 'dp_vars' : dp_vars,
           'sector_lons' : sector_lons}
nc_kw = {'format' : nc_fmt, 'engine' : nc_eng}

# Read data and concatenate
for year in years:
    dailyfiles = collections.defaultdict(list)
    for month in months:
        url_dict = get_url_dict(year, month, version, vargroups)
        days = range(1, atm.days_this_month(year, month) + 1)
        jdays = atm.season_days(atm.month_str(month), atm.isleap(year))
        for day, jday in zip(days, jdays):
            files = read_groups(url_dict, vargroups, datadir, year, month, day,
                                jday, calc_kw, nc_kw)
            for nm in files:
                dailyfiles[nm] += [files[nm]]
    # Consolidate daily files into yearly files and delete daily files
    for nm in dailyfiles:
        data = atm.load_concat(dailyfiles[nm], concat_dim='day')
        for varnm in data.data_vars:
            filenm = get_filename(data[varnm], version, datadir, year)
            print('Saving to ' + filenm)
            atm.save_nc(filenm, data[varnm])
        print('Deleting daily files'):
        for filenm in dailyfiles[nm]:
            print(filenm)
            os.remove(filenm)



# ----------------------------------------------------------------------
# Read data and concatenate
# for year in years:
#     for varnm in varnms:
#         for month in months:
#             url_dict = merra.get_urls(year, month, version, varnm)
#             days = range(1, atm.days_this_month(year, month) + 1)
#             jdays = atm.season_days(atm.month_str(month), atm.isleap(year))
#             urls = [url_dict['%d%02d%02d' % (year, month, day)] for day in days]
#             func_kw = get_kw(jdays, latlon, plevs, sector_lons)
#             data = atm.load_concat(urls, varnm, concat_dim='day', func=calc_data,
#                                    func_kw=func_kw)
#             if isinstance(data, xray.DataArray):
#                 data = data.to_dataset()
#
#             # Save monthly files
#
#             for nm in data.data_vars:
#                 var = data[nm]
#                 filenm = get_filename(var, version, datadir, year, month)
#
#                 print('Saving to ' + filenm)
#                 ds = var.to_dataset()
#                 ds.to_netcdf(filenm, format=nc_fmt, engine=nc_eng)
#
#                 # Check if output is corrupted
#                 with xray.open_dataset(filenm) as ds_check:
#                     print(ds.dims.keys())
#                     print(ds.data_vars.keys())
#                     if len(ds.data_vars.keys()) > 1:
#                         raise ValueError('Corrupted monthly output file')
#
#         # Consolidate monthly files into yearly files
#         for nm in data.data_vars:
#             files = [get_filename(data[nm], version, datadir, year, m) for m in months]
#             var = atm.load_concat(files, nm, concat_dim='day')
#             filenm = get_filename(var, version, datadir, year)
#             print('Saving to ' + filenm)
#             var.to_dataset().to_netcdf(filenm)
#             print('Deleting monthly files:')
#             for filenm in files:
#                 print(filenm)
#                 os.remove(filenm)
