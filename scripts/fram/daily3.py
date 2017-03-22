"""
3-D variables:
--------------
Instantaneous:
['U', 'V', 'OMEGA', 'T', 'QV', 'H']

Time-average:
['DUDTANA']

2-D variables:
--------------
Time-average surface fluxes:
['PRECTOT', 'EVAP', 'EFLUX', 'HFLUX', 'QLML', 'TLML']

Time-average vertically integrated fluxes:
['UFLXQV', 'VFLXQV', 'VFLXCPT', 'VFLXPHI']

Instantaneous vertically integrated fluxes:
['TQV']

Single-level atmospheric variables:
['PS', 'SLP']
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

version = 'merra'
years = [2000]
#version = 'merra2'
#years = np.arange(1980, 2016)

datadir = atm.homedir() + 'datastore/' + version + '/daily/'
months = range(7, 13)

varnms = ['U', 'V', 'OMEGA']

latlon=(-90, 90, 40, 120)
plevs = [1000,925,850,775,700,600,500,400,300,250,200,150,100,70,50,30,20]
#plevs=(850, 200)
sector_lons=(60, 100)
#dp_vars = ['U', 'OMEGA']
dp_vars = []


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

def sector_mean(var, lon1, lon2):
    """Return the sector mean of a variable."""
    name = var.name
    lonstr = atm.latlon_str(lon1, lon2, 'lon')
    if (lon2 - lon1) == 360:
        lon1, lon2 = None, None
        name_out = name + '_ZON'
    else:
        name_out = name + '_SEC'
    varbar = atm.dim_mean(var, 'lon', lon1, lon2)
    varbar.name = name_out
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
        atm.disptime()
        data.to_netcdf(filenm, **nc_kw)
        files[nm] = filenm
    return files

def read_url(url, varnms, datadir, year, month, day, jday, calc_kw, nc_kw):
    """Open url and process selected variables."""
    # Number of times to attempt opening url (in case of server problems)
    NMAX = 3
    # Wait time (seconds) between attempts
    WAIT = 5

    atm.disptime()
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
nc_kw = { 'merra2' : {'format' : 'NETCDF4_classic', 'engine' : 'netcdf4'},
          'merra' : {'format' : None, 'engine' : None}}[version]

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
    # for nm in dailyfiles:
    #     data = atm.load_concat(dailyfiles[nm], concat_dim='day')
    #     for varnm in data.data_vars:
    #         var = data[varnm]
    #         filenm = get_filename(var, version, datadir, year)
    #         var.name = var.attrs.get('varnm', varnm)
    #         print('Saving to ' + filenm)
    #         atm.save_nc(filenm, var)
    #     print('Deleting daily files')
    #     for filenm in dailyfiles[nm]:
    #         print(filenm)
    #         os.remove(filenm)
