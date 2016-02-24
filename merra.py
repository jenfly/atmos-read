"""
Functions to read and save MERRA reanalysis data.

Convention for function names
  Starts with - read_ : Read from OpenDAP url(s)
              - load_ : Load from locally saved files
"""

from __future__ import division
import numpy as np
import xray
import collections
import os
import pandas as pd
import urllib2
from bs4 import BeautifulSoup

import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
import atmos as atm
from atmos import print_if

# ----------------------------------------------------------------------
def get_varname(var_id):
    """Return the variable name in MERRA naming convention.

    Parameters
    ----------
    var_id : {'u', 'v', 'omega', 'hgt', 'T', 'q', 'ps', 'evap', 'precip'}
    """

    var_dict = {'u' : 'U', 'v' : 'V', 'omega' : 'OMEGA', 'hgt' : 'H',
                'T' : 'T', 'q' : 'QV', 'ps' : 'PS', 'evap' : 'EVAP',
                'precip' : 'PRECTOT'}

    if var_id in var_dict:
        return var_dict[var_id]
    else:
        return var_id


# ----------------------------------------------------------------------
def get_url_opts(var_id, version='merra'):

    varnm = get_varname(var_id)

    p_res = {'merra' : 'C', 'merra2' : 'N'}[version]

    optlist = {'int_t' : ('X', 'N', 'T', 'INT'),
               'int_i' : ('X', 'N', 'I', 'INT'),
               'flx_t' : ('X', 'N', 'T', 'FLX'),
               'slv_t' : ('X', 'N', 'T', 'SLV'),
               'asm_i' : ('P', p_res, 'I', 'ASM')}

    optkeys = {}
    for nm in ['UFLXQV', 'VFLXQV', 'UFLXCPT', 'VFLXCPT', 'UFLXPHI', 'VFLXPHI',
               'DQVDT_ANA']:
        optkeys[nm] = 'int_t'
    for nm in ['TQV']:
        optkeys[nm] = 'int_i'
    for nm in ['PRECTOT', 'EVAP', 'EFLUX', 'HFLUX', 'ULML', 'VLML', 'QLML',
               'TLML', 'HLML']:
        optkeys[nm] = 'flx_t'
    for nm in ['PS', 'SLP']:
        optkeys[nm] = 'slv_t'
    for nm in ['U', 'V', 'OMEGA', 'T', 'QV', 'H']:
        optkeys[nm] = 'asm_i'

    vertical, res, time_kind, kind = optlist[optkeys[varnm]]
    opts = {'vertical' : vertical, 'res' : res, 'time_kind' : time_kind,
            'kind' : kind}

    return opts


# ----------------------------------------------------------------------
def get_dataset(var_id, time_res='daily', default='p'):
    """Return the dataset ID corresponding to the variable.

    Parameters
    ----------
    var_id : str
        Variable name.  Can be a generic ID as input to get_varname(),
        or a specific name from MERRA data files.
    time_res : {'daily', 'monthly'}
        Time resolution of dataset.
    default : {'p', 'sfc'}
        If the variable is in both pressure-level and surface flux
        data, then default to this dataset type.

    Returns
    -------
    dataset : {'p_monthly', 'p_daily', 'sfc_monthly', 'sfc_daily'}
        Name of the dataset containing the variable (pressure-level
        or surface fluxes), at the specified time resolution.
    """

    var = get_varname(var_id)

    p_vars = [u'SLP', u'PS', u'PHIS', u'H', u'O3', u'QV', u'QL', u'QI', u'RH',
              u'T', u'U', u'V', u'EPV', u'OMEGA', u'Cov_U_V', u'Cov_U_T',
              u'Cov_V_T', u'Cov_U_H', u'Cov_V_H', u'Cov_U_QV', u'Cov_V_QV',
              u'Cov_U_QL', u'Cov_V_QL', u'Cov_U_QI', u'Cov_V_QI', u'Cov_U_EPV',
              u'Cov_V_EPV', u'Cov_U_O3', u'Cov_V_O3', u'Cov_OMEGA_U',
              u'Cov_OMEGA_V', u'Cov_OMEGA_T', u'Cov_OMEGA_QV', u'Cov_OMEGA_QL',
              u'Cov_OMEGA_QI', u'Cov_OMEGA_O3', u'vsts', u'Var_SLP', u'Var_PS',
              u'Var_PHIS', u'Var_H', u'Var_O3', u'Var_QV', u'Var_QL', u'Var_QI',
              u'Var_RH', u'Var_T', u'Var_U', u'Var_V', u'Var_EPV', u'Var_OMEGA']

    sfc_vars = [u'EFLUX', u'EVAP', u'HFLUX', u'TAUX', u'TAUY', u'TAUGWX',
                u'TAUGWY', u'PBLH', u'DISPH', u'BSTAR', u'USTAR', u'TSTAR',
                u'QSTAR', u'RI', u'Z0H', u'Z0M', u'HLML', u'TLML', u'QLML',
                u'ULML', u'VLML', u'RHOA', u'SPEED', u'CDH', u'CDQ', u'CDM',
                u'CN', u'TSH', u'QSH', u'FRSEAICE', u'PRECANV', u'PRECCON',
                u'PRECLSC', u'PRECSNO', u'PRECTOT', u'PGENTOT', u'Var_EFLUX',
                u'Var_EVAP', u'Var_HFLUX', u'Var_TAUX', u'Var_TAUY',
                u'Var_TAUGWX', u'Var_TAUGWY', u'Var_PBLH', u'Var_DISPH',
                u'Var_BSTAR', u'Var_USTAR', u'Var_TSTAR', u'Var_QSTAR',
                u'Var_RI', u'Var_Z0H', u'Var_Z0M', u'Var_HLML', u'Var_TLML',
                u'Var_QLML', u'Var_ULML', u'Var_VLML', u'Var_RHOA',
                u'Var_SPEED', u'Var_CDH', u'Var_CDQ', u'Var_CDM', u'Var_CN',
                u'Var_TSH', u'Var_QSH', u'Var_PRECANV', u'Var_PRECCON',
                u'Var_PRECLSC', u'Var_PRECSNO', u'Var_PRECTOT', u'Var_PGENTOT']

    if var in p_vars and var in sfc_vars:
        dataset = default + '_' + time_res
    elif var in p_vars:
        dataset = 'p_' + time_res
    elif var in sfc_vars:
        dataset = 'sfc_' + time_res
    else:
        raise ValueError('var_id ' + var_id + ' not found.')

    return dataset


# ----------------------------------------------------------------------
def read_daily(var_ids, year, month, days=None, concat_dim='TIME',
               subset_dict=None, verbose=True):
    """Return MERRA daily pressure-level data for selected variable(s).

    Reads daily MERRA data from OpenDAP urls and concatenates into a
    single DataArray or Dataset for the selected days of the month.

    Parameters
    ----------
    var_ids : str or list of str
        Variable ID(s).  Can be generic ID from the list below, in which
        case get_varname() is called to get the specific ID for MERRA. Or
        var_id can be the exact name as it appears in MERRA data files.
        Generic IDs:
          {'u', 'v', 'omega', 'hgt', 'T', 'q', 'ps', 'evap', 'precip'}
    year, month : int
        Numeric year and month (1-12).
    days : list of ints, optional
       Subset of days to read. If None, all days are included.
    concat_dim : str, optional
        Name of dimension for concatenation.
    subset_dict : dict of 2-tuples, optional
        Dimensions and subsets to extract.  Each entry in subset_dict
        is in the form {dim_name : (lower_or_list, upper)}, where:
        - dim_name : string
            Name of dimension to extract from.
            The dimension name can be the actual dimension name
            (e.g. 'XDim') or a generic name (e.g. 'lon') and get_coord()
            is called to find the specific name.
        - lower_or_list : scalar or list of int or float
            If scalar, then used as the lower bound for the   subset range.
            If list, then the subset matching the list will be extracted.
        - upper : int, float, or None
            Upper bound for subset range. If lower_or_list is a list,
            then upper is ignored and should be set to None.
    verbose : bool, optional
        If True, print updates while processing files.

    Returns
    -------
    data : xray.DataArray or xray.Dataset
        Daily data (3-hourly or hourly) for the month or a selected
        subset of days.
    """

    var_ids = atm.makelist(var_ids)
    var_nms = [get_varname(var_id) for var_id in var_ids]
    dataset = get_dataset(var_ids[0], 'daily')
    urls = url_list(dataset)

    if days is None:
        # All days in the month
        dates = ['%d%02d' % (year, month)]
    elif isinstance(days, int):
        # Single day
        dates = ['%d%02d%02d' % (year, month, days)]
    else:
        # Subset of days
        dates = ['%d%02d%02d' % (year, month, d) for d in days]

    paths = []
    for date in dates:
        paths.extend([urls[key] for key in urls.keys() if date in key])

    data = atm.load_concat(paths, var_nms, concat_dim, subset_dict,
                           verbose)
    return data


# ----------------------------------------------------------------------
def read_daily_eta(var_id, level, year, month, days=None, concat_dim='TIME',
                   xsub='[330:2:450]', ysub='[60:2:301]', verbose=True):
    """Return MERRA daily eta-level data for a single variable.

    Reads a single eta level of daily MERRA data from OpenDAP urls and
    concatenates into a DataArray for the selected days of the month.

    Parameters
    ----------
    var_id : str
        Variable ID.  Can be generic ID from the list below, in which
        case get_varname() is called to get the specific ID for MERRA. Or
        var_id can be the exact name as it appears in MERRA data files.
        Generic IDs:
          {'u', 'v', 'omega', 'hgt', 'T', 'q', 'ps', 'evap', 'precip'}
    level : int
        Eta level to extract (0-71).  Level 71 is near-surface and level 0
        is the top of atmosphere.
    year, month : int
        Numeric year and month (1-12).
    days : list of ints, optional
       Subset of days to read. If None, all days are included.
    concat_dim : str, optional
        Name of dimension for concatenation.
    xsub, ysub : str, optional
        Indices of longitude and latitude subsets to extract.
    verbose : bool, optional
        If True, print updates while processing files.

    Returns
    -------
    data : xray.DataArray or xray.Dataset
        Daily data (3-hourly or hourly) for the month or a selected
        subset of days.
    """

    varnm = get_varname(var_id)
    tsub = '[0:1:3]'
    zsub = '[%d:1:%d]' % (level, level)

    def datafile(year, mon, day, varnm, xsub, ysub, zsub, tsub):
        basedir = ('http://goldsmr3.sci.gsfc.nasa.gov:80/opendap/MERRA/'
                   'MAI6NVANA.5.2.0/')
        url = ('%s%d/%02d/MERRA100.prod.assim.inst6_3d_ana_Nv.%d%02d%02d.hdf'
               '?%s%s%s%s%s,XDim%s,YDim%s,Height%s,TIME%s') % (basedir, year,
               mon, year, mon, day, varnm, tsub, zsub, ysub, xsub, xsub, ysub,
               zsub, tsub)
        return url

    if days is None:
        days = range(1, atm.days_this_month(year, month) + 1)
    urls = [datafile(year, month, day, varnm, xsub, ysub, zsub, tsub) for day
            in atm.makelist(days)]

    var = atm.load_concat(urls, varnm, concat_dim, verbose=verbose)

    return var


# ----------------------------------------------------------------------
def load_daily_season(pathstr, year, season='ann', var_ids=None,
                      lat1=-90, lat2=90, lon1=0, lon2=360,
                      verbose=True, concat_dim=None):
    """Return daily data for a selected year, season and lat-lon subset.

    Loads daily data from locally saved files and concatenates it into
    a single DataArray or Dataset for that year and season.

    Parameters
    ----------
    pathstr : str
       Beginning of path for each data file, where each file name is in
       the format *yyyymm.nc.
       e.g. pathstr = '~/datastore/merra/daily/u200_'
    year : int
       Year to load.
    season : str, optional
       Season to load. Valid values are as listed in atm.season_months()
       e.g. 'jul', 'jja', 'ann'
       Default is entire year ('ann')
    var_ids : str or list of str, optional
       Variable(s) to extract. If omitted, all variables in the data are
       included and the output is a Dataset.
    lat1, lat2, lon1, lon2 : floats, optional
        Lat-lon subset to extract.
    concat_dim : str, optional
        Name of time dimension for concatenation. If None, then
        atm.get_coord() is called to get the name from the data file.
    verbose : bool, optional
        If True, print updates while processing files.

    Returns
    -------
    data : xray.DataArray or xray.Dataset
    """

    months = atm.season_months(season)
    paths = []
    for m in months:
        datestr = '%d%02d' % (year, m)
        paths.append(pathstr + datestr + '.nc')

    # Make sure longitude range is consistent with data
    with xray.open_dataset(paths[0]) as ds:
        lonmax = atm.lon_convention(atm.get_coord(ds, 'lon'))
        if concat_dim is None:
            concat_dim = atm.get_coord(ds, 'time', 'name')
    if lon2 - lon1 == 360:
        if lonmax < lon2:
            offset = -180
        elif lonmax > lon2:
            offset = 180
        else:
            offset = 0
        lon1, lon2 = lon1 + offset, lon2 + offset
    print(lon1, lon2, lonmax)

    # Load daily data
    if var_ids is None:
        var_nms = None
    else:
        var_nms = [get_varname(var_id) for var_id in atm.makelist(var_ids)]
    subset_dict = {'lat' : (lat1, lat2), 'lon' : (lon1, lon2)}
    data = atm.load_concat(paths, var_nms, concat_dim, subset_dict, verbose)

    return data


# ----------------------------------------------------------------------
def calc_fluxes(year, month,
                var_ids=['u', 'q', 'T', 'theta', 'theta_e', 'hgt'],
                concat_dim='TIME', scratchdir=None, keepscratch=False,
                verbose=True):
    """Return the monthly mean of MERRA daily fluxes.

    Reads MERRA daily data from OpenDAP urls, computes fluxes, and
    returns the monthly mean of the daily variable and its zonal and
    meridional fluxes.

    Parameters
    ----------
    year, month : int
        Numeric year and month (1-12).
    var_ids : list of str, optional
        IDs of variables to include.
    concat_dim : str, optional
        Name of dimension for concatenation.
    scratchdir : str, optional
        Directory path to store temporary files while processing data.
        If omitted, the current working directory is used.
    keepscratch : bool, optional
        If True, scratch files are kept in scratchdir. Otherwise they
        are deleted.
    verbose : bool, optional
        If True, print updates while processing files.

    Returns
    -------
    data : xray.Dataset
        Mean of daily data and the mean of the daily zonal fluxes
        (u * var) and meridional fluxes (v * var), for each variable
        in var_ids.
    """

    nms = [get_varname(nm) for nm in atm.makelist(var_ids)]
    u_nm, v_nm = get_varname('u'), get_varname('v')
    nms.extend([u_nm, v_nm])
    if 'theta' in nms:
        nms.append(get_varname('T'))
    if 'theta_e' in nms:
        nms.extend([get_varname('T'), get_varname('q')])
    nms = set(nms)

    days = range(1, atm.days_this_month(year, month) + 1)

    def scratchfile(nm, k, year, month, day):
        filestr = '%s_level%d_%d%02d%02d.nc' % (nm, k, year, month, day)
        if scratchdir is not None:
            filestr = scratchdir + '/' + filestr
        return filestr

    # Read metadata from one file to get pressure-level array
    dataset = 'p_daily'
    url = url_list(dataset, return_dict=False)[0]
    with xray.open_dataset(url) as ds:
        pname = atm.get_coord(ds, 'plev', 'name')
        plev = atm.get_coord(ds, 'plev')
        # Pressure levels in Pa for theta/theta_e calcs
        p_units = atm.pres_units(ds[pname].units)
        pres = atm.pres_convert(plev, p_units, 'Pa')

    # Get daily data (raw and calculate extended variables)
    def get_data(nms, pres, year, month, day, concat_dim, subset_dict, verbose):
        # Lists of raw and extended variables
        ids = list(nms)
        ext = []
        for var in ['theta', 'theta_e']:
            if var in ids:
                ext.append(var)
                ids.remove(var)

        # Read raw data and calculate extended variables
        data = read_daily(ids, year, month, day, concat_dim=concat_dim,
                          subset_dict=subset_dict, verbose=verbose)
        if 'theta' in ext:
            print_if('Computing potential temperature', verbose)
            T = data[get_varname('T')]
            data['theta'] = atm.potential_temp(T, pres)
        if 'theta_e' in ext:
            print_if('Computing equivalent potential temperature', verbose)
            T = data[get_varname('T')]
            q = data[get_varname('q')]
            data['theta_e'] = atm.equiv_potential_temp(T, pres, q)

        return data

    # Iterate over vertical levels
    for k, p in enumerate(plev):
        subset_dict = {pname : (p, p)}
        print_if('Pressure-level %.1f' % p, verbose)

        files = []

        for day in days:
            # Read data for this level and day
            ds = get_data(nms, pres[k], year, month, day, concat_dim,
                           subset_dict, verbose)

            # Compute fluxes
            print_if('Computing fluxes', verbose)
            u = ds[get_varname('u')]
            v = ds[get_varname('v')]
            for nm in var_ids:
                var = ds[get_varname(nm)]
                varname, attrs, _, _ = atm.meta(var)
                u_var = u * var
                v_var = v * var

                u_var.name = get_varname(u_nm) + '*' +  var.name
                units = var.attrs['units'] + ' * ' + u.attrs['units']
                u_var.attrs['units'] = units
                v_var.name = get_varname(v_nm) + '*' +  var.name
                v_var.attrs['units'] = units
                ds[u_var.name] = u_var
                ds[v_var.name] = v_var

            # Save to temporary scratch file
            filenm = scratchfile('fluxes', k, year, month, day)
            files.append(filenm)
            print_if('Saving to scratch file ' + filenm, verbose)
            ds.to_netcdf(filenm)

        # Concatenate daily scratch files
        ds = atm.load_concat(files)

        if not keepscratch:
            for f in files:
                os.remove(f)

        # Compute monthly means
        print_if('Computing monthly means', verbose)
        if k == 0:
            data = ds.mean(dim=concat_dim)
        else:
            data = xray.concat([data, ds.mean(dim=concat_dim)], dim=pname)

    for var in data.data_vars:
        data[var].attrs = ds[var].attrs

    return data


# ======================================================================
# Lists of OpenDAP urls for data files
# ======================================================================

# ----------------------------------------------------------------------
def scrape_url(url, ending='.hdf.html', cut='.html'):
    """Scrape url for links with specified ending string."""

    soup = BeautifulSoup(urllib2.urlopen(url))
    links = []
    for link in soup.find_all('a'):
        links.append(link.get('href'))

    links = list(set([s for s in links if s.endswith(ending)]))
    links = [s.split(cut)[0] for s in links]
    return links


# ----------------------------------------------------------------------
def extract_date(filename, width, ending='.hdf'):
    """Extract yyyymmdd or yyyymm from file name."""
    s = filename.split(ending)[0]
    date = s[-width:]
    return date


# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
def get_urls(years, months=None, version='merra', vertical='P', res='C',
               time_kind='I', kind='ASM'):
    """Return dict of OpenDAP urls for MERRA and MERRA-2 daily data.

    Parameters
    ----------
    years : list or np.ndarray
        List of years to extract urls for.
    months: list or np.ndarray, optional
        List of months to extract urls for.  If None, then all months (1-12)
        are extracted.
    version : {'merra', 'merra2'}, optional
        Select MERRA or MERRA-2 data.
    vertical : {'P', 'X', 'V', 'E'}, optional
        Vertical location : on pressure levels (P), 2-D (X), model
        layers (V), or model layer edges (E).
    res : {'N', 'C'}, optional
        Horizontal resolution: native (N) or coarse (C).
    time_kind : {'I', 'T'}, optional
        Instantaneous (I) or time-averaged (T) diagnostics.
    kind : {'ASM', 'SLV', 'FLX', 'RAD'}, optional
        Type of dataset: assimilated 3-d (ASM), atmospheric single-level
        (SLV), surface turbulent fluxes (FLX), or surface and TOA
        radiation fluxes (RAD).

    Returns
    -------
    urls : dict of date:url for each date in the dataset
    """

    # Housekeeping
    version = version.lower()
    time_kind = time_kind.upper()
    res = res.upper()
    vertical = vertical.upper()
    kind = kind.upper()

    # Make dicts of years and months
    yearvals = atm.makelist(years)
    years = {y : '%d' % y for y in yearvals}
    if months is None:
        monthvals = range(1, 13)
    else:
        monthvals = atm.makelist(months)
    months = {m : '%02d' % m for m in monthvals}

    urlstr = 'http://goldsmr%d.sci.gsfc.nasa.gov/opendap/%s'
    servers = {'merra_X' : urlstr % (2, 'MERRA/MA'),
               'merra' : urlstr % (3, 'MERRA/MA'),
               'merra2_X' : urlstr % (4, 'MERRA2/M2'),
               'merra2' : urlstr % (5, 'MERRA2/M2')}
    version_num = {'merra' : '.5.2.0/', 'merra2' : '.5.12.4/'}
    fmt = {'merra' : '.hdf', 'merra2' : '.nc4'}

    if vertical == 'X':
        time_res = '1'
        server_key = version + '_X'
    else:
        time_res = '3'
        server_key = version

    try:
        basedir = servers[server_key]
        vnum = version_num[version]
    except KeyError:
        raise ValueError('Invalid version %s.  Options are: merra, merra2.' %
                         version)

    basedir = basedir + time_kind + time_res + res + vertical + kind + vnum
    print('Scraping filenames from ' + basedir)

    # Helper function to make daily urls
    def daily_urls(basedir, years, months, fmt):
        url_dict = collections.OrderedDict()
        for y in years:
            for m in months:
                print(years[y] + months[m])
                dirname = basedir + years[y] + '/' + months[m] + '/'
                files = scrape_url(dirname + 'contents.html',
                                   ending=fmt + '.html')
                files.sort()
                dates = [extract_date(nm, width=8, ending=fmt) for nm in files]
                for date, nm in zip(dates, files):
                    url_dict[date] = dirname + nm
        return url_dict

    # Extract urls
    urls = daily_urls(basedir, years, months, fmt[version])

    return urls
