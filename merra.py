import numpy as np
import xray
import collections
import pandas as pd
import urllib2
from bs4 import BeautifulSoup

import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
import atmos as atm
from atmos import print_if

# ----------------------------------------------------------------------
def get_varname(var_id):
    """Return the variable name in MERRA data file.

    Parameters
    ----------
    var_id : {'u', 'v', 'omega', 'hgt', 'T', 'q', 'ps', 'evap', 'precip'}
    """

    var_dict = {'u' : 'U',
                'v' : 'V',
                'omega' : 'OMEGA',
                'hgt' : 'H',
                'T' : 'T',
                'q' : 'QV',
                'ps' : 'PS',
                'evap' : 'EVAP',
                'precip' : 'PRECTOT'}

    if var_id in var_dict:
        return var_dict[var_id]
    else:
        return var_id


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
def load_daily(year, month, var_id, concat_dim='TIME',
               subset1=(None, None, None), subset2=(None, None, None),
               verbose=True):
    """Return daily data for selected year and month for a single variable.

    Parameters
    ----------
    year, month : int
        Numeric year and month (1-12).
    var_id : {'u', 'v', 'omega', 'hgt', 'T', 'q', 'ps', 'evap', 'precip'},
             or str
        Variable ID.  Can be generic ID from the list above, in which
        case get_varname() is called to get the specific ID for MERRA. Or
        var_id can be the exact name as it appears in MERRA data files.
    concat_dim : str, optional
        Name of dimension for concatenation.
    subset1, subset2 : (str, float(s), float(s)), optional
        Tuple to indicate subset(s) to extract, in the form:
        (dim_name, lower_or_list, upper)
        e.g. subset1 = ('lon', 0, 120)
             subset2 = ('lat', -45, 45)
        e.g. subset1 = ('plev', 200, 200)
        The dimension name can be the actual dimension name
        (e.g. 'XDim') or a generic name (e.g. 'lon') and get_coord()
        is called to find the specific name.
    verbose : bool, optional
        If True, print updates while processing files.

    Returns
    -------
    data : xray.DataArray
        Daily data (3-hourly or hourly) for the month.
    """

    var = get_varname(var_id)
    date = '%d%02d' % (year, month)
    dataset = get_dataset(var_id, 'daily')
    urls = url_list(dataset)

    paths = [urls[key] for key in urls.keys() if date in key]

    data = atm.load_concat(paths, var, concat_dim, subset1, subset2, verbose)
    return data


# ----------------------------------------------------------------------
def monthly_from_daily(year, month, var_id, fluxes=False, fluxvars=('u', 'v'),
                       concat_dim='TIME', verbose=True):
    """Return the monthly mean of daily data.

    Parameters
    ----------
    year, month : int
        Numeric year and month (1-12).
    var_id : {'u', 'v', 'omega', 'hgt', 'T', 'q', 'ps', 'evap', 'precip'},
              or str
        Variable ID.  Can be generic ID from the list above, in which
        case get_varname() is called to get the specific ID for MERRA. Or
        var_id can be the exact name as it appears in MERRA data files.
    fluxes : bool, optional
        If True, return the mean fluxes u*var, v*var, along with var.
        If False, return only the specified variable.
    fluxvars : tuple of str, optional
        Names of u and v to extract from data to calculate fluxes. Only
        used if fluxes is True.
    concat_dim : str, optional
        Name of dimension for concatenation.
    verbose : bool, optional
        If True, print updates while processing files.

    Returns
    -------
    var_bar[, u_var_bar, v_var_bar] : xray.DataArray(s)
        Mean of daily data for variable var_id, and the mean of the
        daily zonal fluxes (u * var) and meridional fluxes (v * var),
        if applicable.

    Examples
    --------
    ubar = monthly_from_daily(1979, 1, 'u')
    qbar, uq_bar, vq_bar = monthly_from_daily(1979, 1, 'q', fluxes=True)
    """

    (u_nm, v_nm) = fluxvars

    # Extended variables (calculated from others)
    ext_vars = {'theta' : {'var_id0' : 'T',
                           'var_id1' : None},
                'theta_e' : {'var_id0' : 'T',
                             'var_id1' : 'q'}}

    if var_id in ext_vars:
        var_id0 = ext_vars[var_id.lower()]['var_id0']
        var_id1 = ext_vars[var_id.lower()]['var_id1']
    else:
        var_id0, var_id1 = var_id, None


    # Read metadata from one file to get pressure-level array
    dataset = get_dataset(var_id0, 'daily')
    if dataset.startswith('p_'):
        url = url_list(dataset, return_dict=False)[0]
        ds = xray.open_dataset(url)
        pname = atm.get_coord(ds, 'plev', 'name')
        plev = atm.get_coord(ds, 'plev')
        # Pressure levels in Pa for theta/theta_e calcs
        p_units = atm.pres_units(ds[pname].units)
        pres = atm.pres_convert(plev, p_units, 'Pa')
        ds.close()
        scale1, scale2 = 0.9999, 1.0001
    else:
        plev = [np.nan]
        pres = None
        if fluxes:
            raise ValueError('Fluxes cannot be calculated for surface data.')

    # Get daily data (raw or calculate extended variables)
    def get_data(var_id, var_id0, var_id1, pres, year, month, concat_dim,
                 subset1, verbose):
        var0 = load_daily(year, month, var_id0, concat_dim=concat_dim,
                          subset1=subset1, verbose=verbose)
        if var_id1 is not None:
            var1 = load_daily(year, month, var_id1, concat_dim=concat_dim,
                              subset1=subset1, verbose=verbose)
        if var_id == var_id0:
            var = var0
        elif var_id.lower() == 'theta':
            var = atm.potential_temp(var0, pres)
        elif var_id.lower() == 'theta_e':
            var = atm.equiv_potential_temp(var0, pres, var1)
        else:
            raise ValueError('Invalid var_id ' + var_id)
        return var

    # Iterate over vertical levels
    for k, p in enumerate(plev):
        if np.isnan(p):
            subset1 = (None, None, None)
            print_if('Surface data', verbose)
        else:
            subset1 = (pname, p * scale1, p * scale2)
            print_if('Pressure-level %.1f' % p, verbose)

        var = get_data(var_id, var_id0, var_id1, pres[k], year, month,
                       concat_dim, subset1, verbose)
        _, attrs, _ = atm.meta(var)

        if k == 0:
            var_bar = var.mean(dim=concat_dim)
            var_bar.attrs = attrs
        else:
            var_bar = xray.concat([var_bar, var.mean(dim=concat_dim)],
                                  dim=pname)

        if fluxes:
            u = load_daily(year, month, u_nm, concat_dim=concat_dim,
                           subset1=subset1, verbose=verbose)
            v = load_daily(year, month, v_nm, concat_dim=concat_dim,
                           subset1=subset1, verbose=verbose)
            u_var = u * var
            u_var.name = get_varname(u_nm) + '_*_' +  var_bar.name
            v_var = v * var
            v_var.name = get_varname(v_nm) + '_*_' +  var_bar.name
            if k == 0:
                u_var_bar = u_var.mean(dim=concat_dim)
                v_var_bar = v_var.mean(dim=concat_dim)
                units = var_bar.attrs['units'] + ' * ' + u.attrs['units']
                u_var_bar.attrs['units'] = units
                v_var_bar.attrs['units'] = units
            else:
                u_var_bar = xray.concat([u_var_bar, u_var.mean(dim=concat_dim)],
                                        dim=pname)
                v_var_bar = xray.concat([v_var_bar, v_var.mean(dim=concat_dim)],
                                        dim=pname)

    if fluxes:
        return var_bar, u_var_bar, v_var_bar
    else:
        return var_bar


# # ----------------------------------------------------------------------
# def monthly_from_daily(year, month, var_id_a, var_id_b=None, concat_dim='TIME',
#                        verbose=True):
#     """Return the monthly mean of daily data.
#
#     Parameters
#     ----------
#     year, month : int
#         Numeric year and month (1-12).
#     var_id_a, var_id_b : {'u', 'v', 'omega', 'hgt', 'T', 'q', 'ps',
#                           'evap', 'precip'},  or str
#         Variable IDs.  Can be generic ID from the list above, in which
#         case get_varname() is called to get the specific ID for MERRA. Or
#         var_id can be the exact name as it appears in MERRA data files.
#     concat_dim : str, optional
#         Name of dimension for concatenation.
#     verbose : bool, optional
#         If True, print updates while processing files.
#
#     Returns
#     -------
#     a_bar[, b_bar, ab_bar] : xray.DataArray(s)
#         Mean of daily data for variable a, variable b (if applicable),
#         and a * b (if applicable).
#
#     Examples
#     --------
#     ubar = monthly_from_daily(1979, 1, 'u')
#     ubar, vbar, uvbar = monthly_from_daily(1979, 1, 'u', 'v')
#     """
#
#     # Read metadata from one file to get pressure-level array
#     dataset = get_dataset(var_id_a, 'daily')
#     if dataset.startswith('p_'):
#         url = url_list(dataset, return_dict=False)[0]
#         ds = xray.open_dataset(url)
#         pname = atm.get_coord(ds, 'plev', 'name')
#         plev = atm.get_coord(ds, 'plev')
#         ds.close()
#         scale1, scale2 = 0.9999, 1.0001
#     else:
#         plev = [np.nan]
#
#     # Iterate over vertical levels
#     for k, p in enumerate(plev):
#         if np.isnan(p):
#             subset1 = (None, None, None)
#             print_if('Surface data', verbose)
#         else:
#             subset1 = (pname, p * scale1, p * scale2)
#             print_if('Pressure-level %.1f' % p, verbose)
#
#         a = load_daily(year, month, var_id_a, concat_dim=concat_dim,
#                        subset1=subset1, verbose=verbose)
#
#         if k == 0:
#             a_bar = a.mean(dim=concat_dim)
#         else:
#             a_bar = xray.concat([a_bar, a.mean(dim=concat_dim)], dim=pname)
#
#         if var_id_b is not None:
#             b = load_daily(year, month, var_id_b, concat_dim=concat_dim,
#                            subset1=subset1, verbose=verbose)
#             ab = a * b
#             ab.name = get_varname(var_id_a) + '_times_' + get_varname(var_id_b)
#             if k == 0:
#                 b_bar = b.mean(dim=concat_dim)
#                 ab_bar = ab.mean(dim=concat_dim)
#             else:
#                 b_bar = xray.concat([b_bar, b.mean(dim=concat_dim)], dim=pname)
#                 ab_bar = xray.concat([ab_bar, ab.mean(dim=concat_dim)],
#                                       dim=pname)
#
#     if var_id_b is None:
#         return a_bar
#     else:
#         return a_bar, b_bar, ab_bar
#

# ======================================================================
# Lists of OpenDAP urls for data files
# ======================================================================

# ----------------------------------------------------------------------
def url_list(dataset, return_dict=True):
    """Return list of OpenDAP urls for MERRA data files.

    This function reads a .csv file with a list of urls previously
    scraped from the MERRA website with merra_urls() and saved with
    save_urls().

    Parameters
    ----------
    dataset : {'p_monthly', 'p_daily', 'sfc_monthly', 'sfc_daily'} or string
        If dataset is in the list above, then csv file is assumed to be
        data/merra_urls_dataset.csv.  If not in the list above, then
        dataset should be the path of the csv file to read.
    return_dict : bool, optional
        If True, then return the urls in an OrderedDict of date:url
        pairs.  If False, then return urls as a list.

    Returns
    -------
    urls : OrderedDict() or list
        Dict of date:url pairs or list of urls.
    """

    if dataset in ['p_monthly', 'p_daily', 'sfc_monthly', 'sfc_daily']:
        filename = 'data/merra_urls_' + dataset + '.csv'
    else:
        filename = dataset

    df = pd.read_csv(filename, index_col=0)
    col_name = df.columns[0]
    dates = df.index.values
    dates = [str(d) for d in dates]
    files = df[col_name].values
    files = [str(f) for f in files]

    if return_dict:
        urls = collections.OrderedDict()
        for date, file in zip(dates, files):
            urls[date] = file
    else:
        urls = files

    return urls


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
def merra_urls(dataset='p_monthly'):
    """Return dict of OpenDAP urls for MERRA data.

    Parameters
    ----------
    dataset : {'p_monthly', 'p_daily', 'sfc_monthly', 'sfc_daily'}

    Returns
    -------
    urls : dict of date:url for each date in the dataset
    """

    dataset = dataset.lower()
    years = {}
    for y in range(1979,2015):
        years[y] = '%d' % y
    months = {}
    for m in range(1,13):
        months[m] = '%02d' % m

    rootdir_3d = 'http://goldsmr3.sci.gsfc.nasa.gov/opendap/'
    rootdir_2d = 'http://goldsmr2.sci.gsfc.nasa.gov/opendap/'

    kinds = []
    basedirs = {}
    urls = {}

    # Info for each dataset
    # ---------------------
    # Monthly pressure-level data
    kind = 'p_monthly'
    kinds.append(kind)
    basedirs[kind] = rootdir_3d + 'MERRA_MONTHLY/MAIMCPASM.5.2.0/'

    # Monthly surface fluxes
    kind = 'sfc_monthly'
    kinds.append(kind)
    basedirs[kind] = rootdir_2d + 'MERRA_MONTHLY/MATMNXFLX.5.2.0/'

    # Daily pressure-level data (3-hourly)
    kind = 'p_daily'
    kinds.append(kind)
    basedirs[kind] = rootdir_3d + 'MERRA/MAI3CPASM.5.2.0/'

    # Daily surface fluxes (hourly)
    kind = 'sfc_daily'
    kinds.append(kind)
    basedirs[kind] = rootdir_2d + 'MERRA/MAT1NXFLX.5.2.0/'

    # Helper function to make monthly urls
    def monthly_urls(basedir):
        url_dict = collections.OrderedDict()
        for y in years:
            print(years[y])
            dirname = basedir + years[y] + '/'
            files = scrape_url(dirname + 'contents.html')
            files.sort()
            dates = [extract_date(nm, width=6) for nm in files]
            for date, nm in zip(dates, files):
                url_dict[date] = dirname + nm
        return url_dict

    # Helper function to make daily urls
    def daily_urls(basedir):
        url_dict = collections.OrderedDict()
        for y in years:
            for m in months:
                print(years[y] + months[m])
                dirname = basedir + years[y] + '/' + months[m] + '/'
                files = scrape_url(dirname + 'contents.html')
                files.sort()
                dates = [extract_date(nm, width=8) for nm in files]
                for date, nm in zip(dates, files):
                    url_dict[date] = dirname + nm
        return url_dict

    # Extract urls
    if dataset in ['p_monthly', 'sfc_monthly']:
        urls = monthly_urls(basedirs[dataset])
    elif dataset in ['p_daily', 'sfc_daily']:
        urls = daily_urls(basedirs[dataset])
    else:
        raise ValueError('Invalid dataset ' + dataset)

    return urls


# ----------------------------------------------------------------------
def save_urls(filestart='data/merra_urls_',
              keys=['p_monthly', 'sfc_monthly', 'p_daily', 'sfc_daily']):
    """Save lists of MERRA urls."""

    for key in keys:
        print('*********************\n' + key)
        filename = filestart + key + '.csv'
        urls = merra_urls(key)
        atm.print_odict(urls)

        dates, files = [], []
        for date in urls:
            dates.append(date)
            files.append(urls[date])

        df = pd.DataFrame(files, index=pd.Series(dates, name='date'),
                          columns=['url'])

        print('Writing urls to ' + filename)
        df.to_csv(filename)
