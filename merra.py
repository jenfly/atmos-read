import numpy as np
import xray
import collections
import pandas as pd
import urllib2
from bs4 import BeautifulSoup

import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
import atmos as atm

# ----------------------------------------------------------------------
def varname(var_id):
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
        return var


# ----------------------------------------------------------------------
def get_dataset(var_id, time_res='daily', default='p'):
    """Return the dataset ID corresponding to the variable.

    Parameters
    ----------
    var_id : str
        Variable name.  Can be a generic ID as input to varname(), or
        a specific name from MERRA data files.
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

    var = varname(var_id)

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
        case varname() is called to get the specific ID for MERRA. Or
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

    #var = varname(var_id)
    date = '%d%02d' % (year, month)

    # if var_id in ['evap', 'precip']:
    #     dataset = 'sfc_daily'
    # else:
    #     dataset = 'p_daily'
    dataset = get_dataset(var, )
    urls = url_list(dataset)

    paths = [urls[key] for key in urls.keys() if date in key]

    data = atm.load_concat(paths, var, concat_dim, subset1, subset2, verbose)
    return data


# ----------------------------------------------------------------------
def monthly_from_daily(year, month, var_id_a, var_id_b=None, concat_dim='TIME',
                       verbose=True):
    """Return the monthly mean of daily data.

    """

# ======================================================================
# Lists of OpenDAP urls for data files
# ======================================================================

# ----------------------------------------------------------------------
def url_list(dataset):
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

    Returns
    -------
    urls : OrderedDict()
        Dict of date:url pairs.
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

    urls = collections.OrderedDict()
    for date, file in zip(dates, files):
        urls[date] = file

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
