import numpy as np
import matplotlib.pyplot as plt
import xray
import collections
import pandas as pd
import urllib2
from bs4 import BeautifulSoup

import atmos as atm




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
