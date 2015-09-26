import numpy as np
import matplotlib.pyplot as plt
import xray
import collections
import urllib2
from bs4 import BeautifulSoup

import atmos as atm


def scrape_url(url, ending='.hdf.html', cut='.html'):
    """Scrape url for links with specified ending string."""

    soup = BeautifulSoup(urllib2.urlopen(url))
    links = []
    for link in soup.find_all('a'):
        links.append(link.get('href'))

    links = list(set([s for s in links if s.endswith(ending)]))
    links = [s.split(cut)[0] for s in links]
    return links


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
    dataset : {'p_monthly', 'p_daily', 'sfc_monthly', sfc_daily}

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
            dates = [extract_date(nm, width=6) for nm in sorted(files)]
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
                dates = [extract_date(nm, width=8) for nm in sorted(files)]
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

# Make this a new function:
urls_all = {}
keys = ['p_monthly', 'sfc_monthly', 'p_daily', 'sfc_daily']
for key in keys:
    print('*********************')
    print(key)
    urls_all[key] = merra_urls(key)
    atm.print_odict(urls_all[key])

# Make this a new function:
for key in keys:
    urls = urls_all[key]
    filename = 'merra_urls_' + key + '.csv'
    print(filename)
    with open(filename, 'w') as f:
        for date in urls:
            f.write(date + ', ' + urls[date] + '\n')



# ----------------------------------------------------------------------

# Test the urls that were written to file
filename = 'merra_urls_p_daily.csv'
urls_in = []
with open(filename, 'rU') as f:
    for line in f:
        urls_in.append(line.split(', ')[1].replace('\n',''))

url = urls_in[100]
ds = xray.open_dataset(url)
print(ds)
# ----------------------------------------------------------------------

url = ('http://goldsmr3.sci.gsfc.nasa.gov/opendap/MERRA_MONTHLY/'
    'MAIMCPASM.5.2.0/1979/MERRA100.prod.assim.instM_3d_asm_Cp.197907.hdf')

url2 = ('http://goldsmr2.sci.gsfc.nasa.gov/opendap/MERRA_MONTHLY/'
    'MATMNXFLX.5.2.0/1979/MERRA100.prod.assim.tavgM_2d_flx_Nx.197907.hdf')

ds = xray.open_dataset(url)
ds2 = xray.open_dataset(url2)

u = ds['U']
v = ds['V']
q = ds['QV']
lat = get_coord(u, 'lat')
lon = get_coord(u, 'lon')
plev = get_coord(u, 'plev')

# Convert from (kg/m^2)/s to mm/day
SCALE = 60 * 60 * 24

precip = ds2['PRECTOT'] * SCALE
evap = ds2['EVAP'] * SCALE

mfc = av.moisture_flux_conv(u*q, v*q)

lon1, lon2 = 0, 150
lat1, lat2 = 0, 50
axlims = (lat1, lat2, lon1, lon2)

plt.figure(figsize=(7,8))
plt.subplot(211)
ap.pcolor_latlon(precip, cmap='hot_r', axlims=axlims)
plt.title('Total precip')
plt.subplot(212)
ap.pcolor_latlon(evap, cmap='hot_r', axlims=axlims)
plt.title('Evap')

plt.figure(figsize=(7,8))
plt.subplot(211)
ap.pcolor_latlon(precip - evap, axlims=axlims)
plt.title('Net precip')
plt.subplot(212)
ap.pcolor_latlon(mfc, axlims=axlims)
plt.title('MFC')

# ----------------------------------------------------------------------
# Sub-daily data (3-hourly and hourly for surface fluxes)

# Hourly precip
# url3 = ('http://goldsmr2.sci.gsfc.nasa.gov/opendap/MERRA/MAT1NXFLX.5.2.0/'
#     '1979/07/MERRA100.prod.assim.tavg1_2d_flx_Nx.19790701.hdf')
# ds3 = xray.open_dataset(url3)
#
# precip = ds3['PRECTOT']
