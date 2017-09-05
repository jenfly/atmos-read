import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')
sys.path.append('/home/jwalker/dynamics/python/monsoon-onset')

import xarray as xray
import numpy as np
import collections
import pandas as pd
import scipy

import atmos as atm

datakind = 'INT'
savefile = 'scripts/merra_urls/merra2_' + datakind.lower() + '_%d.txt'
years = np.arange(1980, 2016)

def prodnum(year):
    prod_dict = {yr: '100' for yr in range(1980, 1992)}
    for yr in range(1992, 2001):
        prod_dict[yr] = '200'
    for yr in range(2001, 2011):
        prod_dict[yr] = '300'
    for yr in range(2011, 2017):
        prod_dict[yr] = '400'
    prod = prod_dict[year]
    return prod


def varlist(datakind):
    vardict = {'INT' : ['UFLXCPT', 'UFLXPHI', 'UFLXQV',  'VFLXCPT', 'VFLXPHI', 'VFLXQV' ],
               'RAD' : ['LWTUP', 'SWGNT', 'LWGNT', 'SWTNT', 'HFLUX', 'EFLUX']}
    return vardict.get(datakind.upper())


def varstr(vlist):
    return '&VARIABLES=' + '%2C'.join(map(str.lower, vlist))


def geostr(lat1, lat2, lon1, lon2):
    sep = '%2C'
    box = ('&FORMAT=bmM0Lw&BBOX=' + str(lat1) + sep + str(lon1) + sep + str(lat2)
           + sep + str(lon2))
    return box


def get_url(year, mon, day, datakind='int', lat1=-90, lat2=90, lon1=-180, lon2=180):

    yrstr = '%d' % year
    mstr = '%02d' % mon
    dstr = '%02d' % day
    datestr = yrstr + mstr + dstr
    prod = prodnum(year)

    shortname = 'M2T1NX' + datakind.upper()
    longname = 'MERRA2_' + prod + '.tavg1_2d_' + datakind.lower() + '_Nx.' + datestr
    variables = varstr(varlist(datakind))
    box = geostr(lat1, lat2, lon1, lon2)

    base_url = ('http://goldsmr4.gesdisc.eosdis.nasa.gov/daac-bin/OTF/'
                'HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2F' + shortname +
                '.5.12.4%2F')

    label = '&LABEL=' + longname + '.SUB.nc4'

    url = (base_url + yrstr + '%2F' + mstr + '%2F' + longname + '.nc4' + box + label
           + '&FLAGS=1&SHORTNAME=' + shortname + '&SERVICE=SUBSET_MERRA2'
           +'&LAYERS=&VERSION=1.02' + variables)

    return url


def writefile(url_list, filenm):
    print('Saving url list to ' + filenm)
    with open(filenm, 'w') as f:
        for url in url_list:
            f.write(url + '\n')


for year in years:
    url_list = []
    for mon in range(1, 13):
        ndays = atm.days_this_month(year, mon)
        for day in range(1, ndays + 1):
            url = get_url(year, mon, day, datakind=datakind)
            url_list.append(url)
    writefile(url_list, savefile % year)
