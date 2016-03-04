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
years = np.arange(1979, 2016)
# version = 'merra2'
# years = np.arange(1980, 2016)

datadir = atm.homedir() + 'datastore/' + version + '/monthly/'
months = np.arange(1, 13)

varnms = ['U', 'V']

nc_kw = { 'merra2' : {'format' : 'NETCDF4_classic', 'engine' : 'netcdf4'},
          'merra' : {'format' : None, 'engine' : None}}[version]


for year in years:
    urls = merra.get_urls(year, version=version, varnm=varnms[0], monthly=True)
    for month in months:
        url = urls['%d%02d' % (year, month)]
        print('Loading ' + url)
        with xray.open_dataset(url) as ds:
            for nm in varnms:
                var = ds[nm]
                filenm = '%s%s_%s_%d%02d.nc' % (datadir, version, nm, year,
                                               month)
                print('Saving to ' + filenm)
                var.to_dataset().to_netcdf(filenm, **nc_kw)
