import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import numpy as np
import xray
import pandas as pd
import matplotlib.pyplot as plt

import atmos as atm
import merra

savedir = atm.homedir() + '/datastore/merra/daily/'
filestr = '%smerra_%s_ML%02d_40-120E_60S-60N_%d%02d.nc'

varnms = ['T', 'PS', 'DELP', 'QV', 'V']
years = range(1979, 2000)
months = [4, 5, 6, 7, 8, 9]

lev = 71
xsub = '[330:2:450]'
ysub = '[60:2:301]'

for year in years:
    for month in months:
        for varnm in varnms:
            savefile = filestr % (savedir, varnm, lev, year, month)
            var = merra.read_daily_eta(varnm, lev, year, month,
                                       concat_dim='TIME', xsub=xsub, ysub=ysub)
            print('Saving to ' + savefile)
            atm.save_nc(savefile, var)
