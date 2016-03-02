import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import xray
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import atmos as atm
import precipdat
import merra

# ----------------------------------------------------------------------
#version, years = 'merra', np.arange(1979, 2015)
version, years = 'merra2', np.arange(1980, 2016)

datadir = atm.homedir() + 'datastore/%s/daily/' % version
months = np.arange(1, 13)
subset = '_40E-120E_90S-90N'

def get_var(datadir, version, varnm, subset, year):
    filenm = '%s%s_%s%s_%d.nc' % (datadir, version, varnm, subset, year)
    with xray.open_dataset(filenm) as ds:
        var = ds[varnm].load()
    return var

for year in years:
    print('Calculating MFC %d' % year)
    uq_int = get_var(datadir, version, 'UFLXQV', subset, year)
    vq_int = get_var(datadir, version, 'VFLXQV', subset, year)
    mfc = atm.moisture_flux_conv(uq_int, vq_int, already_int=True)
    mfc.attrs['long_name'] = mfc.name
    mfc.name = 'MFC'
    savefile = datadir + '%s_MFC%s_%d.nc' % (version, subset, year)
    print('Saving MFC to ' + savefile)
    atm.save_nc(savefile, mfc)
