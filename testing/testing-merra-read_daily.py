import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import numpy as np
import matplotlib.pyplot as plt
import merra
import atmos as atm
from merra import read_daily

year = 1979
month = 7

lon1, lon2 = 60, 100
lat1, lat2 = -60, 60

days=None
subset_dict = {'lat' : (lat1, lat2), 'lon' : (lon1, lon2)}
u = read_daily('u', year, month, days, subset_dict=subset_dict)

days=5
u2 = read_daily('u', year, month, days, subset_dict=subset_dict)

days = np.arange(10,13)
u3 = read_daily('u', year, month, days, subset_dict=subset_dict)

# Multiple variables
days = [10, 11]
ds = read_daily(['u', 'v', 'T'], year, month, days, subset_dict=subset_dict)
