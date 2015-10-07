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
u = read_daily('u', year, month, days, subset1=('lat', lat1, lat2),
               subset2=('lon', lon1, lon2))

days=5
u2 = read_daily('u', year, month, days, subset1=('lat', lat1, lat2),
                subset2=('lon', lon1, lon2))

days = np.arange(10,13)
u3 = read_daily('u', year, month, days, subset1=('lat', lat1, lat2),
                subset2=('lon', lon1, lon2))

# Multiple variables
days = [10, 11]
ds = read_daily(['u', 'v', 'T'], year, month, days,
                subset1=('lat', lat1, lat2), subset2=('lon', lon1, lon2))
