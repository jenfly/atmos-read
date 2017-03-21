import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import xarray as xray
import pydap
import numpy as np
import collections
import atmos as atm

url = 'http://iridl.ldeo.columbia.edu/SOURCES/.OSU/.PRISM/.monthly/dods'

ds = xray.open_dataset(url, decode_times=False, engine='pydap')

ds2 = pydap.client.open_url(url)