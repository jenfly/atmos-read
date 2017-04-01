import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import xarray as xray
import pydap
import numpy as np
import collections
import atmos as atm

from pydap.client import open_url
from pydap.cas.urs import setup_session
netrc_file = atm.homedir() + '.netrc'

url_3d = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/'
         'M2I3NVASM.5.12.4/2015/08/MERRA2_400.inst3_3d_asm_Nv.20150802.nc4')

url_2d = ('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/'
          'M2I1NXASM.5.12.4/2015/08/MERRA2_400.inst1_2d_asm_Nx.20150802.nc4')

def get_login(netrc_file):
    with open(netrc_file) as f:
        lines = f.readlines()
    username = lines[1].replace('\n','').split()[1]
    password = lines[2].replace('\n','').split()[1]
    auth = {'username' : username, 'password' : password}
    return auth

auth = get_login(netrc_file)
username, password = auth['username'], auth['password']

def read_file(url, username, password):
    session = setup_session(username, password, check_url=url)
    ds = open_url(url, session=session)
    return ds

ds_3d = read_file(url_3d, username, password)
ds_2d = read_file(url_2d, username, password)


print(ds_3d['T'].shape)
print(ds_3d['T'][0, 0, 0, 0])
print(ds_2d['TS'].shape)
print(ds_2d['TS'][0, 0, 0])

#url = 'http://iridl.ldeo.columbia.edu/SOURCES/.OSU/.PRISM/.monthly/dods'
url = 'http://test.opendap.org/dap/data/nc/coads_climatology.nc'
ds = xray.open_dataset(url, engine='pydap', decode_times=False, decode_coords=False)

url3 = url_3d
session = setup_session(username, password, check_url=url3)
store = xray.backends.PydapDataStore(url3, session=session)
ds3 = xray.open_dataset(store, decode_times=False, decode_coords=False)