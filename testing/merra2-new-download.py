import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import xarray as xray
import pydap
import numpy as np
import collections
import atmos as atm

url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/'
       'M2I3NVASM.5.12.4/2015/08/MERRA2_400.inst3_3d_asm_Nv.20150802.nc4')

url2 = ('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/'
        'M2I1NXASM.5.12.4/2015/08/MERRA2_400.inst1_2d_asm_Nx.20150802.nc4')

from pydap.client import open_url
from pydap.cas.urs import setup_session
netrc_file = atm.homedir() + '.netrc'

def get_login(netrc_file):
    with open(netrc_file) as f:
        lines = f.readlines()
    username = lines[1].replace('\n','').split()[1]
    password = lines[2].replace('\n','').split()[1]
    auth = {'username' : username, 'password' : password}
    return auth

auth = get_login(netrc_file)

session = setup_session(auth['username'], auth['password'], check_url=url)
session2 = setup_session(auth['username'], auth['password'], check_url=url2)
ds = open_url(url, session=session)
ds2 = open_url(url2, session=session2)


print(ds['T'].shape)
print(ds['T'][0, 0, 0, 0])
print(ds2['TS'].shape)
print(ds2['TS'][0, 0, 0])

