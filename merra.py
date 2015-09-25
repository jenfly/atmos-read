import numpy as np
import matplotlib.pyplot as plt
import xray

import atmos.utils as utils
import atmos.plots as ap
import atmos.data as dat
import atmos.variables as av
from atmos.constants import const as constants
from atmos.data import get_coord

# ----------------------------------------------------------------------
# Read data
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
# Hourly precip
#
# url3 = ('http://goldsmr2.sci.gsfc.nasa.gov/opendap/MERRA/MAT1NXFLX.5.2.0/'
#     '1979/07/MERRA100.prod.assim.tavg1_2d_flx_Nx.19790701.hdf')
# ds3 = xray.open_dataset(url3)
#
# precip = ds3['PRECTOT']
