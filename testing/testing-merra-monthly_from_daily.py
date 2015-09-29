import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-read')
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
import merra
import atmos as atm

print('Testing')
year = 1979
month = 1

# def filename(varname):
#     return '/home/jwalker/eady/scratch/' + varname + '.nc'
#
# # One variable, surface
# precip = merra.monthly_from_daily(year, month, 'precip')
# atm.save_nc(filename('precip'), precip)
#
# # One variable, pressure-level
# u = merra.monthly_from_daily(year, month, 'u')
# atm.save_nc(filename('u'), u)
#
# # Two variables, surface
# precip, evap, precip_evap = merra.monthly_from_daily(year, month, 'precip',
#                                                      'evap')
# atm.save_nc(filename('precip_evap'), precip, evap, precip_evap)
#
# # Two variables, pressure-level
# u, v, uv = merra.monthly_from_daily(year, month, 'u', 'v')
# atm.save_nc(filename('uv'), u, v, uv)
