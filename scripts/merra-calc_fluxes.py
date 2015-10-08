import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import atmos as atm
import merra
from merra import calc_fluxes

years = [1979]
months = [1,2,3]

scratchdir = '/net/eady/data1/jwalker/datastore/scratch/'
savedir = '/net/eady/data1/jwalker/datastore/merra/monthly/'
var_ids=['u', 'q', 'T', 'theta', 'theta_e', 'hgt']

def filename(varname, year, month):    
    datestr = '_%d%02d.nc' % (year, month)
    filen = savedir + varname + datestr
    print('Saving to ' + filen)
    return filen

def extract(ds, var_id):
    """Return variable plus its fluxes from dataset"""
    nm = merra.get_varname(var_id)
    uvar_nm = 'U*' + nm
    vvar_nm = 'V*' + nm
    ds_out = ds[nm].to_dataset()
    ds_out[uvar_nm] = ds[uvar_nm]
    ds_out[vvar_nm] = ds[vvar_nm]
    return ds_out

for year in years:
    for month in months:
        ds = calc_fluxes(year, month, var_ids=var_ids, scratchdir=scratchdir)
        for var_id in var_ids:
            ds_var = extract(ds, var_id)
            ds_var.to_netcdf(filename(var_id, year, month))
