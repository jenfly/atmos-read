import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')
sys.path.append('/home/jwalker/dynamics/python/monsoon-onset')

import xarray as xray
import numpy as np
import collections
import pandas as pd
import scipy

import atmos as atm

datalist = 'troubleshooting/data-probs.log'
probdata = pd.read_csv(datalist, skipinitialspace=True)
savefile = 'data/merra_urls_ubudget.txt'

def get_prodnum(year):
    prod_dict = {yr: '100' for yr in range(1980, 1992)}
    for yr in range(1992, 2001):
        prod_dict[yr] = '200'
    for yr in range(2001, 2011):
        prod_dict[yr] = '300'
    for yr in range(2011, 2017):
        prod_dict[yr] = '400'
    prod = prod_dict[year]
    return prod


def get_url(year, mon, day, ana=False, lat1=-90, lat2=90, lon1=40, lon2=120):

    yrstr = '%d' % year
    mstr = '%02d' % mon
    dstr = '%02d' % day
    datestr = yrstr + mstr + dstr
    prod = get_prodnum(year)

    if ana:
        # DUDTANA
        shortname = 'M2T3NPUDT'
        longname = 'MERRA2_' + prod + '.tavg3_3d_udt_Np.' + datestr
        variables = '&VARIABLES=dudtana'
    else:
        # Regular pressure-level data
        shortname = 'M2I3NPASM'
        longname = 'MERRA2_' + prod + '.inst3_3d_asm_Np.' + datestr
        variables = '&VARIABLES=omega%2Ch%2Cu%2Cv'

    base_url = ('http://goldsmr5.gesdisc.eosdis.nasa.gov/daac-bin/OTF/' +
                'HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2F' + shortname +
                '.5.12.4%2F')
    sep = '%2C'
    box = ('&FORMAT=bmM0Lw&BBOX=' + str(lat1) + sep + str(lon1) + sep + str(lat2)
           + sep + str(lon2))
    label = '&LABEL=' + longname + '.SUB.nc4'
    layers = ('&LAYERS=' + 'LAYER_1%2C4%2C7%2C10%2C13%2C15%2C17%2C19'
              '%2C21%2C22%2C23%2C24%2C25%2C26%2C27%2C29%2C30')

    url = (base_url + yrstr + '%2F' + mstr + '%2F' + longname + '.nc4' + box + label
           + '&FLAGS=1&SHORTNAME=' + shortname + '&SERVICE=SUBSET_MERRA2' + layers
           +'&VERSION=1.02' + variables)

    return url

url_list = []

for i in range(len(probdata)):
    varnm = probdata['varnm'][i]
    if varnm == 'DUDTANA':
        ana = True
    else:
        ana = False
    year = probdata['year'][i]
    jday = probdata['jday'][i]
    mon, day = atm.jday_to_mmdd(jday, year)
    #print(probdata['yyyymmdd'][i], year, mon, day, varnm, ana)
    url = get_url(year, mon, day, ana)
    if url not in url_list:
        url_list.append(url)

print('Writing url list to ' + savefile)
with open(savefile, 'w') as f:
    for url in url_list:
        f.write(url + '\n')

# url_ana = ('http://goldsmr5.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?' +
#            'FILENAME=%2Fdata%2FMERRA2%2FM2T3NPUDT.5.12.4%2F1997%2F07%2F' +
#            'MERRA2_200.tavg3_3d_udt_Np.19970710.nc4&FORMAT=bmM0Lw&BBOX=' +
#            '-90%2C-180%2C90%2C180&LABEL=MERRA2_200.tavg3_3d_udt_Np.19970710.SUB.nc4' +
#            '&FLAGS=1&SHORTNAME=M2T3NPUDT&SERVICE=SUBSET_MERRA2&LAYERS=' +
#            'LAYER_1%2C2%2C3%2C4&VERSION=1.02&VARIABLES=dudtana')


#
# url_pres = ('http://goldsmr5.gesdisc.eosdis.nasa.gov/daac-bin/OTF/'
#             'HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2I3NPASM.5.12.4%2F'
#             '1997%2F06%2FMERRA2_200.inst3_3d_asm_Np.19970610.nc4&'
#             'FORMAT=bmM0Lw&BBOX=-90%2C40%2C90%2C120&LABEL='
#             'MERRA2_200.inst3_3d_asm_Np.19970610.SUB.nc4&FLAGS=1&'
#             'SHORTNAME=M2I3NPASM&SERVICE=SUBSET_MERRA2&LAYERS='
#             'LAYER_1%2C4%2C7%2C10%2C13%2C15%2C17%2C19%2C21%2C22%2C23%2C24%2C25'
#             '%2C26%2C27%2C29%2C30&VERSION=1.02&VARIABLES=omega%2Ch%2Cu%2Cv')

#
# years = range(1980, 2016)
# savedir = atm.homedir() + 'dynamics/python/atmos-read/scripts/merra_urls/'
# savefiles = {yr : savedir + 'merra2_rad_%d.txt' % yr for yr in years}
#
# def get_url(year, mon, day, lat1='-65', lat2='65', lon1='40', lon2='120'):
#     yrstr = '%d' % year
#     mstr = '%02d' % mon
#     dstr = '%02d' % day
#
#     prod_dict = {yr : '100' for yr in range(1980, 1992)}
#     for yr in range(1992, 2001):
#         prod_dict[yr] = '200'
#     for yr in range(2001, 2011):
#         prod_dict[yr] = '300'
#     for yr in range(2011, 2016):
#         prod_dict[yr] = '400'
#     prod = prod_dict[year]
#
#     url = ('http://goldsmr4.gesdisc.eosdis.nasa.gov/daac-bin/OTF/' +
#            'HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2T1NXRAD.5.12.4%2F' +
#             yrstr + '%2F' + mstr + '%2FMERRA2_' + prod + '.tavg1_2d_rad_Nx.' +
#             yrstr + mstr + dstr + '.nc4&FORMAT=bmM0Lw&BBOX=' + lat1 +
#             '%2C' + lon1 + '%2C' + lat2 + '%2C' + lon2 + '&LABEL=MERRA2_' +
#             prod + '.tavg1_2d_rad_Nx.' + yrstr + mstr + dstr +
#             '.SUB.nc4&FLAGS=1&SHORTNAME=M2T1NXRAD&SERVICE=SUBSET_MERRA2' +
#             '&LAYERS=&VERSION=1.02&VARIABLES=lwgnt%2Clwtup%2Cswgnt%2Cswtnt')
#     return url
#
# for year in years:
#     print(year)
#     with open(savefiles[year], 'w') as f:
#         for mon in range(1, 13):
#             print('  %d' % mon)
#             days = range(1, atm.days_this_month(year, mon) + 1)
#             for day in days:
#                 print('    %d' % day)
#                 url = get_url(year, mon, day)
#                 f.write(url + '\n')
