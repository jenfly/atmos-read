"""Functions to read various precipitation datasets.

- CMAP
- GPCP
- TRMM
"""

from __future__ import division

import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')

import xray
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import atmos as atm

# ----------------------------------------------------------------------
def read_cmap(datafile, yearmin=None, yearmax=None, pentad_day=3):
    """Read CMAP pentad data for selected years.

    Parameters
    ----------
    datafile : str
        File path for CMAP data in NetCDF format.
    yearmin, yearmax : ints, optional
        Min and max years (inclusive) to extract.
    pentad_day : int, optional
        Which day of pentad to use for day coordinate.

    Returns
    -------
    precip : xray.DataArray
        CMAP pentad data for selected years, reshaped to
        [year, pentad, lat, lon].
    """

    NPYR = 73       # Number of pentads per year
    YRMIN = 1979    # First year of CMAP data

    if yearmin is None:
        yearmin = YRMIN

    # Read data
    with xray.open_dataset(datafile) as ds:
        precip = ds['precip'].load()

    # Discard incomplete year at end
    ny = precip.shape[0] // NPYR
    precip = precip[:ny*NPYR]

    # Split into individual years
    years = np.arange(YRMIN, YRMIN + ny)
    pentads = np.arange(1, NPYR + 1)
    precip = atm.split_timedim(precip, NPYR, time0_name='year',
                               time0_vals=years, time1_name='pentad',
                               time1_vals=pentads)

    # Extract selected years
    if yearmax is None:
        years = np.arange(yearmin, YRMIN + ny)
    else:
        years = np.arange(yearmin, yearmax + 1)
    precip = precip.sel(year=years)
    precip.coords['day'] = atm.pentad_to_jday(precip['pentad'], pmin=1,
                                              day=pentad_day)

    precip = precip.swap_dims({'pentad' : 'day'})
    return precip
