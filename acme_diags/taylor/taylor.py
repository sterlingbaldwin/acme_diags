# AMWG Diagnostics, plot set 14.
# Here's the title used by NCAR:
# DIAG Set 14 - Taylor diagrams

import logging

import cdms2, MV2
import cdutil.times
import vcs

from acme_diags.metrics.compute_rmse import compute_rmse
from acme_diags.metrics.regrid_to_common_grid import regrid_to_common_grid
from acme_diags.metrics.reductions import reduce2scalar_seasonal_zonal
from parameter import *

logger = logging.getLogger(__name__)

seasonsyr=cdutil.times.Seasons('JFMAMJJASOND')

def get_data(var_file, varid, season):
    f = cdms2.open(var_file)

    try:
        if varid in f.variables.keys():
            var = f(varid)(squeeze=1)
    except:
        f.close()
        errmsg = "no data for " + varid + " in " + var_file
        logger.error(errmsg)
        return errmsg
    f.close()
    return var

def plot(taylor_data):
    colors = ['red', 'green', 'blue']
    cnvs = vcs.init()
    td = cnvs.createtaylordiagram('Taylor diagram')
    n = len(taylor_data)
    td.Marker.color = colors[0:n]
    td.Marker.size = n*[3]
    cnvs.plot(taylor_data, td)

    #plot the legend
    lx = .75
    ly = .95
    for i, ltitle in enumerate(test_names):
        text = cnvs.createtext()
        text.string = str(i) + '  ' + ltitle
        text.x = lx
        text.y = ly
        text.height = 14
        cnvs.plot(text, bg=1)
        ly -= .025

#from genutil.statistics import correlation
data = {}
data[ref_name] = get_data(ref_file, varid, season)
for test_name, test_file in zip(test_names, test_files):
    var = get_data(test_file, varid, season)
    if var is str:
        continue
    data[test_name] = var

#compute mean and std for each test
moments = {}
for name, var in data.items():
    mean = MV2.array(reduce2scalar_seasonal_zonal(var, gw=None))
    centered_var = var - mean
    std  = MV2.array(centered_var.std())
    moments[name] = (mean, std)

#compute the data to be plotted
taylor_data = []
ref_mean, ref_std = moments[ref_name]
ref_data = data[ref_name] - ref_mean
for test_name, test_file in zip(test_names, test_files):
    mean, std = moments[test_name]
    test_data = data[test_name] - mean
    ref_data, test_data = regrid_to_common_grid(ref_data, test_data, regridMethod=regridMethod, regridTool=regridTool)
    rmse, corr = compute_rmse(ref_data, test_data)
    taylor_data += [[std.item()/ref_std.item(), corr]]

plot(taylor_data)