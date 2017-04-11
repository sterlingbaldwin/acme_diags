import logging, cdms2, MV2, numpy, cdutil, pdb
from numbers import Number
from math import log10
from unidata import udunits
from acme_diags.metrics import rmse, corr, acme_averager
from defines import all_regions
from acme_diags.derivations import acme
from table_parameter import *
from table_row_spec import table_row_specs
logger_fmt = "%(levelname)s: %(message)s"
logging.basicConfig(level=logging.ERROR, filename=None, format=logger_fmt)
logger = logging.getLogger(__name__)

def findfile(directory, prefix, season):
    import os
    DIRECTORY = os.listdir(directory)
    for file in DIRECTORY:
        if file.startswith(prefix) and season in file:
            return file
    return None

def reconcile_units( mv1, mv2, preferred_units=None ):
    """Changes the units of variables (instances of TransientVariable) mv1,mv2 to be the same,
    mv1 is a TransientVariable or an axis.  mv2 may be a TransientVariable or axis, or it may
    be a udunits object.
    If preferred units are specified, they will be used if possible."""

    # First, if there are no units, take a guess.  I'm reluctant to do this because it will surely
    # be wrong sometimes.  But usually it is correct.
    if not hasattr(mv2,'units') and hasattr(mv2,'lunits'):
        # Fix an apparent typo in two variables in CERES-EBAF data
        mv2.units = mv2.lunits
    if not hasattr(mv1,'units') or mv1.units == 'none':
        logger.warning("Variable %s has no units, will use units=1.", getattr(mv1,'id',''))
        mv1.units = '1'
    if not hasattr(mv2,'units') or mv2.units == 'none':
        logger.warning("Variable %s has no units, will use units=1.", getattr(mv2,'id',''))
        mv2.units = '1'

    # For QFLX and LHFLX variables, call dedicated functions instead.
    # TODO: The conditions for calling this function could perhaps be
    # better-defined.  I primarily wanted to ensure that it would
    # still work for the cases that the previous code worked for.
    #first check if they have an id

    if hasattr(mv1, 'id') and hasattr(mv2, 'id'):
        if mv1.id.find('_QFLX_')>=0 or mv1.id.find('_LHFLX_')>=0 or mv2.id.find('_QFLX_')>=0 \
                or mv2.id.find('_LHFLX_')>=0 or (mv1.units=='kg/m2/s' and mv2.units=='mm/day'):
            #from metrics.packages.amwg.derivations import qflx_lhflx_conversions as flxconv
            #mv1,mv2 = flxconv.reconcile_energyflux_precip(mv1, mv2, preferred_units)
            mv1, mv2 = reconcile_energyflux_precip(mv1, mv2, preferred_units)
            return mv1, mv2

    # I'm still checking for the units attribute here because maybe some time I'll be able to
    # delete the above lines which set it to a default value.
    if hasattr(mv1,'units') and hasattr(mv2,'units') and\
            (preferred_units is not None or mv1.units!=mv2.units):
        fix_troublesome_units( mv1 )
        fix_troublesome_units( mv2 )
        adhoc_convert_units( mv1, mv2.units )
        adhoc_convert_units( mv2, mv1.units )

        if preferred_units is None:
            # Look for whichever of mv1.units, mv2.units gives numbers more O(1).
            try:
                mv1min = mv1.min()
                mv1max = mv1.max()
            except AttributeError:
                # axes don't have min(),max() but they have getData() which variables don't have
                mv1min = mv1.getData().min()
                mv1max = mv1.getData().max()
            try:
                mv2min = mv2.min()
                mv2max = mv2.max()
            except AttributeError:
                # axes don't have min(),max() but they have getData() which variables don't have
                try:
                    mv2min = mv2.getData().min()
                    mv2max = mv2.getData().max()
                except AttributeError:
                    # udunits doesn't have min,getData, but has a scalar value
                    mv2min = mv2.value
                    mv2max = mv2.value
            mag1 = log10(abs(0.5*(mv1min+mv1max)))
            mag2 = log10(abs(0.5*(mv2min+mv2max)))
            if abs(mag1)<=abs(mag2):
                target_units = mv1.units
                changemv1 = False
                changemv2 = True
            else:
                target_units = mv2.units
                changemv1 = True
                changemv2 = False
        else:
            target_units = preferred_units
            changemv1 = True
            changemv2 = True
        if changemv1:
            tmp = udunits(1.0,mv1.units)
            try:
                s,i = tmp.how(target_units)
            except Exception as e:
                # conversion not possible.
                logger.error("Could not convert from %s to %s", mv1.units, target_units)
                if target_units == preferred_units:
                    pair = ("preferred units=", preferred_units)
                else:
                    pair = ("variable mv2=", getattr(mv2, 'id', '(not known)'))

                logger.error("units are from variable mv1=%s and %s %s", getattr(mv1, 'id', '(not known)'), pair[0], pair[1])
                raise e
            if hasattr(mv1,'id'):  # yes for TransientVariable, no for udunits
                mv1id = mv1.id
            if not ( numpy.allclose(s,1.0) and numpy.allclose(i,0.0) ):
                # The following line won't work if mv1 be an axis.
                mv1 = s*mv1 + i
                if hasattr(mv1,'id'):
                    mv1.id = mv1id
                if hasattr(mv1,'mean') and isinstance(mv1.mean,Number):
                    mv1.mean = s*mv1.mean + i
            mv1.units = target_units
        if changemv2:
            tmp = udunits(1.0,mv2.units)
            try:
                s,i = tmp.how(target_units)
            except Exception as e:
                #  conversion not possible
                logger.error("Could not convert from %s to %s",mv2.units, target_units)
                if target_units==preferred_units:
                    pair = "preferred units=",preferred_units
                else:
                    pair = "variable mv1=",getattr(mv1,'id','(not known)')
                logger.error("units are from variable mv2=%s and %s%s", getattr(mv2,'id','(not known)'),pair[0], pair[1])
                raise e
            if hasattr(mv2,'id'):  # yes for TransientVariable, no for udunits
                mv2id = mv2.id
            if not ( numpy.allclose(s,1.0) and numpy.allclose(i,0.0) ):
                # The following line won't work if mv1 be an axis.
                mv2 = s*mv2 + i
                if hasattr(mv2,'id'):
                    mv2.id = mv2id
                if hasattr(mv2,'mean') and isinstance(mv2.mean,Number):
                    mv2.mean = s*mv2.mean + i
            mv2.units = target_units
    return mv1, mv2

def select_lev( mv, slev ):
    """Input is a level-dependent variable mv and a level slev to select.
    slev is an instance of udunits - thus it has a value and a units attribute.
    This function will create and return a new variable mvs without a level axis.
    The values of mvs correspond to the values of mv with level set to slev.
    Interpolation isn't done yet, but is planned !"""
    levax = mv.getLevel()
    if levax is None:
        return None
    # Get ig, the first index for which levax[ig]>slev
    # Assume that levax values are monotonic.
    dummy,slev = reconcile_units( levax, slev )  # new slev has same units as levax
    if levax[0]<=levax[-1]:
        ids = numpy.where( levax[:]>=slev.value )    # assumes levax values are monotonic increasing
    else:
        ids = numpy.where( levax[:]<=slev.value )    # assumes levax values are monotonic decreasing
    if ids is None or len(ids)==0:
        ig = len(levax)-1
    else:
        ig = ids[0][0]
    # Crude first cut: don't interpolate, just return a value
    if levax == mv.getAxisList()[0]:
        mvs = cdms2.createVariable( mv[ig:ig+1,...], copy=1 )  # why ig:ig+1 rather than ig?  bug workaround.
    elif levax == mv.getAxisList()[1]:
        mvs = cdms2.createVariable( mv[:,ig:ig+1,...], copy=1 )
    else:
        logger.error("select_lev() does not support level axis except as first or second dimensions")
        return None
    # createVariable will have copied all non-'_' attributes from mv to mvs; the mv mean is wrong for mvs.
    mvs.mean = None
    if hasattr(mv, 'units'): mvs.units = mv.units
    mvs = mvs(squeeze=1)
    return mvs

def convert_units(mv, units):
   """Returns a new variable like mv but with the specified units"""
   if type(mv) == float:
      return mv
   if not hasattr(mv,'units') and hasattr(mv,'lunits'):
       # Fix an apparent typo in two variables in CERES-EBAF data
       mv.units = mv.lunits
   fix_troublesome_units(mv)
   adhoc_convert_units( mv, units )
   if not hasattr(mv, 'units'):
      mv.units = '1'
   if mv.units == units or mv.units == 'none':
      return mv
   if mv.units=='gC/m^2/s' and units == 'PgC/y': # land set 5 table stuff.
      mv = mv * (60*60*24*365)
      mv = mv / (1.0e15)
      mv.units = units
      return mv

   tmp = udunits(1.0,mv.units)
   try:
      s,i = tmp.how(units)
   except Exception as e:
      # conversion not possible.
      logger.error("Could not convert from %s to %s",mv.units,units)
      return mv
   if hasattr(mv,'id'):  # yes for TransientVariable, no for udunits
      mvid = mv.id
   if not ( numpy.allclose(s,1.0) and numpy.allclose(i,0.0) ):
      # The following line won't work if mv1 is an axis.
      mvmean = mv.mean
      mv = s*mv + i
      try:
          mv.mean = s*mvmean + 1
          if hasattr(mv,'_mean'):
              mv._mean = mv.mean
      except:
          pass  # probably mv.mean is None, or a function
      if hasattr(mv,'id'):
         mv.id = mvid
   mv.units = units
   return mv

def fix_troublesome_units( mv ):
    """This handles many special cases where a variable mv has units which are not understood
    by udunits or otherwise troublesome.  It replaces the units with ones which work better."""
        # Very ad-hoc, but more general would be less safe:
        # BES - set 3 does not seem to call it rv_QFLX. It is set3_QFLX_ft0_climos, so make this just a substring search
    if not hasattr(mv,'units'):
        return mv
    if mv.units == "gpm":
        mv.units="m"
    if mv.units=='mb':
        mv.units = 'mbar' # udunits uses mb for something else
    if mv.units=='mb/day':
        mv.units = 'mbar/day' # udunits uses mb for something else
    if mv.units == '(0 - 1)': # as in ERAI obs
        mv.units = '1'
    if mv.units == '(0-1)':   # as in ERAI obs
        mv.units = '1'
    if mv.units == 'fraction' or mv.units=='dimensionless':
        mv.units = '1'
    if mv.units == 'mixed':  # could mean anything...
        mv.units = '1'       #... maybe this will work
    if mv.units == 'unitless':  # could mean anything...
        mv.units = '1'       #... maybe this will work
    if mv.units == 'W/m~S~2~N~' or mv.units == 'W/m~S~2~N':
        mv.units = 'W/m^2'
    if hasattr(mv,'filetable') and mv.filetable.id().ftid == 'ERA40' and\
            mv.id[0:5]=='rv_V_' and mv.units=='meridional wind':
        # work around a silly error in ERA40 obs
        mv.units = 'm/s'
    return mv

def adhoc_convert_units( mv, units ):
    """Some ad-hoc units conversions"""
    if not hasattr(mv,'units'):
        return mv
    if mv.units =='fraction' and (units == 'percent' or units == '%'):
        mv = 100*mv
        mv.units=units
    if mv.units=='kg/m2' and units=='mm':
        mv.units = 'mm' # [if 1 kg = 10^6 mm^3 as for water]
    if mv.units[0:7]=='kg/(m^2' and units[0:3]=='mm/':
        mv.units = 'mm/('+mv.units[7:] # [if 1 kg = 10^6 mm^3 as for water]
    if mv.units[0:3]=='mm/' and units[0:8]=='kg/(m^2 ':
        mv.units = 'kg/(m^2 '+mv.units[3:]+')' # [if 1 kg = 10^6 mm^3 as for water]
    return mv

def reconcile_energyflux_precip(mv1, mv2, preferred_units=None):
    # To compare LHFLX and QFLX, need to unify these to a common variables
    # e.g. LHFLX (latent heat flux in W/m^2) vs. QFLX (evaporation in mm/day)
    #
    # If preferred_units is not provided, the units of mv2 will be
    # assumed to be the preferred units.
    #
    # This function is used by the derived_variable definitions in
    # amwg_plot_plan's standard_variables (within amwg.py).
    #
    # Author: S.M. Burrows, 9 Feb 2015.

    # If preferred_units is not provided, assume units of mv2 assumed
    # to be the preferred units.

    if hasattr(mv1,'units') and hasattr(mv2,'units'):

        # First, set preferred_units if needed
        if preferred_units is None:
            if ('_QFLX_' in mv2.id) or ('_QFLX' in mv1.id):
                logger.info("Setting preferred_units='mm/day'")
                preferred_units='mm/day'
            if ('_LHFLX_' in mv2.id) or ('_LHFLX' in mv1.id):
                logger.info("Setting preferred_units='W/m^2'")
                preferred_units='W/m^2'
            if preferred_units is None:
                logger.info("Setting preferred_units to mv.units=%s", mv2.units)
                preferred_units = mv2.units

        # syntax correction (just in case)
        if preferred_units=='W/m2':
            preferred_units='W/m^2'

        # Now do conversions to preferred_units (only if needed)
        if mv1.units!=preferred_units:
            mv1 = convert_energyflux_precip(mv1, preferred_units)
        if mv2.units!=preferred_units:
            mv2 = convert_energyflux_precip(mv2, preferred_units)
    else:
        logger.error("missing units in arguments to reconcile_energyflux_precip.")
        exit

    return mv1,mv2

def convert_energyflux_precip(mv, preferred_units):

    # The latent heat of vaporization for water is 2260 kJ/kg
    lhvap = 2260. # 'kJ/kg'
    secondsperday = 86400.
    kJperday = 86.4 # 'kJ/day'

    if hasattr(mv,'id'):
        mvid = mv.id

    # syntax correction (just in case)
    mv.units = mv.units.replace(' m-2','/m^2')
    mv.units = mv.units.replace(' s-1','/s')
    if  mv.units=='W/m2':
        mv.units='W/m^2'
    if mv.units=='mm/d':
        mv.units = 'mm/day'
    # LHFLX
    if mv.units=="W/m~S~2~N":
        logger.info("Arbitrarily decided that W/m~S~2~N is W/m^2 for %s", mv.id)
        mv.units="W/m^2"
    if mv.units=="W/m~S~2~N~":
        logger.info("Arbitrarily decided that W/m~S~2~N~ is W/m^2 for %s", mv.id)
        mv.units="W/m^2"

    if mv.units==preferred_units:
        return mv

    # convert precip between kg/m2/s and mm/day
    if ( mv.units=="kg/m2/s" or mv.units=="kg/m^2/s" or mv.units=="kg/s/m2" or\
             mv.units=="kg/s/m^2") and preferred_units=="mm/day":
        mv = mv * secondsperday # convert to kg/m2/s [= mm/s]
        mv.units="mm/day"         # [if 1 kg = 10^6 mm^3 as for water]

    elif mv.units=='mm/day' and preferred_units=="kg/m2/s":
        mv = mv / secondsperday # convert to mm/sec [= kg/m2/s]
        mv.units="kg/m2/s"      # [if 1 kg = 10^6 mm^3 as for water]

    # convert between energy flux (W/m2) and water flux (mm/day)
    elif mv.units=="kg/m2/s" and preferred_units=="W/m^2":
        mv = mv * kJperday * secondsperday * lhvap
        mv.units = 'W/m^2'

    elif mv.units=='mm/day' and preferred_units=='W/m^2':
        # 1 W = 86.4 kJ / day
        mv = mv * lhvap / kJperday
        mv.units = 'W/m^2'

    elif mv.units=='W/m^2' and preferred_units=='mm/day':
        mv = mv * kJperday / lhvap
        mv.units = 'mm/day'

    else:
        tmp = udunits(1.0,mv.units)
        try:
            s,i = tmp.how(preferred_units)
        except Exception as e:
            # conversion not possible.
            logger.error("could not convert from %s to %s" ,mv.units, preferred_units)
            raise e
        if not ( numpy.allclose(s,1.0) and numpy.allclose(i,0.0) ):
            mv = s*mv + i
        mv.units = preferred_units
    mv.id = mvid # reset variable id

    return mv

def regrid_to_lower_res(mv1, mv2, regrid_tool, regrid_method):
    """regrid transient variable toward lower resolution of two variables"""

    axes1 = mv1.getAxisList()
    axes2 = mv2.getAxisList()

    # use nlat to decide data resolution, higher number means higher data
    # resolution. For the difference plot, regrid toward lower resolution
    if len(axes1[1]) <= len(axes2[1]):
        mv_grid = mv1.getGrid()
        mv1_reg = mv1
        mv2_reg = mv2.regrid(mv_grid, regridTool=regrid_tool,
                             regridMethod=regrid_method)
    else:
        mv_grid = mv2.getGrid()
        mv2_reg = mv2
        mv1_reg = mv1.regrid(mv_grid, regridTool=regrid_tool,
                             regridMethod=regrid_method)
    return mv1_reg, mv2_reg

def verticalize( T, hyam, hybm, ps ):
    """
    For data T with CAM's hybrid level coordinates, interpolates to
    the more standard pressure level coordinates and returns the results.
    The input arguments hyam, hybm, ps are the usual CAM veriables by that
    name.  Order of dimensions must be (lev,lat,lon).
    """
    # constants as in functions_vertical.ncl, lines 5-10:
    plvlO = numpy.array([30., 50., 70., 100., 150., 200., 250., 300., 400., 500.,
                         600., 700., 775., 850., 925., 1000.])  # in mb

    p0 = 1000.   # mb
    if T.getLevel() is None:
        return None

    # Convert p0 to match ps.  Later, we'll convert back to mb.  This is faster than
    # converting ps to millibars.
    if ps.units=='mb':
        ps.units = 'mbar' # udunits uses mb for something else
    tmp = udunits(1.0,'mbar')
    s,i = tmp.how(ps.units)
    p0 = s*p0 + i
    levels_orig = cdutil.vertical.reconstructPressureFromHybrid( ps, hyam, hybm, p0 )

    # At this point levels_orig has the same units as ps.  Convert to mbar
    tmp = udunits(1.0,ps.units)
    s,i = tmp.how('mbar')
    levels_orig = s*levels_orig + i
    levels_orig.units = 'mbar'
    newT = cdutil.vertical.logLinearInterpolation( T, levels_orig, plvlO )
    return newT

def convert_variable( var, target_units ):
    """Converts a variable (cdms2 MV) to the target units (a string) if possible, and returns
    the variable modified to use the new units."""

    unit_synonyms = { 'mb': 'millibar', 'pa': 'pascal', 'gpm': 'meters', 'C': 'Celsius' }
    last_ditch_conversions = { 'fraction': ('percent', 100.0, 0.0) }

    if not hasattr( var, 'units' ):
        return var
    if target_units==None or target_units=='':
        return var
    if var.units == target_units:
        return var
    if var.units in unit_synonyms.keys():
        var.units = unit_synonyms[var.units]
    if var.units in last_ditch_conversions.keys():
        u,s,o = last_ditch_conversions[var.units]
        var = s*var + o
        var.units = u
    try:
        s,o = udunits(1.0,var.units).how(target_units)
    except TypeError as e:
        logging.warning("Could not convert units from %s to %s",var.units,target_units)
        return var
    var = s*var + o
    var.units = target_units
    return var

def create_metrics(ref, test, gw):
    """ Creates the mean, max, min, rmse, corr in a dictionary """

    metrics_dict = {
        'model_mean': acme_averager(test, '(lat)(lon)', weights=gw).item(),
        'obs_mean': acme_averager(ref, '(lat)(lon)', weights=gw).item(),
        'rmse': rmse(test, ref),
        'corr': corr(test, ref)
    }

    return metrics_dict

def ncar_mask(x, y):
    """Apply the mask of x to y and y to x. This seems to apply only for SST and HadISST"""
    X = MV2.masked_where(y.mask, x)
    Y = MV2.masked_where(x.mask, y)
    return X, Y

def get_data(var_file, varid, season):
    global first_model_data_read, hybrid, hyam, hybm, PS
    vars = []
    f = cdms2.open(var_file)

    try:

        #get the data for the requested variable; it may be a derived variable
        try:
            if varid in f.variables.keys():
                var = f(varid)(squeeze=1)
                vars.append(var)
            else:
                var = acme.process_derived_var( varid, acme.derived_variables, f, None)
                if var is not None:
                    vars.append(var)
                else:
                    raise
        except:
            raise

        #get gaussian weights if present
        try:
            gw = f('gw')(squeeze=1)
            vars.append(gw)
        except:
            vars.append(None)

        #get the hybrid variables only once if present
        if first_model_data_read:
            first_model_data_read = False
            try:
                hyam = f('hyam')(squeeze=1)
                hybm = f('hybm')(squeeze=1)
                PS = f('PS')(squeeze=1)

                hybrid = True
            except:
                hybrid = False

    except:
        f.close()
        errmsg = "no data for " + varid + " in " + var_file
        logger.error(errmsg)
        return errmsg
    f.close()
    return vars

def print_table(rows):
    def fpfmt( num ):
        """No standard floating-point format (e,f,g) will do what we need, so this switches between f and e"""
        if num is undefined:
            return '        '
        if num>10001 or num<-1001:
            return format(num,"10.4e")
        else:
            return format(num,"10.3f")

    undefined = -numpy.infty
    name = 'Tables of Global, tropical, and extratropical DJF, JJA, ANN means and RMSE'
    title = ' '.join(['AMWG Diagnostics Set 1', season, 'means', region]) + '\n'
    subtitles = [
        ' '.join(['Test Case:', model_file]) + '\n',
        'Control Case: various observational data\n',
        'Variable                 Test Case           Obs          Test-Obs           RMSE           Correlation\n']

    print title
    print ' '.join(subtitles)
    for row in rows:
        rowname, values = row[0], row[1:]
        rowpadded = (rowname + 10 * ' ')[:20]
        output = [rowpadded]
        for v in values:
            if type(v) is not str:
                output.append( fpfmt(v) )
            else:
                output.append( v )
        print '\t'.join(output)

def compute_row(spec):
    global hybrid, hyam, hybm, PS

    #currently unused
    latmin, latmax, lonmin, lonmax = all_regions[region]

    varid = spec['var']

    obs_prefix = spec['obs']
    obs_fn = findfile(obs_path, obs_prefix, season)
    obs_file = obs_path + obs_fn

    level = spec.get('lev', None)
    if level:
        ulevel = udunits(float(level), 'mbar')

    units = spec.get('units', None)

    rowname = varid + '_'+obs_prefix
    if level:
        rowname += '_' + level

    #get the model and obs data
    model_data = get_data(model_file, varid, season)
    if type(model_data) is str:
        return [rowname, model_data]
    obs_data = get_data(obs_file, varid, season)
    if type(obs_data) is str:
        return [rowname, obs_data]

    #prepare model data
    model, weights = model_data
    gw = None
    if use_weights:
        gw = weights
    if level:
        if hybrid:
            #make sure the level axes are the same for the model and PS
            level_src = model.getLevel()
            PS, level_src = reconcile_units(PS, level_src, preferred_units='mbar')
            level_src = convert_variable(level_src.getValue(), PS.units)
            level_src = cdms2.createAxis(level_src)
            level_src.units = PS.units
            model.setAxis(model.getAxisIndex('lev'), level_src)
            model = verticalize(model, hyam, hybm, PS)
        model = select_lev(model, ulevel)
    if units:
        model = convert_units(model, units)

    #prepare obs data
    obs, dummy1 = obs_data #obs rarely has gw if ever
    if level:
        obs = select_lev(obs, ulevel)
    if units:
        obs = convert_units(obs, units)

    #put model and obs on the same grid
    model_regrid, obs_regrid = regrid_to_lower_res( model, obs, regridTool, regridMethod )

    #special masking for SST
    #Take note. If the model data file does not have SST but does have TS and OCNFRAC(this is common) then
    #SST is derived according to these parameters by masking.  So using the mcart mask is a second mask.
    if varid == 'SST' and obs_prefix.startswith('HadISST') and use_ncar_mask:
        model_new, obs_new = ncar_mask( model_regrid, obs_regrid)

    #the data preparation is complete. perform the calculations
    metrics_dict = create_metrics( obs_regrid, model_regrid, gw )

    return [rowname, metrics_dict['model_mean'], metrics_dict['obs_mean'],
            metrics_dict['model_mean'] - metrics_dict['obs_mean'],
            metrics_dict['rmse'], metrics_dict['corr']]

#one reference to hybrid coordinates so they are read once
first_model_data_read = True
hybrid = False
hyam = None
hybm = None
PS = None

rows = []
for spec in table_row_specs:
    row = compute_row(spec)
    rows.append( row )

print_table( rows )