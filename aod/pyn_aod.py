import numpy as np

# Assumptions:
#
#  -- Arrays are passed as velocity, flux
#  -- Velocity spacing is constant (enough)


def fix_unwriteable_spec(spec):
    # FIX NON-WRITEABLE ARRAYS due to discontiguous memory
    for kkk in spec.keys():
        if isinstance(spec[kkk],(np.ndarray)):
            spec[kkk] = spec[kkk].copy()

    return spec

def xlimit(x, limits):

   def _ret():  return idx1, idx2

   idx1 = (np.abs(x - limits[0])).argmin()
   idx2 = (np.abs(x - limits[1])).argmin()

   return _ret()

def integrate_column(velocity, flux, flux_err,
                continuum, continuum_err, wavc, fval,
                integration_limits = [-100,100]):
    # TODO: Asymmetrical error bars?

    # Some constants and flags
    column_factor = 2.654e-15
    flag_sat = False

    # TODO: Adjust integration limits to reflect whole pixels used.
    # Define the limits of the integration:
    if not integration_limits:
        integration_limits = [spec['v1'],spec['v2']]
        idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,[spec['v1'],spec['v2']])
        idx[xlim1:xlim2] = True
    else:
        idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,integration_limits)
        idx[xlim1:xlim2] = True

    # An array of delta v:
    delv = velocity[1:]-velocity[:-1]
    delv = np.concatenate((delv,[delv[-1]]))

    # Test for clearly saturated pixels:
    idx_saturation = (flux <= 0.)
    if idx_saturation.sum() > 0: flag_sat = True

    # Fix saturation if it's present.
    flux[idx_saturation] = np.abs(flux[idx_saturation])
    flux[(flux==0)] = 2.*flux_err[(flux==0)]

    # Create an optical depth array and its error
    tau_array = np.log(continuum / flux)
    tau_array_err = np.sqrt((flux_err/flux)**2)

    tau_int = np.sum(tau_array[idx]*delv[idx])
    tau_int_err = \
     np.sqrt(np.sum((tau_array_err[idx]*delv[idx])**2))

    # Create an apparent column density array
    nav_array = tau_array/(wavc*fval*column_factor)
    nav_err = tau_array_err/(wavc*fval*column_factor)

    # Integrate the apparent column density profiles
    column = tau_int/(wavc*fval*column_factor)

    # Error in the column
    column_err = tau_int_err/(wavc*fval*column_factor)

    # Continuum error: errors are correlated, so don't add in quadrature.
    column_err_cont = \
        np.sum(((continuum_err[idx]/continuum[idx])*delv[idx])) /\
         (wavc*fval*column_factor)

    # Background uncertainty
    z_eps = 0.01  # Fractional bg error
    yc1 = continuum[idx]*(1.-z_eps)
    y1  = flux[idx]-continuum[idx]*z_eps
    tau1 = np.sum(np.log(yc1/y1)*delv[idx])
    col1 = tau1 / (wavc*fval*column_factor)
    column_err_zero = np.abs(col1-column)

    # Combine errors
    column_err_total = np.sqrt(column_err**2 \
        +column_err_cont**2)

    return column, column_err_total, flag_sat


def pyn_column(spec_in,integration_limits = None):

    spec = spec_in.copy()

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)

    # Some constants and flags
    column_factor = 2.654e-15
    flag_sat = False

    velocity = spec['vel'].copy()
    flux = spec['flux'].copy()
    flux_err = spec['eflux'].copy()
    wavc=spec['wavc'].copy()
    fval=spec['fval'].copy()

    # Deal with the continuum:
    if "contin" in spec.keys():
        continuum=spec['contin'].copy()
    else:
        try:
            continuum=spec['cont'].copy()
        except:
            continuum=spec['ycon'].copy()

    if "contin_err" in spec.keys():
        continuum_err = spec['contin_err'].copy()
    else:
        try:
            continuum_err = spec['econt'].copy()
        except:
            continuum_err = spec['ycon_sig'].copy()


    # Define the limits of the integration:
    if not integration_limits:
        integration_limits = [spec['v1'],spec['v2']]
        idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,[spec['v1'],spec['v2']])
        idx[xlim1:xlim2] = True
    else:
        idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,integration_limits)
        idx[xlim1:xlim2] = True

    # An array of delta v:
    delv = velocity[1:]-velocity[:-1]
    delv = np.concatenate((delv,[delv[-1]]))

    # Test for clearly saturated pixels:
    idx_saturation = (flux <= 0.)
    if idx_saturation.sum() > 0: flag_sat = True

    # Fix saturation if it's present.
    flux[idx_saturation] = np.abs(flux[idx_saturation])
    flux[(flux==0)] = 2.*flux_err[(flux==0)]

    # Create an optical depth array and its error
    tau_array = np.log(continuum / flux)
    tau_array_err = np.sqrt((flux_err/flux)**2)

    tau_int = np.sum(tau_array[idx]*delv[idx])
    tau_int_err = \
     np.sqrt(np.sum((tau_array_err[idx]*delv[idx])**2))

    # TODO: Include an AOD array in output w/continuum errors.
    # Create an apparent column density array
    nav_array = tau_array/(wavc*fval*column_factor)
    nav_err = tau_array_err/(wavc*fval*column_factor)

    # Integrate the apparent column density profiles
    column = tau_int/(wavc*fval*column_factor)

    # Error in the column
    column_err = tau_int_err/(wavc*fval*column_factor)

    # Continuum error: errors are correlated, so don't add in quadrature.
    column_err_cont = \
        np.sum(((continuum_err[idx]/continuum[idx])*delv[idx])) /\
         (wavc*fval*column_factor)

    # Background uncertainty
    z_eps = 0.01  # Fractional bg error
    yc1 = continuum[idx]*(1.-z_eps)
    y1  = flux[idx]-continuum[idx]*z_eps
    tau1 = np.sum(np.log(yc1/y1)*delv[idx])
    col1 = tau1 / (wavc*fval*column_factor)
    column_err_zero = np.abs(col1-column)

    # Combine errors
    column_err_total = np.sqrt(column_err**2 \
        +column_err_cont**2)

    spec['v1'] = integration_limits[0]
    spec['v2'] = integration_limits[1]
    spec['vaod1'] = integration_limits[0]
    spec['vaod2'] = integration_limits[1]

    spec['ncol'] = np.log10(column)
    spec['ncol_err_lo'] = -column_err_total/column*np.log10(np.e)
    spec['ncol_err_hi'] = column_err_total/column*np.log10(np.e)

    spec['flag_sat'] = flag_sat

    try:
        del spec['ncole1']
        del spec['ncole2']
        del spec['ncolez']
    except:
        pass

    return spec


def pyn_eqwidth(spec_in,integration_limits = None):

    spec = spec_in.copy()

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)

    # Some constants and flags
    column_factor = 2.654e-15
    ew_factor = 1.13e17
    lightspeed = 2.998e5 # km/s

    velocity = spec['vel'].copy()
    flux = spec['flux'].copy()
    flux_err = spec['eflux'].copy()
    wavc=spec['wavc'].copy()
    fval=spec['fval'].copy()

    # Deal with the continuum
    if "contin" in spec.keys():
        continuum=spec['contin'].copy()
    else:
        try:
            continuum=spec['cont'].copy()
        except:
            continuum=spec['ycon'].copy()

    if "contin_err" in spec.keys():
        continuum_err = spec['contin_err'].copy()
    else:
        try:
            continuum_err = spec['econt'].copy()
        except:
            continuum_err = spec['ycon_sig'].copy()


    # Define the limits of the integration:
    if not integration_limits:
        integration_limits = [spec['v1'],spec['v2']]
        idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,[spec['v1'],spec['v2']])
        idx[xlim1:xlim2] = True
    else:
        idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,integration_limits)
        idx[xlim1:xlim2] = True

    # An array of delta v:
    delv = velocity[1:]-velocity[:-1]
    delv = np.concatenate((delv,[delv[-1]]))

    # Create the wavelength array
    try:
        wave=spec['wave'].copy()
    except:
        wave = spec['wavc']*(velocity/lightspeed)+spec['wavc']
        spec['wave'] = wave

    # An array of delta wavelength:
    delw = wave[1:]-wave[:-1]
    delw = np.concatenate((delw,[delw[-1]]))

    # Calculate the equivalent width
    eqw_int = np.sum((1.-flux[idx]/continuum[idx])*delw[idx])

    # Random flux errors
    eqw_stat_err = \
        np.sqrt(np.sum((flux_err[idx]/continuum[idx]*delw[idx])**2))
    # Continuum errors
    eqw_cont_err = \
        np.sum(continuum_err[idx]*(flux[idx]/continuum[idx]**2)*delw[idx])
    # Zero point error
    # TODO: Check this calculation
    z_eps = 0.01
    eqw_zero_err = z_eps*eqw_int

    # Combine errors
    eqw_err = np.sqrt(eqw_stat_err**2 \
        +eqw_cont_err**2 + eqw_zero_err**2)

    spec['v1'] = integration_limits[0]
    spec['v2'] = integration_limits[1]

    # Store the EW in milliAngstrom
    spec['EW'] = eqw_int*1000.
    spec['EW_err'] = eqw_err*1000.
    spec['EW_err_stat'] = eqw_stat_err*1000.
    spec['EW_err_cont'] = eqw_cont_err*1000.
    spec['EW_err_zero'] = eqw_zero_err*1000.

    # Add the cumulative EW
    spec['EW_cumulative'] = \
      np.cumsum((1.-flux[idx]/continuum[idx])*delw[idx])*1000.

    # Calculate linear column density and error.
    linear_ncol = \
      ew_factor*spec['EW']/(spec['fval']*spec['wavc']**2)
    linear_ncol2sig = 2.0* \
      ew_factor*spec['EW_err']/(spec['fval']*spec['wavc']**2)
    linear_ncol3sig = 3.0* \
      ew_factor*spec['EW_err']/(spec['fval']*spec['wavc']**2)

    # Fill the output column densities
    spec['ncol_linearCoG'] = np.round(np.log10(linear_ncol),4)
    spec['ncol_linear2sig'] = \
        np.round(np.log10(linear_ncol2sig),4)
    spec['ncol_linear3sig'] = \
        np.round(np.log10(linear_ncol3sig),4)

    # Is the line detected at 2, 3 sigma?
    if spec['EW'] >= 2.*spec['EW_err']:
        spec['detection_2sig'] = True
    else:
        spec['detection_2sig'] = False

    if spec['EW'] >= 3.*spec['EW_err']:
        spec['detection_3sig'] = True
    else:
        spec['detection_3sig'] = False


    # Delete the old versions of the EW quantities
    try:
        del spec['w']
        del spec['w_es']
        del spec['w_ec']
        del spec['w_et']
        del spec['w_ez']
        del spec['col2sig']
        del spec['col3sig']
    except:
        pass

    return spec


def pyn_istat(spec_in,integration_limits = None):

    spec = spec_in.copy()

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)

    # Some constants and flags
    column_factor = 2.654e-15
    ew_factor = 1.13e17
    lightspeed = 2.998e5 # km/s

    velocity = spec['vel'].copy()
    flux = spec['flux'].copy()
    flux_err = spec['eflux'].copy()
    wavc=spec['wavc'].copy()
    fval=spec['fval'].copy()

    # Deal with the continuum
    if "contin" in spec.keys():
        continuum=spec['contin'].copy()
    else:
        try:
            continuum=spec['cont'].copy()
        except:
            continuum=spec['ycon'].copy()

    if "contin_err" in spec.keys():
        continuum_err = spec['contin_err'].copy()
    else:
        try:
            continuum_err = spec['econt'].copy()
        except:
            continuum_err = spec['ycon_sig'].copy()


    # Define the limits of the integration:
    if not integration_limits:
        integration_limits = [spec['v1'],spec['v2']]
        idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,[spec['v1'],spec['v2']])
        idx[xlim1:xlim2] = True
    else:
        idx = np.full(np.size(velocity),False)
        xlim1, xlim2 = \
            xlimit(velocity,integration_limits)
        idx[xlim1:xlim2] = True

    # An array of delta v:
    delv = velocity[1:]-velocity[:-1]
    delv = np.concatenate((delv,[delv[-1]]))

    # TODO: Saturation?
    # Calculate the zeroth moment
    tau = np.log(np.abs(continuum/flux))
    tau_tot = np.sum(tau[idx]*delv[idx])


    # M1
    # Calculate the first moment (average velocity)
    a = np.sum(tau[idx]*velocity[idx]*delv[idx])
    m1 = a/tau_tot

    # M1 error
    dadi = -1./flux[idx] * velocity[idx] * delv[idx]
    dadc = 1 / continuum[idx] * velocity[idx] * delv[idx]
    dwdi = -1./flux[idx] * delv[idx]
    dwdc = 1 / continuum[idx] * delv[idx]

    dm1di = (tau_tot * dadi - a * dwdi) / tau_tot**2
    dm1dc = (tau_tot * dadc - a * dwdc) / tau_tot**2

    q1 = np.sqrt(np.sum(flux_err[idx]**2 * dm1di**2))
    q2 = np.sum(np.sqrt(continuum_err[idx]**2 * dm1dc**2))

    m1err = np.sqrt(q1**2 + q2**2)


    # M2
    # Calculate the second moment (width)
    b = np.sum(tau[idx]*(velocity[idx] - m1)**2*delv[idx])
    m2 = np.sqrt(b/tau_tot)

    # M2 error
    dbdi = -1./flux[idx] * (velocity[idx] - m1)**2 * delv[idx]
    dbdc = 1./ continuum[idx] * (velocity[idx] - m1)**2 * delv[idx]
    dbdm1 = -2. * tau[idx] * (velocity[idx] - m1) * delv[idx]

    dm2di = (tau_tot * dbdi - b * dwdi) / tau_tot**2
    dm2dc = (tau_tot * dbdc - b * dwdc) / tau_tot**2
    dm2dm1 = dbdm1 / tau_tot

    q1 = np.sqrt(np.sum(flux_err[idx]**2 * dm2di**2))
    q2 = np.sum(np.sqrt(continuum_err[idx]**2 * dm2dc**2))
    q3 = np.sqrt(np.sum(m1err**2 * dm2dm1**2))

    m2err = np.sqrt(q1**2 + q2**2 + q3**2)
    m2err = m2err / (2. * m2)

    # M3
    # Calculate the third moment (skewness)
    c = np.sum(tau[idx]*((velocity[idx] - m1)/m2)**3*delv[idx])
    m3 = c/tau_tot

    # M3 error
    dfdi = -1./flux[idx] * ((velocity[idx] - m1) / m2)**3 * delv[idx]
    dfdc = 1./continuum[idx] * ((velocity[idx] - m1) / m2)**3 * delv[idx]
    dfdm1 = tau[idx]*3.*((velocity[idx] - m1) / m2)**2 * (-1./ m2) * delv[idx]
    dfdm2 = tau[idx]*3.*((velocity[idx] - m1) / m2)**2 * (m1 - velocity[idx]) / m2**2 * delv[idx]

    dm3di = (tau_tot * dfdi - c * dwdi) / tau_tot**2
    dm3dc = (tau_tot * dfdc - c * dwdc) / tau_tot**2
    dm3dm1 = dfdm1 / tau_tot
    dm3dm2 = dfdm2 / tau_tot

    q1 = np.sqrt(np.sum(flux_err[idx]**2 * dm3di**2))
    q2 = np.sum(np.sqrt(continuum_err[idx]**2 * dm3dc**2))
    q3 = np.sqrt(np.sum(m1err**2 * dm3dm1**2))
    q4 = np.sqrt(np.sum(m2err**2 * dm3dm2**2))

    m3err = np.sqrt(q1**2 + q2**2 + q3**2 + q4**2)


    # Calculate the extent (same as m2, except that m1 is assumed to be 0)
    b4 = np.sum(tau[idx]*(velocity[idx] - 0)**2*delv[idx])
    m4 = np.sqrt(b4/tau_tot)

    # M4 error
    dbdi = -1./flux[idx] * (velocity[idx] - 0.0)**2 * delv[idx]
    dbdc = 1./ continuum[idx] * (velocity[idx] - 0.0)**2 * delv[idx]

    dm4di = (tau_tot * dbdi - b4 * dwdi) / tau_tot**2
    dm4dc = (tau_tot * dbdc - b4 * dwdc) / tau_tot**2

    q1 = np.sqrt(np.sum(flux_err[idx]**2 * dm4di**2))
    q2 = np.sum(np.sqrt(continuum_err[idx]**2 * dm4dc**2))

    m4err = np.sqrt(q1**2 + q2**2)
    m4err = m4err / (2. * m4)


    # Velocities at 5% and 95% of total optical depth as dv90
    tau_cum = np.cumsum(tau[idx]*delv[idx])
    # 5% limit
    v90a = (np.abs(tau_cum/tau_tot-0.05)).argmin()
    v90a = velocity[idx][v90a]
    # 95% limit
    v90b = (np.abs(tau_cum/tau_tot-0.95)).argmin()
    v90b = velocity[idx][v90b]

    # Calculate dv90:
    dv90 = np.abs(v90b - v90a)


    # Fill the spec output
    spec['va'] = m1
    spec['va_err'] = m1err
    spec['ba'] = m2
    spec['ba_err'] = m2err
    spec['m3'] = m3
    spec['m3_err'] = m3err

    spec['dv90'] = dv90
    spec['v90a'] = v90a
    spec['v90b'] = v90b

    # Get rid of old-style keywords
    try:
        del spec['vaerr']
        del spec['baerr']
        del spec['m3err']
    except:
        pass

    return spec



def pyn_batch(spec_in,integration_limits = None, verbose = True):

    spec = spec_in.copy()

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)


    spec = pyn_column(spec,integration_limits)
    spec = pyn_eqwidth(spec,integration_limits)
    spec = pyn_istat(spec,integration_limits)

    dashes = '--------------------------------------------'
    if verbose:

        try:
            spec['ion'] = spec['ion'].decode('utf-8')
            spec['wni'] = spec['wni'].decode('utf-8')
        except:
            pass
        finally:
            print('********** '+spec['ion']+' '+spec['wni']+' **********')

        print('pyn_batch: Wavelength = {0:0.3f}'.format(spec['wavc']))
        print('pyn_batch: f-value = {0:0.3f}'.format(spec['fval']))
        print('\nVelocity range of integration: {0:0.1f} <= v <= {1:0.1f}'.format(spec['v1'],spec['v2']))

        print('\n'+dashes)
        # Print column densities:
        if spec['flag_sat']:
            print('***** WARNING: SATURATION IS PRESENT! *****')
            print(dashes)
            print('log N > {0:0.3f}'.format(spec['ncol']))
            print(dashes)
        else:
            print('log N = {0:0.3f} ({1:+0.3f}, {2:+0.3f})'.format(\
                    spec['ncol'],spec['ncol_err_lo'],spec['ncol_err_hi']))
        print(dashes)

        print('\n'+dashes)
        print('<v>       = {0:6.2f}  +/- {1:0.2f}'.format(spec['va'], spec['va_err']))
        print('<b>       = {0:6.2f}  +/- {1:0.2f}'.format(spec['ba'],spec['ba_err']))
        print('dv90      = {0:6.2f}  +/- {1:0.2f}'.format(spec['dv90'], spec['va_err']*np.sqrt(2)))
        print('Skew      = {0:6.2f}  +/- {1:0.2f}'.format(spec['m3'],spec['m3_err']))
        print(dashes)

        ew3sigma = 3.*spec['EW_err']

        print('\n'+dashes)
        print('EW           = {0:0.2f}'.format(spec['EW']))
        print('Stat Error   = {0:0.2f}'.format(spec['EW_err_stat']))
        print('Cont Error   = {0:0.2f}'.format(spec['EW_err_cont']))
        print('Tot Error    = {0:0.2f}'.format(spec['EW_err']))
        print('3sigma EW    < {0:0.2f}'.format(ew3sigma))
        if spec['EW'] < ew3sigma:
            print('***** WARNING: LINE NOT DETECTED AT 3 SIGMA! *****')

        print('\n')
        print('Linear COG N = {0:0.2f}'.format(spec['ncol_linearCoG']))
        print('2sigma N     < {0:0.2f}'.format(spec['ncol_linear2sig']))
        print('3sigma N     < {0:0.2f}'.format(spec['ncol_linear3sig']))
        print(dashes)

    return spec
