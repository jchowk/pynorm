import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
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

def integration_weights(x,limits):
    # Calculate the weighting for each pixel in the column density integration.
    #  Includes partial pixel weighting for edge effects.

    # Pixel number array
    pix_num_array = np.arange(len(x))

    # Pixel spacing
    delx = np.median(np.roll(x,-1)-x)

    # Find the pixels that are fully within our integration range.
    idx = ((x-delx/2 >= limits[0]) & (x+delx/2 < limits[1]))
    weights = idx*1.0

    # Limits of complete pixels in indices.
    pix_limits = [np.min(pix_num_array[idx]),
                    np.max(pix_num_array[idx])]

    # Calculate fractions of edge pixels in the integration:
    # Identify edge pixels
    lo_pix = pix_limits[0]-1
    hi_pix = pix_limits[1]+1

    # Fraction of edge pixels contained within integration ranges
    lo_frac = ((x[lo_pix]+delx/2)-limits[0])/(delx)
    hi_frac = (limits[1]-(x[hi_pix]-delx/2))/(delx)

    # Assign fractional weights to the edge pixels
    weights[lo_pix] = lo_frac
    weights[hi_pix] = hi_frac

    return weights


def xlimit(x, limits):
    # Determine the integration limits for integer pixel integrations.

    def _ret():  return idx1, idx2

    # To be consistent with iNorm, we use x>=xmin, x<xmax
    idx1 = np.asarray(x >= limits[0]).nonzero()[0][0]
    idx2 = np.asarray(x > limits[1]).nonzero()[0][0]-1

    return _ret()

def integrate_column(velocity, flux, flux_err,
                continuum, continuum_err, wavc, fval,
                integration_limits = [-100,100],
                partial_pixels = True):
    # TODO: Asymmetrical error bars?

    # Some constants and flags
    column_factor = 2.654e-15
    flag_sat = False

    # Work out partial pixel weighting
    weights = integration_weights(velocity,integration_limits)
    # Uniform weighting if not partial pixel weighting
    if not partial_pixels:
        weights = np.zeros_like(velocity)
        xlim1, xlim2 = xlimit(velocity,integration_limits)
        weights[xlim1:xlim2+1] = 1.0

    # An array of delta v:
    delv = velocity[1:]-velocity[:-1]
    delv = np.concatenate((delv,[delv[-1]]))

    # Test for clearly saturated pixels:
    #   -- If the idx_saturation is already filled, use the results:
    try:
        idx_saturation = ((flux <= 0.) | (idx_saturation == True))
    except:
        idx_saturation = (flux <= 0.)

    # Fix saturation if it's present.
    flux[idx_saturation] = np.abs(flux[idx_saturation])
    flux[(flux==0)] = 2.*flux_err[(flux==0)]

    # Set overall saturation flag for saturation in the integration range
    if (idx_saturation*weights).sum() > 0:
        flag_sat = True

    # Create an optical depth array and its error
    tau_array = np.log(continuum / flux)
    tau_array_err = np.sqrt((flux_err/flux)**2)

    # If optical depth NaN, set to zero.
    bd = np.isnan(tau_array)
    tau_array[bd] = 0.
    tau_array_err[bd] = np.inf

    # Integrate the apparent optical depth
    tau_int = np.sum(tau_array*delv*weights)
    tau_int_err = \
     np.sqrt(np.sum((tau_array_err*delv*weights)**2))

    # Create an apparent column density array
    nav_array = tau_array/(wavc*fval*column_factor)
    nav_err = tau_array_err/(wavc*fval*column_factor)

    # Integrate the apparent column density profiles
    column = tau_int/(wavc*fval*column_factor)

    # Error in the column
    column_err = tau_int_err/(wavc*fval*column_factor)

    # Continuum error: errors are correlated, so don't add in quadrature.
    if partial_pixels:
        column_err_cont = \
            np.sum(((continuum_err/continuum)*delv*weights)) /\
             (wavc*fval*column_factor)
    else:
        column_err_cont = \
            np.sum(((continuum_err/continuum)*delv*weights)) /\
             (wavc*fval*column_factor)


    # Background uncertainty -- Not applied
    # z_eps = 0.01  # Fractional bg error
    # yc1 = continuum*weights*(1.-z_eps)
    # y1  = flux*weights-continuum*weights*z_eps
    # tau1 = np.sum(np.log(yc1/y1)*delv*weights)
    # col1 = tau1 / (wavc*fval*column_factor)
    # column_err_zero = np.abs(col1-column)

    # Combine errors
    column_err_total = np.sqrt(column_err**2 \
        +column_err_cont**2)

    return column, column_err_total, flag_sat


def pyn_column(spec_in, integration_limits = None,
                partial_pixels = True):

    spec = spec_in.copy()

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)

    # Make sure there are integration limits:
    if integration_limits is None:
        integration_limits = [spec['v1'],spec['v2']]

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

    # Work out partial pixel weighting
    weights = integration_weights(velocity,integration_limits)
    # Uniform weighting if not partial pixel weighting
    if not partial_pixels:
        weights = np.zeros_like(velocity)
        xlim1, xlim2 = xlimit(velocity,integration_limits)
        weights[xlim1:xlim2+1] = 1.0

    # An array of delta v:
    delv = velocity[1:]-velocity[:-1]
    delv = np.concatenate((delv,[delv[-1]]))

    # Test for clearly saturated pixels:
    #   -- If the idx_saturation is already filled, use the results:
    try:
        idx_saturation = ((flux <= 0.) | (idx_saturation == True))
    except:
        idx_saturation = (flux <= 0.)

    # Fix saturation if it's present.
    flux[idx_saturation] = np.abs(flux[idx_saturation])
    flux[(flux==0)] = 2.*flux_err[(flux==0)]

    # Set overall saturation flag for saturation in the integration range
    if (idx_saturation*weights).sum() > 0:
        flag_sat = True

    # Create an optical depth array and its error
    tau_array = np.log(continuum / flux)
    tau_array_err = np.sqrt((flux_err/flux)**2)

    # If optical depth NaN, set to zero.
    # This happens when continuum < 0.
    bd = np.isnan(tau_array)
    tau_array[bd] = 0.
    tau_array_err[bd] = 0.

    # TODO: Include an AOD array in output w/continuum errors.
    # Integrate the apparent optical depth
    tau_int = np.sum(tau_array*delv*weights)
    tau_int_err = \
     np.sqrt(np.sum((tau_array_err*delv*weights)**2))

    # Create an apparent column density array
    nav_array = tau_array/(wavc*fval*column_factor)
    nav_err_stat = tau_array_err/(wavc*fval*column_factor)
    nav_err_cont = (continuum_err/continuum)/(wavc*fval*column_factor)
    nav_err_tot = np.sqrt(nav_err_stat**2 + nav_err_cont**2)

    # Integrate the apparent column density profiles
    column = tau_int/(wavc*fval*column_factor)

    # Error in the column
    column_err = tau_int_err/(wavc*fval*column_factor)

    # Continuum error: errors are correlated, so don't add in quadrature.
    column_err_cont = \
        np.sum(((continuum_err/continuum)*delv*weights)) /\
         (wavc*fval*column_factor)

    # Background uncertainty -- Not applied
    # z_eps = 0.01  # Fractional bg error
    # yc1 = continuum*weights*(1.-z_eps)
    # y1  = flux*weights-continuum*weights*z_eps
    # tau1 = np.sum(np.log(yc1/y1)*delv*weights)
    # col1 = tau1 / (wavc*fval*column_factor)
    # column_err_zero = np.abs(col1-column)

    # Combine errors
    column_err_total = np.sqrt(column_err**2 \
        +column_err_cont**2)

    log_n_err_lo = np.log10(column-column_err_total) - np.log10(column)
    log_n_err_hi = np.log10(column+column_err_total) - np.log10(column)

    spec['v1'] = integration_limits[0]
    spec['v2'] = integration_limits[1]
    spec['ncol'] = np.log10(column)
    # Symmetrical errors:
    # spec['ncol_err_lo'] = -column_err_total/column*np.log10(np.e)
    # spec['ncol_err_hi'] = column_err_total/column*np.log10(np.e)
    # Asymmetrical errors:
    spec['ncol_err_lo'] = log_n_err_lo
    spec['ncol_err_hi'] = log_n_err_hi

    spec['flag_sat'] = flag_sat

    # Fill the Na(v) arrays
    spec['Nav'] = nav_array
    spec['Nav_err'] = nav_err_tot
    spec['Nav_sat'] = idx_saturation

    # Add the weights of the pixel integrations
    spec['integration_weights'] = weights

    # ---> Do this in the EW code.
    # Is the line detected at 2, 3 sigma?
    # if column >= 2.*column_err_total:
    #     spec['detection_2sig'] = True
    # else:
    #     spec['detection_2sig'] = False
    #
    # if column >= 3.*column_err_total:
    #     spec['detection_3sig'] = True
    # else:
    #     spec['detection_3sig'] = False


    if 'efnorm' in spec.keys():
        spec['fnorm_err'] = spec['efnorm']
        spec['fnorm_err_contin'] = spec['efnorm']*0.
        spec['fnorm_err_stat'] = spec['efnorm']

        del spec['efnorm']
        del spec['efnorm1']
        del spec['efnorm2']


    try:
        del spec['ncole1']
        del spec['ncole2']
        del spec['ncolez']
    except:
        pass

    return spec


def pyn_eqwidth(spec_in,integration_limits = None,
                partial_pixels = True):

    spec = spec_in.copy()

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)

    # Make sure there are integration limits:
    if integration_limits is None:
        integration_limits = [spec['v1'],spec['v2']]

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

    # Work out partial pixel weighting
    weights = integration_weights(velocity,integration_limits)
    # Uniform weighting if not partial pixel weighting
    if not partial_pixels:
        weights = np.zeros_like(velocity)
        xlim1, xlim2 = xlimit(velocity,integration_limits)
        weights[xlim1:xlim2+1] = 1.0

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
    eqw_int = np.sum((1.-flux/continuum)*delw*weights)
    # Random flux errors
    eqw_stat_err = \
        np.sqrt(np.sum((flux_err/continuum*delw*weights)**2))
    # Continuum errors
    eqw_cont_err = \
        np.sum(continuum_err*(flux/continuum**2)*delw*weights)

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
      np.cumsum((1.-flux/continuum)*delw*weights)*1000.

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

    # Pixel weighting factors
    spec['integration_weights'] = weights

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


def pyn_istat(spec_in,integration_limits = None,
                partial_pixels = True):

    spec = spec_in.copy()
    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)

    # Make sure there are integration limits:
    if integration_limits is None:
        integration_limits = [spec['v1'],spec['v2']]

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

    # Work out partial pixel weighting
    weights = integration_weights(velocity,integration_limits)
    # Uniform weighting if not partial pixel weighting
    if not partial_pixels:
        weights = np.zeros_like(velocity)
        xlim1, xlim2 = xlimit(velocity,integration_limits)
        weights[xlim1:xlim2+1] = 1.0

    # An array of delta v:
    delv = velocity[1:]-velocity[:-1]
    delv = np.concatenate((delv,[delv[-1]]))

    # TODO: Saturation?
    # Calculate the zeroth moment
    tau = np.log(np.abs(continuum/flux))
    tau_tot = np.sum(tau*delv*weights)

    # M1
    # Calculate the first moment (average velocity)
    a = np.sum(tau*velocity*delv*weights)
    m1 = a/tau_tot

    # M1 error
    dadi = -1./flux*velocity*delv*weights
    dadc = 1 / continuum*velocity*delv*weights
    dwdi = -1./flux*delv*weights
    dwdc = 1 / continuum*delv*weights

    dm1di = (tau_tot * dadi - a * dwdi) / tau_tot**2
    dm1dc = (tau_tot * dadc - a * dwdc) / tau_tot**2

    q1 = np.sqrt(np.sum((flux_err*weights)**2 * dm1di**2))
    q2 = np.sum(np.sqrt((continuum_err*weights)**2 * dm1dc**2))

    m1err = np.sqrt(q1**2 + q2**2)

    # M2
    # Calculate the second moment (width)
    bsqared = np.sum(tau*(velocity-m1)**2*delv*weights)
    m2 = np.sqrt(bsqared/tau_tot)

    # M2 error
    dbdi = -1./flux*(velocity-m1)**2 * delv*weights
    dbdc = 1./ continuum*(velocity-m1)**2 * delv*weights
    dbdm1 = -2.*tau*(velocity-m1) * delv*weights

    dm2di = (tau_tot * dbdi - bsqared*dwdi) / tau_tot**2
    dm2dc = (tau_tot * dbdc - bsqared*dwdc) / tau_tot**2
    dm2dm1 = dbdm1 / tau_tot

    q1 = np.sqrt(np.sum((flux_err*weights)**2 * dm2di**2))
    q2 = np.sum(np.sqrt((continuum_err*weights)**2 * dm2dc**2))
    q3 = np.sqrt(np.sum(m1err**2 * dm2dm1**2))

    m2err = np.sqrt(q1**2 + q2**2 + q3**2)
    m2err = m2err / (2.*m2)

    bvalue = m2*np.sqrt(2.)
    bvalue_err = m2err*np.sqrt(2.)

    # M3
    # Calculate the third moment (skewness)
    c = np.sum(tau*weights*((velocity*weights - m1)/m2)**3*delv*weights)
    m3 = c/tau_tot

    # M3 error
    dfdi = -1./flux*((velocity*weights - m1) / m2)**3 * delv*weights
    dfdc = 1./continuum*((velocity*weights - m1) / m2)**3 * delv*weights
    dfdm1 = tau*weights*3.*((velocity*weights - m1) / m2)**2 * (-1./ m2) * delv*weights
    dfdm2 = tau*weights*3.*((velocity*weights - m1) / m2)**2 * (m1 - velocity*weights) / m2**2 * delv*weights

    dm3di = (tau_tot * dfdi - c * dwdi) / tau_tot**2
    dm3dc = (tau_tot * dfdc - c * dwdc) / tau_tot**2
    dm3dm1 = dfdm1 / tau_tot
    dm3dm2 = dfdm2 / tau_tot

    q1 = np.sqrt(np.sum((flux_err*weights)**2 * dm3di**2))
    q2 = np.sum(np.sqrt((continuum_err*weights)**2 * dm3dc**2))
    q3 = np.sqrt(np.sum(m1err**2 * dm3dm1**2))
    q4 = np.sqrt(np.sum(m2err**2 * dm3dm2**2))

    m3err = np.sqrt(q1**2 + q2**2 + q3**2 + q4**2)

    # Calculate the extent (same as m2, except that m1 is assumed to be 0)
    b4 = np.sum(tau*(velocity*weights - 0)**2*delv*weights)
    m4 = np.sqrt(b4/tau_tot)

    # M4 error
    dbdi = -1./flux*(velocity*weights - 0.0)**2 * delv*weights
    dbdc = 1./ continuum*(velocity*weights - 0.0)**2 * delv*weights

    dm4di = (tau_tot * dbdi - b4 * dwdi) / tau_tot**2
    dm4dc = (tau_tot * dbdc - b4 * dwdc) / tau_tot**2

    q1 = np.sqrt(np.sum((flux_err*weights)**2 * dm4di**2))
    q2 = np.sum(np.sqrt((continuum_err*weights)**2 * dm4dc**2))

    m4err = np.sqrt(q1**2 + q2**2)
    m4err = m4err / (2. * m4)


    # Velocities at 5% and 95% of total optical depth as dv90
    tau_cum = np.cumsum(tau*delv*weights)
    # 5% limit
    v90a = (np.abs(tau_cum/tau_tot-0.05)).argmin()
    v90a = velocity[v90a]
    # 95% limit
    v90b = (np.abs(tau_cum/tau_tot-0.95)).argmin()
    v90b = velocity[v90b]

    # Calculate dv90:
    dv90 = np.abs(v90b - v90a)


    # Fill the spec output
    spec['va'] = m1
    spec['va_err'] = m1err
    # spec['ba'] = m2
    # spec['ba_err'] = m2err
    spec['ba'] = bvalue
    spec['ba_err'] = bvalue_err
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
# Added by saloni 
# to turn off blemish_correction set it to False in both pyn_batch and read_rbcodes
def pyn_blemish(spec_in,blemish_correction):
    spec = spec_in.copy()
    if blemish_correction:
        print('***** Will correct for blemishes, if present. *****')
    for i in range(len(spec['vel'])):
        if (((spec['eflux'][i]==-1.)|(spec['eflux'][i]>1))):
            if ((spec['vel'][i]>=spec['v1'])&(spec['vel'][i]<=spec['v2'])):
                spec['flag_blemish']=True #check
    
            if blemish_correction:
                x = (spec['vel'][i-20:i+21]).tolist(); y = (spec['flux'][i-20:i+21]).tolist(); z = (spec['eflux'][i-20:i+21]).tolist()
                ind = np.where((np.array(y)==0.0)|(np.array(z)>1))[0]
                for k in sorted(ind,reverse=True):
                    del x[k]; del y[k]
                if (len(x)<2):                     # unable to interpolate if blemish is on the leftmost edge
                    continue
                if (x[len(x)-1]<spec['vel'][i]):   # does not correct for blemishes on the edges, skips it
                    continue
                else:
                    ff=interpolate.interp1d(x,y,fill_value='extrapolate')
                    if (ff(spec['vel'][i])>=0):
                        spec['flux'][i]=ff(spec['vel'][i])
                    else:
                        spec['flux'][i]=0.0
                    spec['eflux'][i]=-0.9
    return spec

def pyn_batch(spec_in,integration_limits = None,
                partial_pixels = True, blemish_correction=True,
                verbose = True):

    spec = spec_in.copy()

    # FIX NON-WRITEABLE ARRAYS due to discontiguous
    # memory in some readsav inputs
    if ~spec['vel'].flags.writeable:
        spec = fix_unwriteable_spec(spec)


    # Make sure there are integration limits:
    if integration_limits is None:
        integration_limits = [spec['v1'],spec['v2']]
    
    spec = pyn_blemish(spec,blemish_correction)
    spec = pyn_eqwidth(spec,integration_limits, partial_pixels)
    spec = pyn_column(spec,integration_limits, partial_pixels)
    spec = pyn_istat(spec,integration_limits, partial_pixels)
    
    dashes = '--------------------------------------------'
    if verbose:

        try:
            spec['ion'] = spec['ion'].decode('utf-8')
            spec['wni'] = spec['wni'].decode('utf-8')
        except:
            pass
        finally:
            # TODO: FIX ME TO LOOK FOR MISSING WNI
            print('********** '+spec['ion']+' '+spec['wni']+' **********')

        print('pyn_batch: Wavelength = {0:0.3f}'.format(spec['wavc']))
        # if spec['fval'] <= 1.e-3:
        #     print('pyn_batch: f-value = {0:0.2g}'.format(spec['fval']))
        # else:
        #     print('pyn_batch: f-value = {0:0.3f}'.format(spec['fval']))
        print('pyn_batch: f-value = {0:0.3e}'.format(spec['fval']))
        print('\nVelocity range of integration: {0:0.1f} <= v <= {1:0.1f}'.format(spec['v1'],spec['v2']))

        print('\n'+dashes)
        # Print column densities:
        if spec['flag_blemish']:
            print('***** WARNING: BLEMISHES PRESENT IN THE INTEGRATION RANGE! *****')

        if spec['flag_sat']:
            print('***** WARNING: SATURATION IS PRESENT! *****')
            print(dashes)
            print('log N > {0:0.3f} ({1:+0.3f}, {2:+0.3f})'.format(\
                    spec['ncol'],spec['ncol_err_lo'],spec['ncol_err_hi']))
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
        if not spec['detection_3sig']:
            print('***** WARNING: LINE NOT DETECTED AT 3 SIGMA! *****')

        print('\n')
        print('Linear COG N = {0:0.2f}'.format(spec['ncol_linearCoG']))
        print('2sigma N     < {0:0.2f}'.format(spec['ncol_linear2sig']))
        print('3sigma N     < {0:0.2f}'.format(spec['ncol_linear3sig']))
        print(dashes)

    return spec
