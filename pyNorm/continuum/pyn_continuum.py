# TODO: Switch to scipy F-test
# TODO: Switch to numpy Legendre fitting / error analysis

def pyn_ftest(nu, p):
    """
    This function calculates the significance of a fit using the
     	"f-test".  See Bevington (1969).
    """
    import numpy as np

    n_params = 2
    def _ret():  return None

    f_array = np.zeros([20, 8], dtype="float32")
    i = np.array([[0.50, 0.25, 0.10, 0.05, 0.025, 0.01, 0.005, 0.001]])
    j = np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20, 24, 30, 40, 60, 120, 1.e5]])

    f_array[0,:] = np.array([[1.00, 5.83, 39.9, 161., 648., 4050., 16200., 406000.]])
    f_array[1,:] = np.array([[0.667, 2.57, 8.53, 18.5, 38.5, 98.5, 198., 998.]])
    f_array[2,:] = np.array([[0.585, 2.02, 5.54, 10.1, 17.4, 34.1, 55.6, 167.]])
    f_array[3,:] = np.array([[0.549, 1.81, 4.54, 7.71, 12.2, 21.2, 31.3, 74.1]])
    f_array[4,:] = np.array([[0.528, 1.69, 4.06, 6.61, 10.0, 16.3, 22.8, 47.2]])

    f_array[5,:] = np.array([[0.515, 1.62, 3.78, 5.99, 8.81, 13.7, 18.6, 35.5]])
    f_array[6,:] = np.array([[0.506, 1.57, 3.59, 5.59, 8.07, 12.2, 16.2, 29.2]])
    f_array[7,:] = np.array([[0.499, 1.54, 3.46, 5.32, 7.57, 11.3, 14.7, 25.4]])
    f_array[8,:] = np.array([[0.494, 1.51, 3.36, 5.12, 7.21, 10.6, 13.6, 22.9]])
    f_array[9,:] = np.array([[0.490, 1.49, 3.28, 4.96, 6.94, 10.0, 12.8, 21.0]])

    f_array[10,:] = np.array([[0.486, 1.47, 3.23, 4.84, 6.72, 9.65, 12.2, 19.7]])
    f_array[11,:] = np.array([[0.484, 1.46, 3.18, 4.75, 6.55, 9.33, 11.8, 18.6]])
    f_array[12,:] = np.array([[0.478, 1.43, 3.07, 4.54, 6.20, 8.68, 10.8, 16.6]])
    f_array[13,:] = np.array([[0.472, 1.40, 2.97, 4.35, 5.87, 8.10, 9.94, 14.8]])
    f_array[14,:] = np.array([[0.469, 1.39, 2.93, 4.26, 5.72, 7.82, 9.55, 14.0]])

    f_array[15,:] = np.array([[0.466, 1.38, 2.88, 4.17, 5.57, 7.56, 9.18, 13.3]])
    f_array[16,:] = np.array([[0.463, 1.36, 2.84, 4.08, 5.42, 7.31, 8.83, 12.6]])
    f_array[17,:] = np.array([[0.461, 1.35, 2.79, 4.00, 5.29, 7.08, 8.49, 12.0]])
    f_array[18,:] = np.array([[0.458, 1.34, 2.75, 3.92, 5.15, 6.85, 8.18, 11.4]])
    f_array[19,:] = np.array([[0.455, 1.32, 2.71, 3.84, 5.02, 6.63, 7.88, 10.8]])

    jj = np.where(np.ravel(j <= nu))[0]
    jj = jj[jj.size - 1]

    ii = np.where(np.ravel(i <= p))[0]
    ii = ii[0]

    return f_array[jj,ii]

def legbasis(x, maxord):
    import numpy as np

    numpix = x.size

    #Form legendre polynomial.
    p = np.zeros([maxord + 1, numpix])
    p[0,:] = 1.
    p[1,:] = x
    for j in np.arange(2, maxord+1):
        p[j,:] = ((2.*j- 1.)*x*p[j-1,:] - (j-1)*p[j-2,:])/j

    return p

def legfit(x,y,nord):
    import numpy as np

    numpix = x.size
    p = legbasis(x, nord)

    ncoeff = nord + 1

    #Form alpha and beta matrices.
    beta = np.zeros([nord + 1])
    alpha = np.zeros([nord + 1, nord + 1], "float32")
    # Fill alpha and beta matrices
    for k in np.arange(0, nord+1):
        beta[k] = np.sum(y * p[k,:])
    for k in np.arange(0, nord+1):
        for j in np.arange(0, nord+1):
            alpha[j, k] = np.sum(p[j,:] * p[k,:])

    # Invert alpha matrix ==> error matrix eps.
    eps = np.linalg.inv(alpha)

    #Calculate coefficients and fit.
    a = np.transpose(np.matmul(np.transpose(beta), np.transpose(eps)))
    yfit = np.zeros([numpix], "float32")

    # Sum the Legendre orders
    for j in np.arange(0, nord+1):
        yfit += a[j]*p[j,:]

    return yfit, a

def legpoly(x,coeff):
    import numpy as np

    numpix = x.size
    nord = coeff.size - 1

    p = legbasis(x, nord)
    yfit = np.zeros_like(x)

    # Sum the Legendre orders
    for j in np.arange(0, nord+1):
        yfit += coeff[j]*p[j,:]

    return yfit


def legerr(x,y,a):
    """Warning: assumes uniform weighting for pixels.
    """
    import numpy as np

    ncoeff = a.size
    nord = a.size - 1

    numpix = x.size
    ix = np.arange(numpix,dtype="float32")

    p = legbasis(x,nord)
    yfit = legpoly(x,a)

    # Degrees of freedum
    nu = numpix - ncoeff - 1
    # Variance in f-test parlance
    variance = np.sum((y - yfit)**2)
    #Calculate chi squared of fit - uniform weighting=1.
    chi2 = variance / (numpix - ncoeff - 1)

    #Form alpha and beta matrices.
    beta = np.zeros([nord + 1])
    alpha = np.zeros([nord + 1, nord + 1], dtype="float32")
    # Fill alpha and beta matrices
    for k in np.arange(0, nord+1):
        beta[k] = np.sum(y * p[k,:])
    for k in np.arange(0, nord+1):
        for j in np.arange(0, nord+1):
            alpha[j, k] = np.sum(p[j,:] * p[k,:])

    # Invert alpha matrix ==> error matrix eps.
    eps = np.linalg.inv(alpha)

    eps1 = chi2 * eps
    tot = np.zeros([numpix], dtype="float32")
    for i in np.arange(0, (numpix - 1)+(1)):
        tot[i] = 0
        for l in np.arange(0, nord+1):
            for k in np.arange(0, nord+1):
                tot[i] += eps1[k,l] * p[l,i] * p[k,i]
    error = np.sqrt(tot)

    return error

def convert_inorm_mask(spec_in):
    import numpy as np

    spec = spec_in.copy()

    # Extract the mask
    try:
        mask = spec['mask_cont']
    except:
        return spec

    # Create a boolean mask from original
    gdcont = (mask == 1)
    # Find the places where the mask changes.
    delta = mask-np.roll(mask,1)

    # Identify where the mask transitions are:
    vstarts = np.where(delta == -1)[0]
    vstops = np.where(delta == 1)[0]
    # If the first data point is fitted, adjust the continuum boundaries.
    if mask[0] == 1:
        strt = vstops
        stps = vstarts
        vstarts = np.concatenate([np.array([0]), strt])
        vstops = np.concatenate([stps, np.array([np.size(mask)-1])])
    # Transform the results to velocity ranges:
    vstarts = spec['vel'][vstarts]
    vstops = spec['vel'][vstops]

    # # Create a boolean mask from velocity ranges
    # gdcont = np.repeat(False,np.size(spec['vel']))
    # for j in np.arange(np.size(vstarts)):
    #     gdnew = (spec['vel'] >= vstarts[j]) & (spec['vel'] <= vstops[j])
    #     gdcont = gdcont | gdnew

    spec['contin_mask_bits'] = mask
    spec['contin_mask_bool'] = gdcont
    spec['contin_v1'] = vstarts
    spec['contin_v2'] = vstops

    try:
        del spec['mask_cont']
    except:
        pass

    return spec

def continuum_fit(spec_in, minord, maxord):
    """
    Fit Legendre polynomial continua
    """
    import numpy as np
    import scipy

    spec = spec_in.copy()

    # Create new-style continuum mask
    if "mask_cont" in spec:
        spec = convert_inorm_mask(spec)
    # Which are the continuum regions
    gd = spec['contin_mask_bool']

    # Misc.
    nflag = 0
    nord = maxord

    # Set the variables for the Legendre fit
    xmax = np.max(np.abs(spec['vel'][gd]))
    x = spec['vel'][gd]/xmax
    y = spec['flux'][gd]

    # Array subscript length and vector.
    numpix = x.size
    idx = np.arange(numpix)

    # Begin loop to do fit for each order.
    for nord in np.arange(minord, maxord+1):
        ncoeff = nord + 1
        yfit, coeff = legfit(x,y,nord)

        # Variance in f-test parlance
        variance = np.sum((y - yfit)**2)
        # Degrees of freedum
        nu = numpix - ncoeff - 1
        #Calculate chi squared of fit - uniform weighting=1.
        chi2 = variance / (numpix - ncoeff - 1)

        #Check chi2 against previous chi2 to see if additional term should
        #be kept.  Check probability for "95% confidence".
        if nord > minord:
            # p_value = scipy.stats.f.cdf(F, df1, df2)
            # F statistic
            f = (variance1 - variance) / chi2

            fcutoff = pyn_ftest(nu, 0.05)
            if f < fcutoff:
                nord = nord - 1
                yfit, coeff = legfit(x,y,nord)
                nflag = 1
                break

        variance1 = variance

    spec['contin'] = legpoly(spec['vel']/xmax,coeff)
    spec['contin_err'] = legerr(spec['vel']/xmax,spec['flux'],coeff)
    spec['contin_order'] = nord
    spec['contin_coeff'] = coeff

    # Normalized continuum
    spec['vnorm'] = spec['vel']
    spec['fnorm'] = spec['flux'] / spec['contin']
    # Continuum error in fnorm
    spec['fnorm_err_contin'] = spec['contin_err'] / spec['contin']
    # Statistical error in fnorm
    try:
        ferr_stat = spec['eflux']
    except:
        ferr_stat = spec['ef']
    spec['fnorm_err_stat'] = ferr_stat / spec['contin']
    # Total error in fnorm
    spec['fnorm_err'] = np.sqrt(spec['fnorm_err_contin']**2 + spec['fnorm_err_stat']**2)
    # Estimate signal-to-noise ratio in continuum regions
    spec['SNR'] = np.median(spec['fnorm'][gd]/spec['fnorm_err'][gd])

    # Remove old-style continuum / normalized flux information
    try:
        del spec['ycon']
        del spec['ycon_sig']
        del spec['efnorm']
        del spec['efnorm1']
        del spec['efnorm2']
        del spec['maxord']
        del spec['sn']
    except:
        pass

    return spec


def continuum_autofit(spec_in,
                vshift=0.,vcont = [-1000,1000],
                vclip = 5.):

    import numpy as np

    spec = spec_in.copy()

    # Make sure velocities are in the right order
    if vcont[0] > vcont[1]: np.reverse(vcont)

    # TODO: Enable a frame of reference flag [LSR v Helio]
    # Shift the velocity.
    spec['vel'] += vshift
    spec['v1'] += vshift
    spec['v2'] += vshift

    # Define fitting region for autocontinuum
    gd = (spec['vel'] >= vcont[0]) & (spec['vel'] <= vcont[1])

    # Extract initial continuum fit region
    vi1 = spec['vel'][gd]
    fi1 = spec['flux'][gd]
    ei1 = spec['eflux'][gd]

    # Set the variables for the Legendre fit
    xmax = np.max(np.abs(spec['vel'][gd]))
    x = spec['vel'][gd]/xmax
    y = spec['flux'][gd]

    # New mask variable
    mask_cont = gd.copy()


    # Cut the spectrum in sub-sections to check regions to be masked.
    # This should help when flux changes by a large amount over the considered
    # velocity interval.
    diffv = (vl2-vl1)/vclip1
    diffv1  = vl1
    diffv2   = vl1 + diffv




    return spec
