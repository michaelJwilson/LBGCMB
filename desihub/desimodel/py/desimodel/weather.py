# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
=================
desimodel.weather
=================

Model of the expected weather conditions at KPNO during the DESI survey.

To generate a random time series of expected FWHM seeing in arcsecs and
atmospheric transparency, use, for example::

    n = 10000
    dt = 300 # seconds
    t = np.arange(n) * dt
    gen = np.random.RandomState(seed=123)
    seeing = sample_seeing(n, dt_sec=dt, gen=gen)
    transp = sample_transp(n, dt_sec=dt, gen=gen)

The resulting arrays are randomly sampled from models of the 1D probability
density and 2-point power spectral density derived from MzLS observations.
See `DESI-doc-3087
<https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=3087>`__
for details.

Used by :mod:`surveysim.weather` for simulations of DESI observing and
survey strategy studies.
"""
from __future__ import print_function, division

import numpy as np

import scipy.interpolate
import scipy.special


def whiten_transforms_from_cdf(x, cdf):
    """
    Calculate a pair of transforms to whiten and unwhiten a distribution.

    The whitening transform is monotonic and invertible.

    Parameters
    ----------
    x : array
        1D array of non-decreasing values giving bin edges for the distribution
        to whiten and unwhiten.
    cdf : array
        1D array of non-decreasing values giving the cummulative probability
        density associated with each bin edge.  Does not need to be normalized.
        Must have the same length as x.

    Returns
    -------
    tuple
        Tuple (F,G) of callable objects that whiten y=F(x) and unwhiten x=G(y)
        samples x of the input distribution, so that y has a Gaussian
        distribution with zero mean and unit variance.
    """
    x = np.asarray(x)
    cdf = np.asarray(cdf)
    if x.shape != cdf.shape:
        raise ValueError('Input arrays must have same shape.')
    if len(x.shape) != 1:
        raise ValueError('Input arrays must be 1D.')
    if not np.all(np.diff(x) >= 0):
        raise ValueError('Values of x must be increasing.')
    if not np.all(np.diff(cdf) >= 0):
        raise ValueError('Values of cdf must be increasing.')
    # Normalize.
    cdf /= cdf[-1]
    # Use linear interpolation for the forward and inverse transforms between
    # the input range and Gaussian CDF values.
    args = dict(
        kind='linear', assume_sorted=True, copy=False, bounds_error=True)
    forward = scipy.interpolate.interp1d(x, cdf, **args)
    backward = scipy.interpolate.interp1d(cdf, x, **args)
    # Add wrappers to convert between CDF and PDF samples.
    root2 = np.sqrt(2)
    forward_transform = (
        lambda x: root2 * scipy.special.erfinv(2 * forward(x) - 1))
    inverse_transform = (
        lambda y: backward(0.5 * (1 + scipy.special.erf(y / root2))))
    return forward_transform, inverse_transform


def whiten_transforms(data, data_min=None, data_max=None):
    """Calculate a pair of transforms to whiten and unwhiten a distribution.

    Uses :func:`desimodel.weather.whiten_transforms_from_cdf`.

    Parameters
    ----------
    data : array
        1D array of samples from the distribution to whiten.
    data_min : float or None
        Clip the distribution to this minimum value, or at min(data) if None.
        Must be <= min(data).
    data_max : float or None
        Clip the distribution to this maximum value, or at max(data) if None.
        Must be >= max(data).

    Returns
    -------
    tuple
        See :func:`desimodel.weather.whiten_transforms_from_cdf`.
    """
    n_data = len(data)
    # Sort the input data with padding at each end for the min/max values.
    sorted_data = np.empty(shape=n_data + 2, dtype=data.dtype)
    sorted_data[1:-1] = np.sort(data)
    if data_min is None:
        sorted_data[0] = sorted_data[1]
    else:
        if data_min > sorted_data[1]:
            raise ValueError('data_min > min(data)')
        sorted_data[0] = data_min
    if data_max is None:
        sorted_data[-1] = sorted_data[-2]
    else:
        if data_max < sorted_data[-2]:
            raise ValueError('data_max < max(data)')
        sorted_data[-1] = data_max
    # Calculate the Gaussian CDF value associated with each input value in
    # sorted order. The pad values are associated with CDF = 0, 1 respectively.
    cdf = np.arange(n_data + 2) / (n_data + 1.)

    return whiten_transforms_from_cdf(sorted_data, cdf)


def _seeing_fit_model(x):
    """Evalute the fit to MzLS seeing described in DESI-doc-3087.
    """
    p = np.array([  0.07511146,   0.44276671,  23.02442192,  38.07691498])
    y = (1 + ((x - p[0]) / p[1]) ** 2) ** (-p[2]) * x ** p[3]
    return y / (y.sum() * np.gradient(x))


def get_seeing_pdf(median_seeing=1.1, max_seeing=2.5, n=250):
    """Return PDF of FWHM seeing for specified clipped median value.

    Note that this is atmospheric seeing, not delivered image quality.
    The reference wavelength for seeing values is 6355A, in the r band,
    and the observed wavelength dependence in Dey & Valdes is closer to
    ``lambda ** (-1/15)`` than the ``lambda ** (-1/5)`` predicted by
    Kolmogorov theory. See DESI-doc-3087 for details.

    Scales the clipped MzLS seeing PDF in order to achieve the requested
    median value.  Note that clipping is applied before scaling, so
    the output PDF is clipped at scale * max_seeing.

    Parameters
    ----------
    median_seeing : float
        Target FWHM seeing value in arcsec. Must be in the range [0.95, 1.30].
    max_seeing : float
        Calculate scaled median using unscaled values below this value.
    n : int
        Size of grid to use for tabulating the returned arrays.

    Returns
    -------
    tuple
        Tuple (fwhm, pdf) that tabulates pdf[fwhm]. Normalized so that
        ``np.sum(pdf * np.gradient(fwhm)) = 1``.
    """
    # Tabulate the nominal (scale=1) seeing PDF.
    fwhm = np.linspace(0., max_seeing, n)
    pdf = _seeing_fit_model(fwhm)
    pdf /= (pdf.sum() * np.gradient(fwhm))
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    # Tabulate the median as a function of FWHM scale.
    scale = np.linspace(0.9, 1.4, 11)
    median = np.empty_like(scale)
    for i, s in enumerate(scale):
        median[i] = np.interp(0.5, cdf, s * fwhm)
    if median_seeing < median[0] or median_seeing > median[-1]:
        raise ValueError('Requested median is outside allowed range.')
    # Interpolate to find the scale factor that gives the requested median.
    s = np.interp(median_seeing, median, scale)
    return fwhm * s, pdf / s


def sample_timeseries(x_grid, pdf_grid, psd, n_sample, dt_sec=180., gen=None):
    """Sample a time series specified by a power spectrum and 1D PDF.

    The PSD should describe the temporal correlations of whitened samples.
    Generated samples will then be unwhitened to recover the input 1D PDF.
    See DESI-doc-3087 for details.

    Uses :func:`whiten_transforms_from_cdf`.

    Parameters
    ----------
    x_grid : array
        1D array of N increasing grid values covering the parameter range
        to sample from.
    pdf_grid : array
        1D array of N increasing PDF values corresponding to each x_grid.
        Does not need to be normalized.
    psd : callable
        Function of frequency in 1/days that returns the power-spectral
        density of whitened temporal fluctations to sample from. Will only be
        called for positive frequencies.  Normalization does not matter.
    n_sample : int
        Number of equally spaced samples to generate.
    dt_sec : float
        Time interval between samples in seconds.
    gen : np.random.RandomState or None
        Provide an existing RandomState for full control of reproducible random
        numbers, or None for non-reproducible random numbers.
    """
    x_grid = np.array(x_grid)
    pdf_grid = np.array(pdf_grid)
    if not np.all(np.diff(x_grid) > 0):
        raise ValueError('x_grid values are not increasing.')
    if x_grid.shape != pdf_grid.shape:
        raise ValueError('x_grid and pdf_grid arrays have different shapes.')
    # Initialize random numbers if necessary.
    if gen is None:
        gen = np.random.RandomState()
    # Calculate the CDF.
    cdf_grid = np.cumsum(pdf_grid)
    cdf_grid /= cdf_grid[-1]
    # Calculate whitening / unwhitening transforms.
    whiten, unwhiten = whiten_transforms_from_cdf(x_grid, cdf_grid)
    # Build a linear grid of frequencies present in the Fourier transform
    # of the requested time series.  Frequency units are 1/day.
    dt_day = dt_sec / (24. * 3600.)
    df_day = 1. / (n_sample * dt_day)
    f_grid = np.arange(1 + (n_sample // 2)) * df_day
    # Tabulate the power spectral density at each frequency.  The PSD
    # describes seeing fluctuations that have been "whitened", i.e., mapped
    # via a non-linear monotonic transform to have unit Gaussian probability
    # density.
    psd_grid = np.empty_like(f_grid)
    psd_grid[1:] = psd(f_grid[1:])
    # Force the mean to zero.
    psd_grid[0] = 0.
    # Force the variance to one.
    psd_grid[1:] /= psd_grid[1:].sum() * df_day ** 2
    # Generate random whitened samples.
    n_psd = len(psd_grid)
    x_fft = np.ones(n_psd, dtype=complex)
    x_fft[1:-1].real = gen.normal(size=n_psd - 2)
    x_fft[1:-1].imag = gen.normal(size=n_psd - 2)
    x_fft *= np.sqrt(psd_grid) / (2 * dt_day)
    x_fft[0] *= np.sqrt(2)
    x = np.fft.irfft(x_fft, n_sample)
    # Un-whiten the samples to recover the desired 1D PDF.
    x_cdf = 0.5 * (1 + scipy.special.erf(x / np.sqrt(2)))
    return np.interp(x_cdf, cdf_grid, x_grid)


def _seeing_psd(freq):
    """Evaluate the 'chi-by-eye' fit of the seeing PSD described in
    DESI-doc-3087.
    """
    N, f0, a0, a1 = 8000, 0.10, 2.8, -1.1
    return (N * (freq/f0)**a0 / (1 + (freq/f0)**a0) *
            (freq/f0) ** a1 / (10 + (freq/f0) ** a1))


def sample_seeing(n_sample, dt_sec=180., median_seeing=1.1, max_seeing=2.5,
                  gen=None):
    """Generate a random time series of FWHM seeing values.

    See DESI-doc-3087 for details. Uses :func:`get_seeing_pdf`,
    :func:`_seeing_psd` and :func:`sample_timeseries`.

    Parameters
    ----------
    n_sample : int
        Number of equally spaced samples to generate.
    dt_sec : float
        Time interval between samples in seconds.
    median_seeing : float
        See :func:`get_seeing_pdf`.
    mex_seeing : float
        See :func:`get_seeing_pdf`.
    gen : np.random.RandomState or None
        Provide an existing RandomState for full control of reproducible random
        numbers, or None for non-reproducible random numbers.

    Returns
    -------
    array
        1D array of randomly generated samples.
    """
    fwhm_grid, pdf_grid = get_seeing_pdf(median_seeing, max_seeing)
    return sample_timeseries(
        fwhm_grid, pdf_grid, _seeing_psd, n_sample, dt_sec, gen)


_transp_pdf_cum = np.array([0.06,0.11,1.0])
_transp_pdf_powers = np.array([0., 2.5, 35.])


def get_transp_pdf(n=250):
    """Return PDF of atmospheric transparency.

    Derived from MzLS observations, but corrected for dust accumulation and
    measurement error.  See DESI-doc-3087 for details.

    Parameters
    ----------
    n : int
        Size of grid to use for tabulating the returned arrays.

    Returns
    -------
    tuple
        Tuple (transp, pdf) that tabulates pdf[transp]. Normalized so that
        ``np.sum(pdf * np.gradient(transp)) = 1``.
    """
    transp = np.linspace(0., 1., n)
    pdf = np.zeros_like(transp)
    last_c = 0.
    for c, p in zip(_transp_pdf_cum, _transp_pdf_powers):
        pdf += (c - last_c) * np.power(transp, p) * (p + 1)
        last_c = c
    pdf /= pdf.sum() * np.gradient(transp)
    return transp, pdf


def _transp_psd(freq):
    """Evaluate the 'chi-by-eye' fit of the transparency PSD described in
    DESI-doc-3087.
    """
    N, f0, a0, a1 = 500, 1.5, 0.0, -1.5
    return (N * (freq/f0)**a0 / (1 + (freq/f0)**a0) *
            (freq/f0) ** a1 / (1 + (freq/f0) ** a1))


def sample_transp(n_sample, dt_sec=180., gen=None):
    """Generate a random time series of atmospheric transparency values.

    See DESI-doc-3087 for details. Uses :func:`get_transp_pdf`,
    :func:`_transp_psd` and :func:`sample_timeseries`.

    Parameters
    ----------
    n_sample : int
        Number of equally spaced samples to generate.
    dt_sec : float
        Time interval between samples in seconds.
    gen : np.random.RandomState or None
        Provide an existing RandomState for full control of reproducible random
        numbers, or None for non-reproducible random numbers.

    Returns
    -------
    array
        1D array of randomly generated samples.
    """
    transp_grid, pdf_grid = get_transp_pdf()
    return sample_timeseries(
        transp_grid, pdf_grid, _transp_psd, n_sample, dt_sec, gen)
