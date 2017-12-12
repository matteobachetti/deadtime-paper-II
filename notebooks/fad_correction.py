import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from stingray.utils import rebin_data_log
from stingray.gti import bin_intervals_from_gtis
from stingray.gti import cross_two_gtis
from stingray import Powerspectrum, Crossspectrum, AveragedPowerspectrum, AveragedCrossspectrum
from hendrics.base import r_in, r_det
from hendrics.io import load_events
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import UnivariateSpline
import copy
import logging
import pytest


def _get_fourier_intv(lc, start_ind, end_ind):
    time = lc.time[start_ind:end_ind]
    counts = lc.counts[start_ind:end_ind]

    fourier = scipy.fftpack.fft(counts)

    freq = scipy.fftpack.fftfreq(len(time), lc.dt)
    good = freq > 0

    return freq[good], fourier[good], fourier[good] * np.sqrt(2 / (np.sum(counts)))


def FAD_power_spectra(lc1, lc2, segment_size, plot=False, smoothing_alg='gauss',
                      smoothing_length=None, verbose=False, tolerance=0.05,
                      strict=False, all_leahy=False):
    """Calculate Frequency Amplitude Difference-corrected (cross) power spectra.

    Parameters
    ----------
    lc1: class:`stingray.ligthtcurve.Lightcurve`
    lc1: class:`stingray.ligthtcurve.Lightcurve`
    segment_size: float
        Length of the segments to be averaged

    Other parameters
    ----------------
    plot : bool
        Plot diagnostics
    smoothing_alg : {'gauss', 'spline'}
        Smoothing algorithm
    smoothing_length : int
        Number of bins to smooth in gaussian window smoothing
    verbose: bool
        Print out information on the outcome of the algorithm (recommended)
    tolerance : float
        Accepted relative error on the FAD-corrected Fourier amplitude, to be
        used as success diagnostics.
        Should be
        ```
        stdtheor = 2 / np.sqrt(n)
        std = (average_corrected_fourier_diff / n).std()
        np.abs((std - stdtheor) / stdtheor) < tolerance
    strict : bool
        Fail if condition on tolerance is not met.

    Returns
    -------
    """
    if smoothing_length is None:
        smoothing_length = segment_size * 3
    gti = cross_two_gtis(lc1.gti, lc2.gti)
    lc1.gti = gti
    lc2.gti = gti
    lc1._apply_gtis()
    lc2._apply_gtis()
    summed_lc = copy.deepcopy(lc1)
    summed_lc.counts = lc1.counts + lc2.counts
    start_inds, end_inds = \
        bin_intervals_from_gtis(gti, segment_size, lc1.time,
                                dt=lc1.dt)
    pds1 = 0
    pds2 = 0
    ptot = 0
    cs = 0
    n = 0
    average_diff = average_diff_uncorr = 0

    if plot:
        plt.figure()

    for start_ind, end_ind in zip(start_inds, end_inds):
        freq, f1, f1_leahy = _get_fourier_intv(lc1, start_ind, end_ind)
        freq, f2, f2_leahy = _get_fourier_intv(lc2, start_ind, end_ind)
        freq, ftot, ftot_leahy = \
            _get_fourier_intv(summed_lc, start_ind, end_ind)
        
        fourier_diff = f1_leahy - f2_leahy

        if smoothing_alg == 'gauss':
            smooth_real = gaussian_filter1d(fourier_diff.real**2,
                                            smoothing_length)
        elif smoothing_alg == 'spline':
            spl = UnivariateSpline(freq, fourier_diff.real**2, s=1)
            smooth_real = spl(freq.astype(np.float64))
        if plot:
            plt.scatter(freq, fourier_diff, s=1)

        if all_leahy:
            f1 = f1_leahy
            f2 = f2_leahy
        p1 = (f1 * f1.conj()).real
        p1 = p1 / smooth_real * 2 
        p2 = (f2 * f2.conj()).real
        p2 = p2 / smooth_real * 2 
        pt = (ftot * ftot.conj()).real
        pt = pt / smooth_real * 2
        
        c = (f2 * f1.conj()).real
        c = c / smooth_real * 2

        if n == 0 and plot:
            plt.plot(freq, smooth_real, zorder=10, lw=3)
            plt.plot(freq, f1_leahy, zorder=5, lw=1)
            plt.plot(freq, f2_leahy, zorder=5, lw=1)

        ptot += pt
        pds1 += p1
        pds2 += p2
        cs += c
        average_diff += fourier_diff / smooth_real ** 0.5 * np.sqrt(2)
        average_diff_uncorr += fourier_diff
        n += 1

    std = (average_diff / n).std()
    stdtheor = 2 / np.sqrt(n)
    stduncorr = (average_diff_uncorr / n).std()
    is_compliant = np.abs((std - stdtheor) / stdtheor) < tolerance
    verbose_string = \
    '''
-------- FAD correction ----------
I smoothed over {smoothing_length} power spectral bins
{n} intervals averaged.
The uncorrected standard deviation of the Fourier
differences is {stduncorr} (dead-time affected!)
The final standard deviation of the FAD-corrected
Fourier differences is {std}. For the results to be
acceptable, this should be close to {stdtheor}
to within {tolerance} %.
In this case, the results ARE {compl}complying.
{additional}
----------------------------------
'''.format(smoothing_length=smoothing_length,
               n=n,
               stduncorr=stduncorr,
               std=std,
               stdtheor=stdtheor,
               tolerance=tolerance * 100,
               compl='NOT ' if not is_compliant else '',
               additional='Maybe something is not right.' if not is_compliant else '')
    logging.info(verbose_string)
    if strict:
        assert is_compliant
    if verbose:
        print(verbose_string)

    results = type('', (), {})()
    results.freq = freq
    results.pds1 = pds1 / n
    results.pds2 = pds2 / n
    results.cs = cs / n
    results.ptot = ptot / n
    results.is_compliant = is_compliant
    results.fad = average_diff / n
    results.fad_delta = (std - stdtheor) / stdtheor
    return results


@pytest.mark.parametrize('ctrate', [0.5, 5, 50, 500])
def test_fad_power_spectrum(ctrate):
    from stingray.lightcurve import Lightcurve
    from hendrics.fake import filter_for_deadtime
    dt = 0.1
    deadtime = 2.5e-3
    length = 51200
    time = np.arange(length) + dt / 2
    gti = np.array([[0, length]])
    segment_size = 512.
    ncounts = np.int(ctrate * length)
    ev1 = np.random.uniform(0, length, ncounts)
    ev2 = np.random.uniform(0, length, ncounts)
    ev1.sort()
    ev2.sort()
    ev1 = filter_for_deadtime(ev1, deadtime)
    ev2 = filter_for_deadtime(ev2, deadtime)
    ncounts1 = len(ev1)
    ncounts2 = len(ev2)

    lc1 = Lightcurve.make_lightcurve(ev1, dt=dt, tstart=0, tseg=length, gti=gti)
    lc2 = Lightcurve.make_lightcurve(ev2, dt=dt, tstart=0, tseg=length, gti=gti)

    lctot = lc1 + lc2

    results = \
        FAD_power_spectra(lc1, lc2, segment_size, plot=False,
                          smoothing_alg='gauss',
                          smoothing_length=segment_size*2,
                          strict=True, verbose=False,
                          tolerance=0.05)

    freq_f = results.freq
    pds1_f = results.pds1
    pds2_f = results.pds2
    cs_f = results.cs
    ptot_f = results.ptot
    is_compliant = results.is_compliant

    n = length / segment_size
    ncounts_per_intv1 = ncounts1 * segment_size / length
    ncounts_per_intv2 = ncounts2 * segment_size / length
    ncounts_per_intvtot = (ncounts1 + ncounts2) * segment_size / length
    ncounts_per_intv_geomav = np.sqrt(ncounts1 * ncounts2) * segment_size / length

    pds_std_theor = 2 / np.sqrt(n)
    cs_std_theor = np.sqrt(2 / n)

    assert np.isclose(pds1_f.std() * 2 / ncounts_per_intv1, pds_std_theor, rtol=0.05)
    assert np.isclose(pds2_f.std() * 2 / ncounts_per_intv2, pds_std_theor, rtol=0.05)
    assert np.isclose(cs_f.std() * 2 / ncounts_per_intv_geomav, cs_std_theor, rtol=0.05)
    assert np.isclose(ptot_f.std() * 2 / ncounts_per_intvtot, pds_std_theor, rtol=0.05)

    plt.show()
