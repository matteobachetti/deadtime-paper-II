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


def _get_fourier_intv(lc, start_ind, end_ind):
    time = lc.time[start_ind:end_ind]
    counts = lc.counts[start_ind:end_ind]

    fourier = scipy.fftpack.fft(counts)

    freq = scipy.fftpack.fftfreq(len(time), lc.dt)
    good = freq > 0
# #     print("Bu", factor, countrate_factor, counts[:10], counts_r_out[:10], fourier[:10])

    return freq[good], fourier[good], fourier[good] * np.sqrt(2 / (np.sum(counts)))


def FAD_power_spectra(lc1, lc2, segment_size, plot=False, smoothing_alg='gauss'):
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

    if plot:
        plt.figure()

    for start_ind, end_ind in zip(start_inds, end_inds):
        freq, f1, f1_leahy = _get_fourier_intv(lc1, start_ind, end_ind)
        freq, f2, f2_leahy = _get_fourier_intv(lc2, start_ind, end_ind)
        freq, ftot, ftot_leahy = _get_fourier_intv(summed_lc, start_ind, end_ind)
        
        fourier_diff = f1_leahy - f2_leahy
        if smoothing_alg == 'gauss':
            smooth_real = gaussian_filter1d(fourier_diff.real**2, segment_size * 3)
        elif smoothing_alg == 'spline':
            spl = UnivariateSpline(freq, fourier_diff.real**2, s=1)
            smooth_real = spl(freq.astype(np.float64))
        elif smoothing_alg == 'savgol':
            spl = UnivariateSpline(freq, fourier_diff.real**2, s=1)
            smooth_real = spl(freq.astype(np.float64))
        if plot:
            plt.scatter(freq, fourier_diff, s=1)

        p1 = (f1 * f1.conj()).real
        p1 = p1 / smooth_real * 2 
        p2 = (f2 * f2.conj()).real
        p2 = p2 / smooth_real * 2 
        pt = (ftot * ftot.conj()).real
        pt = pt / smooth_real * 2
        
        c = (f2 * f1.conj()).real
        c = c / smooth_real * 2 

        if n == 0 and plot:

#           plt.plot(freq, fourier_diff, zorder=10, lw=1)
            plt.plot(freq, smooth_real, zorder=10, lw=3)

            plt.plot(freq, f1_leahy, zorder=5, lw=1)
            plt.plot(freq, f2_leahy, zorder=5, lw=1)

        ptot += pt
        pds1 += p1
        pds2 += p2
        cs += c
        n += 1
        
#         break
    
    return freq, pds1 / n, pds2 / n, cs / n, ptot / n

