'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import os

import h5py
import numpy as np
from scipy.stats import mannwhitneyu

from taco.lib.base import Strand, TacoError
from taco.lib.cChangePoint import mse as mse_cython


__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def run_changepoint(a, size_cutoff=20, cp_func=mse_cython,
                    smooth_window="hanning", window_len=75):
    '''
    a: array with expression data
    size_cutoff: minimum length of interval to allow recursive change point
                 searches
    cp_func: distance function to compute change point location
    smooth_window: numpy smoothing window type
    window_len: size of smoothing window

    returns list of changepoints where each is a tuple (i, p, j, k, sign)
        i: index of change point within array 'a'
        p: p-value of change point (mannwhitneyu)
        j: distance from i to left slope boundary
        k: distance from i to right slope boundary
        sign: direction of change
    '''
    if a.shape[0] < window_len:
        return []
    s_a = np.gradient(smooth(a, window_len=window_len, window=smooth_window))
    cps = bin_seg_slope(a, s_a, cp_func, size_cutoff=size_cutoff)
    return cps


def get_data(rundir, chrom, start, end, strand):
    filename = os.path.join(rundir, 'expression.h5')
    f = h5py.File(filename, 'r')
    if strand == '+':
        strand_idx = Strand.POS
    elif strand == '-':
        strand_idx = Strand.NEG
    else:
        strand_idx = Strand.NA
    a = f[chrom][strand_idx, start:end]
    return a


def slope_extract(slope_a, i):
    slope_a = np.sign(slope_a)
    sign = slope_a[i]
    if slope_a[i] == 0:
        j = 0
        k = 0
    else:
        j = 0
        while (i - j) >= 0 and slope_a[i - j] == slope_a[i]:
            j += 1
        k = 0
        while (i + k) < len(slope_a) and slope_a[i + k] == slope_a[i]:
            k += 1
    return (j, k, sign)


def mwu_ediff(a, i):
    a1 = a[:i]
    a2 = a[i:]
    a1 = a1[np.ediff1d(a1).nonzero()[0]]
    a2 = a2[np.ediff1d(a2).nonzero()[0]]
    if len(a1) == 0 or len(a2) == 0:
        return (None, 1)
    elif (len(np.unique(a1)) == 1 and
          np.array_equal(np.unique(a1), np.unique(a2))):
        return (None, 1)
    else:
        U, p = mannwhitneyu(a1, a2)
        return (U, p)


def bin_seg_slope(a, s_a, cp_func=mse_cython, PVAL=0.05, cps=None, offset=0,
                  size_cutoff=20):
    '''
    a: expr data vector
    s_a: slope vector
    cp_func: function to choose change point
    PVAL: threshold to call change points
    cps: (recursion) list of change points
    offset: (recursion) offset into vectors
    size_cutoff: stop searching for change points when vector
                 length < size_cutoff
    '''
    if a.shape[0] < size_cutoff:
        return cps
    if cps is None:
        cps = []
    stat, i = cp_func(a)
    U, p = mwu_ediff(a, i)
    if i <= 1:
        return cps
    elif p < PVAL:
        j, k, sign = slope_extract(s_a, offset + i)
        if j != 0 and k != 0:
            # save changepoint
            print 'CHANGE', i, p, j, k, sign
            cps.append((i + offset, p, j, k, sign))
            # test left
            if (offset+i-j) > offset:
                b1 = a[:i-j]
                cps = bin_seg_slope(b1, s_a, cp_func, cps=cps, offset=offset)
            # test right
            if (offset+i+k) < offset+len(a):
                b2 = a[i+k:]
                cps = bin_seg_slope(b2, s_a, cp_func, cps=cps,
                                    offset=(offset + i + k))
    return cps


def find_change_points(a, start=0, threshold=0):
    '''
    a - numpy array of expression data (all values >= 0)
    start - genomic start of 'a'
    threshold - expression threshold to detect changes
    '''
    if a.shape[0] == 0:
        return []
    if a.shape[0] == 1:
        return []
    change_points = []
    for i in xrange(1, a.shape[0]):
        if (a[i-1] > threshold) and (a[i] <= threshold):
            change_points.append(start + i)
        elif (a[i-1] <= threshold) and (a[i] > threshold):
            change_points.append(start + i)
    return change_points


def mse(x):
    change_points = (np.ediff1d(x) != 0).nonzero()[0] + 1
    if len(change_points) == 0:
        return 0.0, -1

    mse_min = None
    mse_i = -1
    for i in change_points:
        x1 = x[:i]
        x2 = x[i:]
        mse1 = np.power(x1 - x1.mean(), 2).sum()
        mse2 = np.power(x2 - x2.mean(), 2).sum()
        mse = mse1 + mse2

        if (mse_i == -1) or (mse < mse_min):
            mse_min = mse
            mse_i = i
    return mse_min, mse_i


def smooth(x, window_len=11, window='hanning'):
    """
    http://scipy-cookbook.readthedocs.org/items/SignalSmooth.html

    smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd
                    integer
        window: the type of window from 'flat', 'hanning', 'hamming',
                'bartlett', 'blackman', flat window will produce a moving
                average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve, scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead
          of a string
    NOTE: length(output) != length(input), to correct this:
          return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is not one of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    # moving average
    if window == 'flat':
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y[(window_len/2-1):-(window_len/2)]
