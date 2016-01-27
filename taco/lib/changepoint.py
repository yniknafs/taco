'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import numpy as np

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


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
    change_points = (np.ediff1d(x) != 0).nonzero()[0]
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
