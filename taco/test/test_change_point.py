'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import numpy as np

from taco.lib.base import Strand
from taco.lib.dtypes import FLOAT_DTYPE
from taco.lib.cChangePoint import mse as mse_cython
from taco.lib.changepoint import mse as mse_python
from taco.lib.changepoint import smooth

import matplotlib.pyplot as plt


def test_mse():
    a = np.zeros(100, dtype=FLOAT_DTYPE)
    a[:50] = 10
    a[50:] = 0
    a += np.random.random(100)
    mse_min_c, mse_i_c = mse_cython(a)
    mse_min_py, mse_i_py = mse_python(a)
    assert mse_i_c == 50
    assert mse_i_c == mse_i_py
    assert np.allclose(mse_min_c, mse_min_py, atol=0.01)


def test_smooth():
    a = np.zeros(1000)
    a[:200] = 100
    a[200:300] = np.linspace(100, 0, num=100)
    a[20] = 8
    a[40] = 8
    a[600:700] = np.linspace(0, 500, num=100)
    a[700:] = 500
    a += 50 * np.random.random(1000)


    for wsize in xrange(11, 101):
        y = smooth(a, window_len=52, window='hanning')
        print len(a), len(y), wsize

    print 'a', len(a)
    print 'y', len(y)

    plt.figure(1)
    plt.subplot(411)
    plt.plot(a)
    plt.subplot(412)
    plt.plot(y)
    plt.subplot(413)
    plt.plot(np.gradient(y))
    plt.subplot(414)
    y = np.ediff1d(a)
    y = np.gradient(a)

    print 'ediff1a', len(y)

    y = smooth(y, window_len=52, window='flat')
    plt.plot(y)

    plt.show()



def get_data1():
    import h5py
    filename = '/Users/mkiyer/Documents/lab/projects/taco/ccle/run1/expression.h5'
    f = h5py.File(filename, 'r')
    chrom = 'chr17'
    start = 40283388
    end = 40287603
    ref_stop = 40284136
    a = f[chrom][Strand.POS, start:end]
    return a


def permute(a, nperms=10):
    d = np.ediff1d(a)
    positions = np.arange(a.shape[0])
    start_values = d[d > 0]
    stop_values = d[d < 0]

    e = np.zeros((nperms, a.shape[0]))
    for i in xrange(nperms):
        print i
        e[i, :] += a[0]
        start_pos = np.random.choice(positions, size=len(start_values), replace=False)
        for j in xrange(len(start_values)):
            e[i, :start_pos[j]] += start_values[j]
        stop_pos = np.random.choice(positions, size=len(stop_values), replace=False)
        for j in xrange(len(stop_values)):
            e[i, stop_pos[j]:] += stop_values[j]

    return e





# def test1():
#     a = get_data1()
#     a = a[866:]
#     mse_obs, mse_i = mse(a)
#
#     nperms = 100
#
#     mse_rand = np.zeros(nperms)
#
#     # plt.figure(1)
#     # plt.subplot(411)
#     # plt.plot(a, 'k')
#
#     e = permute(a, nperms)
#     for i in xrange(nperms):
#         mse1 = mse(e[i])[0]
#         mse2 = mse_fast(e[i])[0]
#         mse_rand[i] = mse2
#
#         #plt.plot(e[i], 'r--')
#         #print 'plot', i
#         print 'hello', i, mse_rand[i]
#
#     print percentileofscore(mse_rand, mse_obs)
#     print mse_obs, mse_i
#     print np.sort(mse_rand)
#     plt.show()


def cusum(a):
    amean = a.mean()
    S = np.empty_like(a)
    S[0] = 0
    Smin = 0
    Smax = 0
    Sdiff = 0
    Sm = 0
    Sm_i = 0
    for i in xrange(1, S.shape[0]):
        S[i] = S[i-1] + (a[i] - amean)
        Smax = max(Smax, S[i])
        Smin = min(Smin, S[i])
        if (Smax - Smin) > Sdiff:
            Sdiff = Smax - Smin
        if np.abs(S[i]) > np.abs(Sm):
            Sm = S[i]
            Sm_i = i
    return S, Sdiff, Sm_i





def plot_cusum(a, ref_pos=None):
    S, Sdiff, Sm_i = cusum(a)
    val, m = mse(a)

    plt.figure(1)
    plt.subplot(211)
    plt.plot(a, 'k')
    plt.subplot(212)
    plt.plot(S, 'r--')

    for i in xrange(10):
        b = a.copy()
        np.random.shuffle(b)
        Sb, Sbdiff, j = cusum(b)
        plt.plot(Sb)

    plt.axvline(x=m, linewidth=1, color='r')
    plt.axvline(x=Sm_i, linewidth=1, color='b')
    if ref_pos is not None:
        plt.axvline(x=ref_pos, linewidth=1, color='g', linestyle='--')

    plt.show()

# def test_blah():
#     import h5py
#     filename = '/Users/mkiyer/Documents/lab/projects/taco/ccle/run1/expression.h5'
#     f = h5py.File(filename, 'r')
#
#     chrom = 'chr17'
#     start = 40283388
#     end = 40287603
#     ref_stop = 40284136
#     a = f[chrom][Strand.POS, start:end]
#     print cusum(a)
#     print mse(a)
#     plot_cusum(a, ref_pos=ref_stop - start)
#
#     chrom = 'chr17'
#     start = 40301912
#     end = 40314155
#     ref_stop = 40304657
#     ref_start = 40309193
#     a = f[chrom][Strand.POS, start:end]
#     print cusum(a)
#     print mse(a)
#     plot_cusum(a, ref_pos=ref_stop - start)
#
#     chrom = 'chr17'
#     start = 45101380
#     end = 45112704
#     ref_stop = 45108966
#     #ref_start = 40309193
#     a = f[chrom][Strand.NEG, start:end]
#     print cusum(a)
#     print mse(a)
#     plot_cusum(a, ref_pos=ref_stop - start)

# def test_cusum2():
#     a = np.zeros(1000)
#     a[:200] = 100
#     a[200:300] = np.linspace(100, 0, num=100)
#     a[20] = 8
#     a[40] = 8
#     a[600:700] = np.linspace(0, 500, num=100)
#     a[700:] = 500
#     a += 50 * np.random.random(1000)
#     # plot_cusum(a)
