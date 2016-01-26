'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import numpy as np

from taco.lib.base import Strand
import matplotlib.pyplot as plt


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


def mse(x):
    mse = None
    mse_m = None
    for i in xrange(1, x.shape[0]):
        x1 = x[:i]
        x2 = x[i:]

        mse1 = np.power(x1 - x1.mean(), 2).sum()
        mse2 = np.power(x2 - x2.mean(), 2).sum()

        if (mse is None) or (mse1 + mse2 < mse):
            mse = mse1 + mse2
            mse_m = i
    return mse, mse_m


def plot_cusum(a, ref_pos=None):
    return

    #a = np.ediff1d(a)
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

def test_cusum2():
    a = np.zeros(1000)
    a[:200] = 100
    a[200:300] = np.linspace(100, 0, num=100)
    a[20] = 8
    a[40] = 8
    a[600:700] = np.linspace(0, 500, num=100)
    a[700:] = 500
    a += 50 * np.random.random(1000)
    # plot_cusum(a)
