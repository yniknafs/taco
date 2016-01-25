'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import numpy as np

from taco.lib.base import TacoError, Strand
from taco.lib.dtypes import FLOAT_DTYPE
from taco.lib.transfrag import Transfrag
from taco.lib.cChangePoint import find_change_points as c_find_change_points
from taco.lib.changepoint import find_change_points
from taco.lib.locus import Locus, StrandedLocus, split_transfrag

from taco.test.base import read_gtf, read_single_locus


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
    import matplotlib.pyplot as plt

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


def test_blah():
    import h5py
    filename = '/Users/mkiyer/Documents/lab/projects/taco/ccle/run1/expression.h5'
    f = h5py.File(filename, 'r')

    chrom = 'chr17'
    start = 40283388
    end = 40287603
    ref_stop = 40284136
    a = f[chrom][Strand.POS, start:end]
    print cusum(a)
    print mse(a)
    plot_cusum(a, ref_pos=ref_stop - start)

    chrom = 'chr17'
    start = 40301912
    end = 40314155
    ref_stop = 40304657
    ref_start = 40309193
    a = f[chrom][Strand.POS, start:end]
    print cusum(a)
    print mse(a)
    plot_cusum(a, ref_pos=ref_stop - start)

    chrom = 'chr17'
    start = 45101380
    end = 45112704
    ref_stop = 45108966
    #ref_start = 40309193
    a = f[chrom][Strand.NEG, start:end]
    print cusum(a)
    print mse(a)
    plot_cusum(a, ref_pos=ref_stop - start)





# def test_cusum():
#     a = np.array([10.7, 13.0, 11.4, 11.5, 12.5, 14.1, 14.8, 14.1, 12.6,
#                   16.0, 11.7, 10.6, 10.0, 11.4, 7.9, 9.5, 8.0, 11.8, 10.5,
#                   11.2, 9.2, 10.1, 10.4, 10.5])
#     #plot_cusum(a)


# def test_cusum2():
#     a = np.zeros(100)
#     a[50:] = 10
#     a[20] = 8
#     a[40] = 8
#     a += np.random.random(100)
#     plot_cusum(a)


# def test_find_boundaries():
#     t_dict, locus = read_single_locus('splice_sites.gtf')
#     transfrags = t_dict.values()
#     splice_sites = set()
#     for t in transfrags:
#         splice_sites.update(t.itersplices())
#     splice_sites = tuple(sorted(splice_sites))
#     assert splice_sites == (100, 200, 250, 300, 400)
#     # aggregate expression
#     slocus = StrandedLocus.create(transfrags)
#     # zero change points
#     zero_sites = tuple(find_change_points(slocus.expr_data, slocus.start))
#     assert zero_sites == (100, 150, 300, 375)
#     # combined boundaries
#     boundaries = tuple(slocus._find_node_boundaries())
#     assert boundaries == (10, 100, 150, 200, 250, 300, 375, 400, 525)
