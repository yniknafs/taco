'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt

from taco.lib.base import Strand
from taco.lib.dtypes import FLOAT_DTYPE
from taco.lib.cChangePoint import mse as mse_cython
from taco.lib.changepoint import mse as mse_python
from taco.lib.changepoint import smooth, run_changepoint


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
    return


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


def plot_slope(a, ref_pos=None):
    import matplotlib.pyplot as plt
    smooth_a = np.array(smooth(a), dtype=FLOAT_DTYPE)
    slope_a = np.gradient(smooth_a)

    plt.figure(1)
    plt.subplot(311)
    plt.plot(a, 'k')
    plt.subplot(312)
    plt.plot(smooth_a, 'g')

    cps = run_changepoint(a)
    print cps
    for cp in cps:
        m, p, j, k, sign = cp
        sign = np.sign(slope_a[m])
        color='k'
        if sign == 1: color = 'r'
        if sign == -1: color = 'b'
        plt.axvline(x=m, linewidth=1, color=color)
        for ref in ref_pos:
            plt.axvline(x=ref, linewidth=1, color='m')
    plt.subplot(313)
    plt.plot(slope_a, 'g')
    for cp in cps:
        m, p, j, k, sign = cp
        plt.axvline(x=m-j, linewidth=1, color='y')
        plt.axvline(x=m+k, linewidth=1, color='y')

    plt.show()


def test_blah():
    tests = [
        ('chr17', 40283388, 40287603, '+'),
        ('chr18', 31042592, 31068212, '-'),
        ('chr17', 40301912, 40314155, '+'),
        ('chr12', 53996031, 54000718, '-'),
        ('chr12', 55992235, 55997163, '+'),
        ('chr12', 56101529, 56104825, '+'),
        ('chr12', 56315821, 56317643, '-'),
        ('chr16', 19551194, 19555731, '+'),
        ('chr18', 31545721, 31550507, '+'),
        ('chr18', 31621550, 31625761, '-'),
        ('chr18', 31630959, 31631146, '-'),
        ('chr18', 31826338, 31830989, '-'),
        ('chr18', 32287031, 32288203, '-'),
        ('chr18', 32089621, 32092235, '+'),
    ]

    def tester(tup, ref):
        chrom, start, end, strand = tup
        rundir = '/Users/mkiyer/Documents/lab/projects/taco/ccle/ccle55_stringtie_changept'
        a = get_data(rundir, chrom, start, end, strand)
        print chrom, start, end
        ref_i = []
        for r in ref:
            ref_i.append(r - start)
        plot_slope(a, ref_i)

    # for test in tests:
    #     tester(test)

    # tester(('chr1',11936,14829,'-'), [14403])
    # tester(('chr1',183489,184971,'-'), [184922, 184924, 184926])
    # tester(('chr1',910143,917720,'-'), [916869])
    # tester(('chr1',9268168,9271563,'+'), [9271337])
    # tester(('chr1',781937,793041,'+'), [])
    # tester(('chr1',778769,787810,'+'), [])
    tester(('chr1',11936, 14829, '-'), [])


def test_smooth():
    return
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
