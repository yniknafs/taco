'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import pytest
import numpy as np
from scipy.stats import mannwhitneyu, percentileofscore
from taco.lib.base import TacoError, Strand
from taco.lib.dtypes import FLOAT_DTYPE
from taco.lib.cChangePoint import mse as mse_cython
from taco.lib.changepoint import mse as mse_python
from taco.lib.changepoint import run_changepoint, smoother
import os

import matplotlib.pyplot as plt

def get_data(rundir, chrom, start, end, strand):
    import h5py
    filename = os.path.join(rundir, 'expression.h5')
    f = h5py.File(filename, 'r')
    if strand == '+': strand_idx = Strand.POS
    elif strand == '-': strand_idx = Strand.NEG
    else: strand_idx = Strand.NA
    a = f[chrom][strand_idx, start:end]
    return a

# def test_mse():
#     a = np.zeros(100, dtype=FLOAT_DTYPE)
#     a[:50] = 10
#     a[50:] = 0
#     a += np.random.random(100)
#     mse_min_c, mse_i_c = mse_cython(a)
#     mse_min_py, mse_i_py = mse_python(a)
#     assert mse_i_c == 50
#     assert mse_i_c == mse_i_py
#     assert np.allclose(mse_min_c, mse_min_py, atol=0.01)


# def test_smooth():
#     a = np.zeros(1000)
#     a[:200] = 100
#     a[200:300] = np.linspace(100, 0, num=100)
#     a[20] = 8
#     a[40] = 8
#     a[600:700] = np.linspace(0, 500, num=100)
#     a[700:] = 500
#     a += 50 * np.random.random(1000)


#     for wsize in xrange(11, 101):
#         y = smooth(a, window_len=52, window='hanning')
#         # print len(a), len(y), wsize

#     print 'a', len(a)
#     print 'y', len(y)

#     # plt.figure(1)
#     # plt.subplot(411)
#     # plt.plot(a)
#     # plt.subplot(412)
#     # plt.plot(y)
#     # plt.subplot(413)
#     # plt.plot(np.gradient(y))
#     # plt.subplot(414)
#     # y = np.ediff1d(a)
#     # y = np.gradient(a)

#     # print 'ediff1a', len(y)

#     y = smooth(y, window_len=52, window='flat')
#     # plt.plot(y)

#     # plt.show()



# def get_data(chrom, start, stop, strand):
#     import h5py
#     filename = '/Users/yashar/Documents/taco/test/expression.h5'
#     f = h5py.File(filename, 'r')
#     print strand, 'po'
#     if strand == '+': strand_idx = Strand.POS
#     elif strand == '-': strand_idx = Strand.NEG
#     else: strand_idx = Strand.NA
#     a = f[chrom][strand_idx, start:stop]
#     return a


# def permute(a, nperms=10):
#     d = np.ediff1d(a)
#     positions = np.arange(a.shape[0])
#     start_values = d[d > 0]
#     stop_values = d[d < 0]

#     e = np.zeros((nperms, a.shape[0]))
#     for i in xrange(nperms):
#         print i
#         e[i, :] += a[0]
#         start_pos = np.random.choice(positions, size=len(start_values), replace=False)
#         for j in xrange(len(start_values)):
#             e[i, :start_pos[j]] += start_values[j]
#         stop_pos = np.random.choice(positions, size=len(stop_values), replace=False)
#         for j in xrange(len(stop_values)):
#             e[i, stop_pos[j]:] += stop_values[j]

#     return e




# # def plot_binseg(a, s_a, ref_pos=None):
# #     import matplotlib.pyplot as plt

# #     #a = np.ediff1d(a)
# #     S, Sdiff, Sm_i = cusum(a)
# #     mse_val, m = mse(a)
# #     md_val, md_m = md(a)
# #     # cps = bin_seg(mse, a, s_a)
# #     cps = bin_seg_slope(mse, a, s_a)
# #     plt.figure(1)
# #     plt.subplot(211)
# #     plt.plot(a, 'k')
# #     plt.subplot(212)
# #     plt.plot(S, 'r--')

# #     for i in xrange(10):
# #         b = a.copy()
# #         np.random.shuffle(b)
# #         Sb, Sbdiff, j = cusum(b)
# #         plt.plot(Sb)

# #     for i,p in cps:
# #         s = np.sign(s_a[i])
# #         if s == 1: 
# #             plt.axvline(x=i, linewidth=1, color='r')
# #         if s == -1: 
# #             plt.axvline(x=i, linewidth=1, color='b')

# #     if ref_pos is not None:
# #         plt.axvline(x=ref_pos, linewidth=1, color='g', linestyle='--')

# #     plt.show()


# def slope_extract(slope_a, i):
#     slope_a = np.sign(slope_a)
#     if slope_a[i] == 0: 
#         j=0
#         k=0
#     else:
#         j=0
#         while slope_a[i-j]==slope_a[i]: 
#             j+=1
            
#         k=0
#         while slope_a[i+k]==slope_a[i] :
#             k+=1
#             if (i+k)==len(slope_a): break
#     return (j, k) 

# def smoother(a, smooth_window="hanning", window_size = 75):
#     return smooth(a, window_len=window_size, window=smooth_window)

def plot_slope(a, ref_pos=None):
    import matplotlib.pyplot as plt

    print 'wangchung'
    smooth_a = np.array(smoother(a), dtype=FLOAT_DTYPE)
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

# def mwu(a, i):
#     a1 = a[:i]
#     a2 = a[i:]
#     if (len(np.unique(a1)) == 1) and np.array_equal(np.unique(a1), np.unique(a2)):
#         return (999999999, 1)
#     else:
#         U, p = mannwhitneyu(a1, a2)
#         # p = p * len(a)
#         return (U, p)

# def mwu_ediff(a, i):
#     a1 = a[:i]
#     a2 = a[i:]
#     a1 = a1[np.ediff1d(a1).nonzero()[0]]
#     a2 = a2[np.ediff1d(a2).nonzero()[0]]
#     if len(a1) == 0 or len(a2) == 0: 
#         return (999999999, 1)
#     elif (len(np.unique(a1)) == 1) and np.array_equal(np.unique(a1), np.unique(a2)):
#         return (999999999, 1)
#     else:
#         U, p = mannwhitneyu(a1, a2)
#         # p = p * len(a)
#         return (U, p)


# def bin_seg_slope(cp_func, a, s_a, PVAL=.05, cps=None, offset=0, size_cutoff=20):
#     print 'boobies', offset, offset+len(a)
#     if a.shape[0] < size_cutoff: return cps
#     if cps is None: 
#         cps = []
#     stat, i = cp_func(a)
#     U,p = mwu_ediff(a, i)
#     print 'p:', p, 'i:', i
#     if i < 2: return cps
#     elif p < PVAL:
#         j, k = slope_extract(s_a, offset+i)
#         if j!=0 and k!=0:
#             stat1, i1 = cp_func(a[:i])
#             stat2, i2 = cp_func(a[i:])
#             print "cost", stat, stat1, stat2 
#             # print i, offset, offset+i, p, offset+i-j, offset+i+k 
#             cps.append((i + offset, p, j, k))
#             #test left
#             if (offset+i-j) > offset: 
#                 b1 = a[:i-j]
#                 cps = bin_seg_slope(cp_func, b1, s_a, cps=cps, offset=offset)
#             #test right
#             if (offset+i+k) < offset+len(a): 
#                 b2 = a[i+k:]
#                 cps = bin_seg_slope(cp_func, b2, s_a, cps=cps, offset=offset + i+k)

#     return cps

# def bin_seg(cp_func, a, s_a, PVAL=1e-5, cps=None, offset=0):
#     if cps is None: 
#         cps = []
#     stat, i = cp_func(a)
#     U,p = mwu_ediff(a, i)
#     if i < 2: return cps
#     elif p < PVAL:
#         print len(a), i, p
#         cps.append((i + offset, p))
#         #test left
#         b1 = a[:i]
#         cps = bin_seg(cp_func, b1, s_a, cps=cps, offset=offset)
#         #test right
#         b2 = a[i:]
#         cps = bin_seg(cp_func, b2, s_a, cps=cps, offset=i)

#     return cps

def test_blah():
    import h5py
    # slope_size = 50

    # a = np.zeros(1000)
    # a[:200] = 200
    # a[200:300] = np.linspace(200, 0, num=100)
    # a[20] = 8
    # a[40] = 8
    # a[600:700] = np.linspace(0, 500, num=100)
    # a[700:] = 500
    # a += 50 * np.random.random(1000)

    # a = np.array(a, dtype=FLOAT_DTYPE)
    # s_a = np.gradient(smoother(a))
    # s_a = np.array(s_a, dtype=FLOAT_DTYPE)
    # # plot_binseg(a, s_a)
    # # plot_slope(a, s_a)


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
        rundir = '/Users/yashar/Documents/taco/test'
        a = get_data(rundir, chrom, start, end, strand)
        s_a = np.gradient(smoother(a))
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
    tester(('chr1',781937,793041,'+'), [])

