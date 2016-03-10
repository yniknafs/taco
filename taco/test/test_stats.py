'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import numpy as np

from scipy.stats import distributions

from taco.lib.scipy.norm_sf import norm_sf

from scipy.stats import mannwhitneyu as scipy_mwu
from taco.lib.stats import mannwhitneyu as mwu

def test_mannwhitneyu():
    x = [1, 2, 3, 4, 5]
    y = [6, 7, 8, 9, 10]
    p1 = scipy_mwu(x, y).pvalue
    p2 = mwu(x, y).pvalue / 2.0
    assert abs(p1 - p2) < 1e-5
