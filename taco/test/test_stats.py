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

    print scipy_mwu(x, y)
    print mwu(x, y)
