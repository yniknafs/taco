'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import numpy as np

from taco.lib.base import Exon, Strand
from taco.lib.bedgraph import array_to_bedgraph

def test_omit_zeros():
    a = np.array([1,0,1,0])
    bg = tuple(array_to_bedgraph(a, omit_zeros=False))
    assert bg == ((0, 1, 1), (1, 2, 0), (2, 3, 1), (3, 4, 0))
    bg = tuple(array_to_bedgraph(a))
    assert bg == ((0, 1, 1), (2, 3, 1))

def test_single_run():
    a = np.ones(5) * 5
    bg = tuple(array_to_bedgraph(a))
    assert bg == ((0, 5, 5),)

def test_zero_sized_array():
    a = np.zeros(0)
    bg = tuple(array_to_bedgraph(a))
    assert bg == tuple()

def test_array_with runs():
    a = np.array([1,2,2,3,3,3,4,4,4,4,5,5,5,5,5])
    bg = tuple(array_to_bedgraph(a))
    correct = ((0, 1, 1), (1, 3, 2), (3, 6, 3), (6, 10, 4), (10, 15, 5))
    assert bg == correct
    a = np.array([1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,0,0,0,0])
    bg = tuple(array_to_bedgraph(a))
    assert bg == correct
