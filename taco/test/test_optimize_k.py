'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from taco.lib.optimize import maximize_bisect

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def test_maximize_bisect():
    def f1(x):
        return x
    x, y = maximize_bisect(f1, 0, 100, 0)
    assert x == 100
    assert y == 100

    def f2(x):
        return -(x - 50) ** 2
    x, y = maximize_bisect(f2, 1, 100, 0)
    assert x == 50
    assert y == 0

    def f3(x):
        a = [1, 5, 11, 14, 16, 50, 100, 10000, 5, 4, 3, 2, 1]
        return a[x]
    x, y = maximize_bisect(f3, 1, 12, 0)
    print x, y
    assert x == 7
    assert y == 10000

    def f4(x):
        return 1
    x, y = maximize_bisect(f4, 1, 1, 0)
    assert x == 1
    assert y == 1
