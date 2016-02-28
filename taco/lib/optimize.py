'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def maximize_bisect(f, xmin, xmax, max_iterations=0):
    if max_iterations == 0:
        max_iterations = xmax

    x = (xmin + xmax) / 2
    y = {}
    iterations = 0
    while iterations < max_iterations:
        if (xmax - xmin) <= 3:
            best_y = None
            best_x = None
            for x1 in xrange(xmin, xmax + 1):
                y1 = y[x1] if x1 in y else f(x1)
                y[x1] = y1
                if (best_y is None) or (y1 > best_y):
                    best_y = y1
                    best_x = x1
            return best_x, best_y

        xlist = ((x + xmin) / 2, x, (x + xmax) / 2)
        best_y = None
        best_i = None
        for i in xrange(3):
            x1 = xlist[i]
            y1 = y[x1] if x1 in y else f(x1)
            y[x1] = y1
            if (best_y is None) or (y1 > best_y):
                best_y = y1
                best_i = i
        if best_i == 0:
            xmax = x
            x = xlist[0]
        elif best_i == 1:
            xmin = xlist[0]
            xmax = xlist[2]
        else:
            xmin = x
            x = xlist[2]
    return None, 0
