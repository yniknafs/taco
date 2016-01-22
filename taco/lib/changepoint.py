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


def find_change_points(a, start=0, threshold=0):
    '''
    a - numpy array of expression data (all values >= 0)
    start - genomic start of 'a'
    threshold - expression threshold to detect changes
    '''
    if a.shape[0] == 0:
        return []
    if a.shape[0] == 1:
        return []
    change_points = []
    for i in xrange(1, a.shape[0]):
        if (a[i-1] > threshold) and (a[i] <= threshold):
            change_points.append(start + i)
        elif (a[i-1] <= threshold) and (a[i] > threshold):
            change_points.append(start + i)
    return change_points
