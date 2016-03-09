'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy
numpy_inc = numpy.get_include()

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.1"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


cython_extensions = [
    Extension('taco.lib.cbedgraph',
              sources=['taco/lib/cbedgraph.pyx'],
              include_dirs=[numpy_inc]),
    Extension('taco.lib.cchangepoint',
              sources=['taco/lib/cchangepoint.pyx'],
              include_dirs=[numpy_inc]),
    Extension('taco.lib.bx.cluster',
              sources=['taco/lib/bx/cluster.pyx',
                       'taco/lib/bx/intervalcluster.c'],
              include_dirs=['taco/lib/bx']),
    Extension('taco.lib.bx.intersection',
              sources=['taco/lib/bx/intersection.pyx']),
    Extension('taco.lib.csuffixarray',
              sources=['taco/lib/csuffixarray.pyx', 'taco/lib/sais.c']),
    Extension('taco.lib.cpathfinder',
              sources=['taco/lib/cpathfinder.pyx']),
    Extension('taco.lib.cbisect',
              sources=['taco/lib/cbisect.pyx', 'taco/lib/bsearch.c']),
    Extension('taco.lib.scipy.norm_sf',
              sources=['taco/lib/scipy/norm_sf.pyx', 'taco/lib/scipy/ndtr.c'],
              include_dirs=[numpy_inc, 'taco/lib/scipy'])
]

extensions = [
    Extension('taco.lib.caggregate', sources = ['taco/lib/caggregate.c'])
]


def main():
    setup(name='taco',
          version=__version__,
          description='transcriptome meta-assembly from rna-seq',
          author='Matthew Iyer, Yashar Niknafs, Balaji Pandian',
          author_email='yniknafs@umich.edu',
          requires=['numpy', 'scipy', 'networkx', 'h5py', 'pysam'],
          license='GPL',
          platforms='Linux',
          url='https://github.com/yniknafs/taco',
          ext_modules=extensions + cythonize(cython_extensions),
          packages=['taco', 'taco.lib', 'taco.lib.bx'],
          scripts=['taco/taco_run.py'])


if __name__ == '__main__':
    main()
