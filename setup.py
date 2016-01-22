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
__version__ = "1.0.1"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


extensions = [
    Extension('taco.lib.cBedGraph',
              sources=['taco/lib/cBedGraph.pyx'],
              include_dirs=[numpy_inc]),
    Extension('taco.lib.cChangePoint',
              sources=['taco/lib/cChangePoint.pyx'],
              include_dirs=[numpy_inc])
]


def main():
    setup(name='taco',
          version='0.4.0',
          description='transcriptome meta-assembly from rna-seq',
          author='Matthew Iyer and Yashar Niknafs',
          author_email='yniknafs@umich.edu',
          requires=['numpy', 'networkx', 'h5py'],
          license='GPL',
          platforms='Linux',
          url='https://github.com/yniknafs/taco',
          ext_modules=cythonize(extensions),
          packages=['taco'],
          scripts=['taco/taco_run.py'])

if __name__ == '__main__':
    main()
