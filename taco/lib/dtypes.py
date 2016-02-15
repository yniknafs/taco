import numpy as np

# define array types
FLOAT_DTYPE = np.float32
INT_DTYPE = np.int_

H5_CHUNKSIZE = (1 << 18) / (3 * np.dtype(FLOAT_DTYPE).itemsize)
