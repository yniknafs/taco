import numpy as np

# define array types
FLOAT_DTYPE = np.float32
H5_CHUNKSIZE = (1 << 18) / np.dtype(FLOAT_DTYPE).itemsize
