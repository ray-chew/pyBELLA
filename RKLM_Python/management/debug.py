import numpy as np
import logging

# taken from: https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    logging.info("val = ", array.flat[idx])
    logging.info("idx = ", idx)
    return array.flat[idx], idx

