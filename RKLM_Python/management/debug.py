import numpy as np

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    print("val = ", array.flat[idx])
    print("idx = ", idx)
    return array.flat[idx], idx

