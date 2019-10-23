import numpy as np

# taken from: https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    print("val = ", array.flat[idx])
    print("idx = ", idx)
    return array.flat[idx], idx

