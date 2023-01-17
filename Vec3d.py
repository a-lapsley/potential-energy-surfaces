import numpy as np

class Vec3d(np.ndarray):

    def __new__(cls, lst = (0, 0, 0)):
        ar = np.array(lst, dtype=float)
        x = np.ndarray.__new__(cls, shape=(3,), dtype=float, buffer=ar)
        return x
    
    def length(self):
        return np.sqrt(np.dot(self, self))