# models.py
import numpy as np

class LinearModel:
    def __init__(self, c, d):
        self.c = float(c)
        self.d = np.asarray(d, dtype=float)

    def log_density(self, x):
        val = self.c + np.dot(self.d, x)
        if val <= 0.0:
            return -np.inf
        return np.log(val)

class ExponentialModel:
    def __init__(self, a, b):
        self.a = float(a)
        self.b = np.asarray(b, dtype=float)

    def log_density(self, x):
        return self.a + np.dot(self.b, x)
