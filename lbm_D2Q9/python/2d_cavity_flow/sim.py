from numba import int64, float64
from numba.experimental import jitclass

# Define the data types for each attribute
spec = [
    ('x', int64),
    ('y', int64),
    ('rho', float64)
]


@jitclass(spec)
class node:
    def __init__(self, x, y, rho):
        self.x = x
        self.y = y
        self.rho = rho
