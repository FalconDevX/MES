import numpy as np

def clean_near_zero(value):
    if np.isclose(value, 0, atol=1e-12):
        value = 0.0
    return value