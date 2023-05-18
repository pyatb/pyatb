import numpy as np
from typing import Union

SQRT_PI = np.sqrt(np.pi)

def gauss(sigma: float, x: Union[np.ndarray, float]) -> Union[np.ndarray, float]:
    temp_x = x * x / sigma**2
    ans = 1.0 / (sigma * SQRT_PI) * np.exp(-temp_x)
    return ans