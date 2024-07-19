# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 13:44:02 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

import numpy as np
from scipy.linalg import solve

keff = 1.11706

mtx = np.array([[0.89357, 0.28779, 0.15006, 0.05032],
                [0.89357, 0.28779, 0.15006, 0.05032],
                [0.12139, 0.28695, 0.38502, 0.19824],
                [0.16864, 0.43211, 0.58580, 0.75710]])

src = np.array([1.0, 1.0, 1.0, 1.0])


detM = np.linalg.det(mtx)

src0 = np.dot(mtx, src)

