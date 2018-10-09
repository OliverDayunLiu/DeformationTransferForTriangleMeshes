import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import os
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

A = scipy.sparse.lil_matrix((10, 10), dtype=np.float128)
print A.toarray()

A[2,:] = 10
A[5,2] = 5
print A.toarray()