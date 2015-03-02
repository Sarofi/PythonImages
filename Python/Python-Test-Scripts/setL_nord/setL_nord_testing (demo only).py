import numpy as np
import math
import helper as h
import hetero as t
import sys
import scipy.linalg
import ty_nnz
import ty_any
import ty_isnumeric
import inspect
import time as tm
#from EMAN2 import EMNumPy
from scipy import io
import scipy.special as s
from sympy.mpmath import besselj, bessely, legendre
from structures import *
from numpy import linalg as LA
from math import pi, sqrt, sin, cos, acos, atan2, factorial
#from EMAN2 import EMData

import funcs
from funcs import xyztosph_vec, plgndr

def setL_nord(Rabc, ll, mm, vk, Htable, map_unique2lp, tilde_b):
    # copied from funcs for debugging and modification
    # omitted here for simplicity
    # returns a matrix L
    pass

def set_Ttable_nord(ll, nn, theta, phi,tilde_b):
    # copied from funcs for debugging and modification
    # omitted here for simplicity
    pass

for cnt in range(1, 1001):
    mat = io.loadmat('setL_nord_testing%d.mat' % cnt)
    
    # load inputs from matrices in .mat file
    Rabc = mat['Rabc']
    ll = mat['ll']
    mm = mat['mm']
    vk = mat['vk']
    Htable = mat['Htable']
    map_unique2lp = mat['map_unique2lp']
    tilde_b = mat['tilde_b']
    
    # get the results in its Python counterpart
    L = setL_nord(Rabc, ll, mm, vk, Htable, map_unique2lp, tilde_b)

    # load the original output in Matlab
    mat_L = mat['L']
    a, b = L.shape

    # threshold testing
    OUT_OF_BOUND = False
    for i in range(a):
        if OUT_OF_BOUND:
            break
        for j in range(b):
            if abs(L[i][j] - mat_L[i][j]) > 1e-7:
                print "error out of bound"
                OUT_OF_BOUND = True
                break
    if OUT_OF_BOUND:
        print "test %d failed" % cnt
