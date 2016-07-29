__author__ = 'Tom Mingasson'

import numpy as np
from math import *

def samplingLogNormale(paramAxons):

    N = paramAxons["nbA"]

    sigma_instru = 1.0;
    x = np.zeros(N)
    x[0] = 1.0;

    for n in range(N-1):

        x_star = np.random.normal(x[n], sigma_instru, 1)
        proba_accept = (funInstru(x[n], x_star, sigma_instru) * funObj(x_star, paramAxons)) / (funInstru(x_star, x[n], sigma_instru) * funObj(x[n], paramAxons))
        u = np.random.uniform(0, 1, 1)

        if u < proba_accept:
            x[n + 1] = x_star
        else:
            x[n + 1] = x[n]
    return x

def funObj(x, paramAxons):
    # Probability pobj(x) where funObj is the function we want to sample
    # Here funObj is a log-normale function

    mean = float(paramAxons["meanA"])
    var = float(paramAxons["varA"])
    threshold = float(paramAxons["thresholdA"])

    # radii distribution
    mu = log(mean) - (1./ 2) * log(1 + var / (mean)**2)
    sigma = sqrt(log(1 + var / (mean)**2))

    # no diameter above threshold
    if x >= threshold or x <= 0:
        f_logn = 0
    else:
        f_logn = 1 / (sigma * x * sqrt(2 * pi)) * exp(-(log(x) - mu) ** 2 / (2 * sigma**2))
    return f_logn

def funInstru(x_star, x_current, sigma_instru):
    q = 1.0 / (sigma_instru * sqrt(2 * pi)) * exp(-(x_star - x_current)**2 / (2 * sigma_instru**2))
    return q

def initPacking(nbAxons, side):
    pts0 = np.zeros((nbAxons,2))
    sqrtTmp = int(sqrt(nbAxons))+1
    Xgrid, Ygrid = np.meshgrid(np.arange(0, side, side/sqrtTmp), np.arange(0, side, side/sqrtTmp))
    Xgrid = np.reshape(Xgrid, sqrtTmp**2); Ygrid = np.reshape(Ygrid, sqrtTmp**2);
    P = np.random.permutation(sqrtTmp**2)
    for i in range(nbAxons):
        pts0[i,:] = [Xgrid[P[i]], Ygrid[P[i]]]
    return pts0