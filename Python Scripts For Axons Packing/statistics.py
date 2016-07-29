__author__ = 'Tom Mingasson'

import numpy as np
from math import *
import matplotlib.pyplot as plt



def computeStatistics(pts, radii, side, resolution=2048, g=0.72):

    # Create a mask of the axons
    xc = pts[:,0]
    yc = pts[:,1]
    xx, yy = np.mgrid[1:resolution+1, 1:resolution+1]

    mask = np.full((resolution,resolution), False, dtype=bool)

    for i in range(len(radii)):
        mask = np.logical_or(mask, np.sqrt((xx - xc[i] * (resolution/side))**2 + (yy - yc[i] * (resolution/side))**2) <= radii[i] * (resolution/side))

    Ls = np.sqrt(np.sum(pi*radii**2))
    Xmin = int((np.mean(pts[:,0]) - 2 * Ls/5) * resolution/side)
    Xmax = int((np.mean(pts[:,0]) + 2 * Ls/5) * resolution/side)
    Ymin = int((np.mean(pts[:,1]) - Ls/3) * resolution/side)
    Ymax = int((np.mean(pts[:,1]) + Ls/3) * resolution/side)

    maskTrunc = mask[Xmin:Xmax, Ymin:Ymax]
    maskTrunc = maskTrunc.astype(int)

    # Display Truncated Mask
    # plt.figure(1)
    # fig1 = plt.gcf()
    # plt.imshow(maskTrunc)
    # plt.show()


    # Compute micro structure characteristics

    Phi = float(np.sum(maskTrunc)) / ((Xmax - Xmin) * (Ymax - Ymin))    # = (fiber + myelin) / background

    Fr = g**2 * Phi / (1 + (g**2 - 1) * Phi)                            # = fiber / (fiber + background)

    AVF = g**2 * Phi                                                    # = fiber / area

    MVF = (1-g**2) * Phi                                                # = myelin / area

    return Phi, Fr, AVF, MVF