__author__ = 'Tom Mingasson'

import numpy as np
from scipy.spatial import distance
from axonsFeatures import *


def processPacking(pts0, radii, paramAxons, side, IterMax):
    pts = pts0
    for i in range(IterMax):
        MyGrad = ComputeGrad(pts, radii + gapA * np.ones(radii.shape), side)
        pts = pts + MyGrad
    return pts


def ComputeGrad(pts, R, side):

    # Parameters
    Kcenter0 = 0.01     # center step coefficient for disks withOUT overlapping
    Kcenter1 = 0        # center step coefficient for disks with overlapping
    KrepMULT = 1.5
    Krep = KrepMULT * Kcenter0  # repulsion step coefficient for disks with overlapping
    PowerRep = 0.5              # repulsion power for disks with overlapping

    x = pts.ravel()

    # Attract axons to center area
    ptsCenter = (side / 2) * np.ones(pts.shape) - pts
    dist2center = np.sqrt(ptsCenter[:, 0] ** 2 + ptsCenter[:, 1] ** 2)
    dist2center[dist2center == 0] = 1
    gradCenter = np.divide(ptsCenter, np.tile(dist2center, (2, 1)).T)

    # Repulsion between axons which are overlapping
    P = distance.squareform(distance.pdist(pts, 'euclidean')) + np.identity(nbA)
    Rsum = (np.tile(R.T, (nbA, 1)) + np.tile(R, (nbA, 1)).T) * (np.tril(np.ones((nbA,nbA)), -1) + np.tril(np.ones((nbA,nbA)), -1).T)
    L = np.divide(Rsum,P) - np.ones((nbA,nbA)); L[L == np.inf]=10
    Linter = L * (L>0)
    Linter2 = np.reshape(np.repeat(Linter,2),(nbA, 2*nbA))
    LinterSum = np.sum(np.abs(Linter2),0); LinterSum[LinterSum == 0] = 1; LinterSum = np.tile(LinterSum,(nbA, 1));
    LinterNorm = np.divide(Linter2,LinterSum)
    LinterWeight = np.abs(LinterNorm)**PowerRep * np.sign(LinterNorm)
    U = np.tile(x,(nbA,1)) - np.tile(pts,(1,nbA))
    gradRepulsion = U*LinterWeight

    # Final Gradient
    LinterBIN1 = np.tile(np.sum(np.triu(L > 0), 1), (2, 1))
    LinterBIN0 = np.ones((2,nbA)) - LinterBIN1
    MyGrad = (Kcenter0*gradCenter*LinterBIN0.T + Kcenter1*gradCenter*LinterBIN1.T).ravel() + Krep * np.sum(gradRepulsion,0).T
    MyGrad = np.reshape(MyGrad, (nbA, 2))

    return  MyGrad



