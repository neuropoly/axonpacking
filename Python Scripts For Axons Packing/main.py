__author__ = 'Tom Mingasson'

import numpy as np
from math import *
import matplotlib.pyplot as plt
import time
from axonsFeatures import *

from setup import *
from processPacking import *
from statistics import *


start_time = time.time()
plt.close('all')

########################################################################################################################
# Axons Characteristics
########################################################################################################################

paramAxons = {'nbA': nbA, 'meanA': meanA, 'varA': varA, 'gapA': gapA, 'thresholdA': thresholdA}
print '\n#############################################################################################'
print '\nPACKING FEATURES'
print 'Number of Axons : ', nbA
print 'Mean diameter : ', meanA, 'um'
print 'Variance diameter : ', varA, 'um'
print 'No diameter above :', thresholdA, 'um'
print 'Gap between axons : ', gapA, 'um'
print 'g-ratio : ', g


########################################################################################################################
#                                                 Setup Packing
########################################################################################################################


# Parameters
# Area
sideCoeff = 4
side = sideCoeff * (meanA + gapA) * sqrt(nbA)

# Create Axons of mean diameter meanTheo and variance diameter variance diameter under a log-normale law distribution
radii = samplingLogNormale(paramAxons)

# Initialize the packing : random permutations on a grid
pts0 = initPacking(nbA, side)


########################################################################################################################
#                                                 Process Packing
# Itermax = 2000 is enough is nbAxons <= 1000
# IterMax = 4000 is enough is nbAxons <= 2000
########################################################################################################################
print  '\n#############################################################################################'

IterMax = 2000

print  '\nPacking in process... ', IterMax, 'iterations to perform'

pts = processPacking(pts0, radii, paramAxons, side, IterMax)

print '\nPacking has finished !'
########################################################################################################################
#                                                 Compute Characteristics
########################################################################################################################
resolution = 2048
Phi, Fr, MVF, AVF = computeStatistics(pts, radii, side, resolution, g)
print '\nRESULTS '
print 'Phi :', Phi
print 'Fr :', Fr
print 'MVF :', MVF
print 'AVF :', AVF

print("\nTotal execution time : %s seconds" % (time.time() - start_time))
########################################################################################################################
#                                                 Display Packing
########################################################################################################################

# Radii histograms
# plt.figure(1)
# fig1 = plt.gcf()
# bins = np.linspace(0, threshold + 2, 100)
# plt.hist(radii, bins)
# plt.show()


# Packing
plt.figure(2)
fig2 = plt.gcf()
for i in range(nbA):
    circle = plt.Circle(pts[i, :], radii[i])
    fig2.gca().add_artist(circle)
plt.axis([0, side, 0, side])
plt.show()