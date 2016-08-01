
# Disks Close Random Packing - Application to the White Matter simulation

author : Tom Mingasson
contact : tom.mingasson@eleves.ec-nantes.fr
institution : University Polytechnique of Montreal, NeuroPoly
date : 2016

## Description 

Here is  a new random disks packing algorithm for numerical simulation of the white matter. 

White matter tissue was modelled with N parallel cylinders spaced from each other by a specific gap of various radii following a lognormal distribution
defined by σ and μ. Due to the geometric consistency along the axonal fibers, the problem was reduced to a 2-dimensional packing of perfectly round and 
non-compressible disks.

The software packing can provide microstructure features (axon density PHI, Myelin Volume Fraction MVF, Axon Volume Fraction AVF, fraction restricted FR) a
ssuming a constant 
g-ratio.

## Scripts

main.m 
axonsSetup.m 
processPacking.m
computeStatistics.m 
progressBar.m 

## How to use it ?

### INPUTS
In ‘main.m’ change the inputs

- nbAxons 	: the number of axons
- meanTheo 	: theoretical mean of axon diameters in um. Can be a vector if you want to create more than one packing. 
- varTheo 	: theoretical variances of axon diameters in um. Has to be a correspond with the size of ‘meanTheo’.
- threshold : no radius above this threshold when the log-normale radii distribution is generated.
- gapTheo 	: gap between axons in um.  Has to be a correspond with the size of ‘meanTheo’ and ‘varTheo’.
- g_ratio 	: g ratio assumed constant. 
- ITERmax 	: number of migrations to perform in the packing process.

#### Help 	
ITERmax = 5000 enough if N = 2000 axons 
ITERmax = 2000 enough if N = 1000 axons

#### Example  	
nbAxons = 1000;                                       
meanTheo = [3 3.5];         
varTheo = [1 1];                                   
threshold = 10;                                   
gapTheo = [0 0.3];                                 
g_ratio = 0.72;                                    
ITERmax = 2000;                             

### OUPUTS
The function ‘computeStatistics.m’ provides MVF, AVF, PHI and FR for each packing image defined by the input combinations. To do that it creates a binary masks
from the packing image. The user can set the resolution of this mask. The default resolution is 2048. 

Outputs are stored in data structures. Data structures outputs :
- in ‘axons.struct’ is stored : axon features (theoretical mean radii, theoretical variance radii, theoretical gap between axons, radii threshold, number of axons, set of the radii after sampling).
- in ‘packing.struct’ is stored : packing results (initial positions, final position, side and ITERmax). 
- in ‘statistics.struct’ is stored : resolution, MVF, AVF, PHI, FR. 

