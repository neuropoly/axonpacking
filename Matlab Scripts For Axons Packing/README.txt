
############################################################################################################################################################

								Disks Close Random Packing
							Application to the White Matter simulation

############################################################################################################################################################

author : Tom Mingasson
contact : tom.mingasson@eleves.ec-nantes.fr
institution : University Polytechnique of Montreal, NeuroPoly
date : 2016


############################################################################################################################################################
									Description 
############################################################################################################################################################

Here is  a new random disks packing algorithm for numerical simulation of the white matter. 

White matter tissue was modelled with N parallel cylinders spaced from each other by a gap d of various radii following a lognormal distributiondefined by σ and μ. Due to the geometric consistency along the axonal fibers, the problem was reduced to a 2-dimensional packing of perfectly round and 
non-compressible disks.

The software packing that can easily provide microstructure features (axon density, Myelin Volume Fraction, Axon Volume Fraction, fraction restricted) 
assuming a constant g-ratio.

############################################################################################################################################################
									Scripts
############################################################################################################################################################
main.m 
axonsSetup.m 
processPacking.m
computeStatistics.m 
ProgressBar.m 
############################################################################################################################################################
									How to use it
############################################################################################################################################################

In ‘main.m’ change the inputs

- N 		: the number of axons
- meanTheo 	: theoretical mean of axon diameters in um. Can be a vector.
- varTheo 	: theoretical variances of axon diameters in um. Has to be a correspond with the size of ‘meanTheo’.
- threshold 	: no radius above this threshold when the log-normale radii distribution is generated.
- gapAxons 	: gap between axons in um.  Has to be a correspond with the size of ‘meanTheo’ and ‘varTheo’.
- g_ratio 	: g ratio assumed constant. 
- ITERmax 	: number of migrations to perform in the packing process.


HELP : 	ITERmax = 5000 enough if N = 2000 axons 
	ITERmax = 2000 enough if N = 1000 axons