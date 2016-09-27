
# AxonPack - Application to the simulation of White Matter microstructure

author : Tom Mingasson    
contact : tom.mingasson@eleves.ec-nantes.fr          
institution : University Polytechnique of Montreal, NeuroPoly   
date : 2016 

<img src="https://github.com/neuropoly/axonpacking/blob/master/img1.jpeg" width="400px" align="middle" />

## Description 

Here is  a new random disks packing algorithm for numerical simulation of the white matter. 

White matter tissue is divided in three compartments:  axons, myelin sheath and extra-axonal space. Axons are assumed to be parallel cylinders, therefore the invariance along the fiber axis makes it possible to consider this problem in 2D. The dense packing of axons is thus equivalent to the generation of random 2-dimensional packing of N perfectly round and non-compressible disks. Axon diameter distributions follow a Gamma distribution (defined by its mean µ and variance σ2). Interestingly the g-ratio is fairly constant across species and white matter regions (31,32) and is dependent mostly on the diameter of the axon according to the relationship presented in (Ikeda M, Oka Y. Brain Behav. 2012):  gratio= 0.220 * log(DIAMETER_unmyelinated) +0.508. 

The different steps to process packing are the following: first, the diameters of the disks are randomly chosen using a gamma distribution parameterized with the mean (µ), standard deviation (σ) and number of axons (N).  Then, the positions of disks are initialized on a grid, and they migrate toward the center of the packing area until the maximum disk density is achieved. 


The software packing provides microstructure features (Fiber Volume Fraction FVF, Myelin Volume Fraction MVF, Axon Volume Fraction AVF, fraction restricted FR).

<img src="https://github.com/neuropoly/axonpacking/blob/master/img2.jpeg" width="1000px" align="middle" />

## Scripts

main.m
axonsSetup.m
processPacking.m
computeStatistics.m
progressBar.m

## How to use it ?

### INPUTS
In ‘main.m’ change the inputs

- numberAxons        % number of axons (N)
- mean_theoretical   % theoretical mean of axon diameters in um. Can be a vector if you want to create more than one packing.
- var_theoretical    % theoretical variances of axon diameters in um. Has to be a correspond with the size of ‘mean_theoretical’.
- gap_theoretical    % gap between axons in um 
- threshold_high     % no diameter above 'threshold_high'
- threshold_low      % no diameter under 'threshold_low'
- ITERmax            % number of iteration i.e migrations to perform. Example: ITERmax = 30000 ok if N = 1000
- ITERfvf            % the disk density i.e Fiber Volume Fraction (FVF) is computed and displayed every 'ITERfvf' iterations

#### Help 	

The disk density increases over the migrations and tends toward a limit value. It is necessary to first launch the algorithm with the packing inputs (N, µ, σ and Δ) and a high number of iterations. MRI metrics such as the disks density e.g. FVF can be calculated every p iterations to assess the sufficient number of iterations to reach a certain degree of precision. p is a user defined integer: p = 250 or 1000 for example. 

When mean_theoretical (μ) closed to 3, var_theoretical (σ2) closed to 1 and numberAxons about 1000, ITERmax = 30000 is sufficient. 


#### Example  	
numberAxons = 25;
mean_theoretical = [3 3]; 
var_theoretical  = [1 1];                    
gap_theoretical  = [0 0.5];                                
threshold_high = 10;                                     
threshold_low = 0.2;                                         
ITERmax = 30000;                       
ITERfvf = 1000;                             

### OUPUTS
The function ‘computeStatistics.m’ provides MVF, AVF, FVF and FR for each packing image defined by the input combinations. To do that it creates a binary masks
from the packing image. The user can set the resolution of this mask. The default resolution is 2048. 

Outputs are stored in data structures. Data structures outputs :
- in ‘axons.mat’ is stored : axon features (theoretical mean diameters, theoretical variance diameters, theoretical gap between axons, diameters threshold low and high, number of axons, set of the diameters after sampling).
- in ‘packing.mat’ is stored : packing results (initial positions, final position, final overlap ratio in the paking, length of the packing area 'side' and ITERmax). 
- in ‘stats.mat’ is stored : resolution, MVF, AVF, FVF, FR, g-ratio. 

