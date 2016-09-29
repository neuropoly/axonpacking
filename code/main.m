% author : Tom Mingasson

close all
clear variables
clc

% addpath('... path to the code folder ...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE BELOW : AXONS FEATURES and ITERmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberAxons = 100;                                           % number of axons 
mean_theoretical = 3;                                        % theoretical mean of axon radii in um
var_theoretical  = 1;                                        % theoretical variances of axon radii in um
gap_theoretical  = 0;                                        % gap between axons in um 
threshold_high = 10;                                         % no diameter above 'threshold_high'
threshold_low = 0.2;                                         % no diameter under 'threshold_low'
ITERmax = 20000;                                             % number of iteration i.e migrations to perform. Example: ITERmax = 30000 ok if numberAxons = 1000
ITERfvf = 1000;                                              % the disk density i.e Fiber Volume Fraction (FVF) is computed and displayed every 'ITERfvf' iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(mean_theoretical)
    %% Process Packing
    % axons features
    axons.mean_theoretical{k} = mean_theoretical(k);
    axons.var_theoretical{k}  = var_theoretical(k);
    axons.gap_theoretical{k}  = gap_theoretical(k);
    axons.threshold_high{k}   = threshold_high;
    axons.threshold_low{k}    = threshold_low;
    axons.numberAxons{k}      = numberAxons;
    
    % axons diameters sampling (under a gamma law or lognormal and initialization of positions 'x0' in a square area of length 'side'
    [D, x0, side] = axonsSetup(axons,'gamma',k);
    
    % packing of the axons
    [finalPositions, overlap, FVF_historic] = processPacking(x0, D, gap_theoretical(k), side, ITERmax, ITERfvf);
    
    % results
    axons.diameters{k}                = D;
    packing.initialPositions{k}       = reshape(x0,2,length(x0)/2);
    packing.side{k}                   = side;
    packing.finalPositions{k}         = finalPositions;
    packing.finalOverlap{k}           = overlap;
    packing.FVF_historic{k}           = FVF_historic;
    packing.ITERmax{k}                = ITERmax;
    
    %% Statistics from the packing
    resolution = 2048;
    [FVF, FR, MVF, AVF, g_ratio] = computeStatistics(axons.diameters{k}, axons.gap_theoretical{k}, packing.finalPositions{k}, packing.side{k}, resolution);
    
    % results
    stats.resolution{k} = resolution;
    stats.FVF{k}        = FVF;
    stats.FR{k}         = FR;
    stats.MVF{k}        = MVF;
    stats.AVF{k}        = AVF;
    stats.g_ratio{k}    = g_ratio;

    
end


%% save results
% save_var  = num2str(var_theoretical);  save_var(save_var == ' ') = '';
% save_mean = num2str(mean_theoretical); save_mean(save_mean == ' ') = '';
% save_gap  = num2str(gap_theoretical);  save_gap(save_gap == ' ') = '';
% save_iter  = num2str(ITERmax);
% saveName  = ['Axons', num2str(numberAxons), '_Mean', save_mean, '_Var', save_var, '_Gap', save_gap, '_Iter',save_iter]
% mkdir('... path for saving ...', saveName)
% 
% cd(['... path ...',saveName])
% save('axons.mat', '-struct', 'axons');
% save('packing.mat', '-struct', 'packing');
% save('stats.mat', '-struct', 'stats');






