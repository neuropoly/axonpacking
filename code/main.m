% author : Tom Mingasson

close all
clear variables
clc

% addpath('... path for other scripts ...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE BELOW : AXONS FEATURES and ITERmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberAxons = 25;                                           % number of axons
mean_theoretical = [3 3]; Lm = length(mean_theoretical);     % theoretical mean of axon radii in um
var_theoretical  = [1 1];                                    % theoretical variances of axon radii in um
gap_theoretical  = [0 0.5];                                  % gap between axons in um 
threshold_high = 10;                                         % no diameter above 'threshold_high'
threshold_low = 0.2;                                         % no diameter under 'threshold_low'
ITERmax = 30000;                                             % number of iteration i.e migrations to perform. Example: ITERmax = 30000 ok if N = 1000
ITERfvf = 1000;                                              % the disk density i.e Fiber Volume Fraction (FVF) is computed and displayed every 'ITERfvf' iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k=1:Lm
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
    statistics.resolution{k} = resolution;
    statistics.FVF{k}        = FVF;
    statistics.FR{k}         = FR;
    statistics.MVF{k}        = MVF;
    statistics.AVF{k}        = AVF;
    statistics.g_ratio{k}    = g_ratio;

    %% Display
    % display the packing with results annotations for FVF and Fr
    set( figure ,'WindowStyle','docked' ); clf
    t = 0:.1:2*pi+0.1;
    for i=1:numberAxons
        title(['Diameter Mean : ',num2str(mean(D(:))),' µm    ','Diameter Variance : ',num2str(var(D(:))),' µm    ','Gap : ',num2str(gap_theoretical(k)),' µm    '],'FontSize',10,'FontWeight','bold');
        patch(D(i)*cos(t) + finalPositions(1,i), D(i)*sin(t) + finalPositions(2,i), [.5 .5 1], 'FaceAlpha', 0.4, 'EdgeColor', [.2 .2 1]);
        xlim([0 side])
        ylim([0 side])
    end
    strPHI = [ 'Phi = ', num2str(roundn(PHI,-4))];
    strFR =  [ 'Fr = ',  num2str(roundn(FR,-4))];
    str = {strPHI, strFR};
    annotation('textbox','String',str,'FontSize',10, 'HorizontalAlignment','center','LineStyle','none','BackgroundColor',[1 1 1]);
    drawnow
    
end


%% save results
% save_var  = num2str(var_theoretical);  save_var(save_var == ' ') = '';
% save_mean = num2str(mean_theoretical); save_mean(save_mean == ' ') = '';
% save_gap  = num2str(gap_theoretical);  save_gap(save_gap == ' ') = '';
% save_iter  = num2str(ITERmax);
% saveName  = ['Axons', num2str(N), '_Mean', save_mean, '_Var', save_var, '_Gap', save_gap, '_Iter',save_iter]
% mkdir('... path for saving ...', saveName)
% 
% cd(['... path ...',saveName])
% save('axons.mat', '-struct', 'axons');
% save('packing.mat', '-struct', 'packing');
% save('optimization.mat', '-struct', 'optimization');
% save('statistics.mat', '-struct', 'statistics');






