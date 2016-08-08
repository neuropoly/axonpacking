% author : Tom Mingasson

close all
clear variables
clc

addpath('... path for other scripts ...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE BELOW : AXONS FEATURES and ITERmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbAxons = 1000;                                       % number of axons
meanTheo = [3 3.5]; Lm = length(meanTheo);          % theoretical mean of axon diameters in um
varTheo = [1 1];                                    % theoretical variances of axon diameters in um
threshold = 10;                                     % no radii above 'threshold'
gapTheo = [0 0.3];                                  % gap between axons in um 
g_ratio = 0.72;                                     % g ratio
ITERmax = 2500;                                     % ITERmax = 2500 ok if N = 1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k=1:Lm
    %% Process Packing
    % axons features
    axons.meanTheo{k}  = meanTheo(k);
    axons.varTheo{k}   = varTheo(k);
    axons.gapTheo{k}   = gapTheo;
    axons.threshold{k} = threshold;           % no diameter above 'threshold' (in micro metres)
    axons.nbAxons{k}   = nbAxons;
    axons.gRatio{k}    = g_ratio;

    % dimension of the square area
    sideCoeff = 4;
    side = sideCoeff*(meanTheo(k)+gapTheo(k))*sqrt(nbAxons);
    
    % axons diameters sampling and initialization of positions x0
    [R,x0] = axonsSetup(axons,side,k);
    
    % packing of the axons
    x = processPacking(x0,R,gapTheo(k),side,ITERmax);
    
    % results
    axons.radii{k}                    = R;
    packing.initialPositions{k}       = x0;
    packing.finalPositions{k}         = reshape(x,2,length(x)/2);
    packing.side{k}                   = side;
    packing.ITERmax{k}                = ITERmax;
    
    %% Statistics from the packing
    resolution = 2048;
    [PHI, FR, MVF, AVF] = computeStatistics( R, x, side, g_ratio, resolution);
    
    % results
    statistics.resolution{k} = resolution;
    statistics.PHI{k} = PHI;
    statistics.FR{k}  = FR;
    statistics.MVF{k} = MVF;
    statistics.AVF{k} = AVF;
    
    %% Display
    % display the packing with results annotations for Phi and Fr
    set( figure ,'WindowStyle','docked' ); clf
    pts = reshape(x,2,length(x)/2);
    t = 0:.1:2*pi+0.1;
    for i=1:nbAxons
        title(['Diameter Mean : ',num2str(mean(R(:))),' µm    ','Diameter Variance : ',num2str(var(R(:))),' µm    ','Gap : ',num2str(gapTheo(k)),' µm    '],'FontSize',10,'FontWeight','bold');
        patch(R(i)*cos(t) + pts(1,i), R(i)*sin(t) + pts(2,i), [.5 .5 1], 'FaceAlpha', 0.4, 'EdgeColor', [.2 .2 1]);
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
% save_var  = num2str(varTheo); save_var(save_var == ' ') = '';
% save_mean = num2str(meanTheo); save_mean(save_mean == ' ') = '';
% save_gap  = num2str(gapTheo); save_gap(save_gap == ' ') = '';
% saveName  = ['Axons', num2str(N), '_Mean', save_mean, '_Var', save_var, '_Gap', save_gap, '_Sim', num2str(NbSim)]
% mkdir('... path for saving ...', saveName)
% 
% cd(['... path ...',saveName])
% save('axons.mat', '-struct', 'axons');
% save('packing.mat', '-struct', 'packing');
% save('optimization.mat', '-struct', 'optimization');
% save('statistics.mat', '-struct', 'statistics');






