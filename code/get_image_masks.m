
function [AVF_mask,FVF_mask]=get_image_masks(sim_params)

% MAIN INPUTS
N = sim_params.N;            % number of axons i.e disks to pack  
d_mean = sim_params.d_mean;         % theoretical mean of axon diameters in um
d_var  = sim_params.d_var;         % theoretical variance of axon diameters in um
Delta  = sim_params.Delta;         % gap between the edge of axons in um 
iter_max = sim_params.iter_max;    % number of iteration i.e migrations to perform. Example: iter_max = 30000 ok if N = 1000

% SECONDARY INPUTS
threshold_high = sim_params.threshold_high;     % no diameter above 'threshold_high'
threshold_low = sim_params.threshold_low;     % no diameter under 'threshold_low'
iter_fvf = sim_params.iter_max/10;  % to study the packing convergence the disk density i.e Fiber Volume Fraction (FVF) can be computed and displayed every 'iter_fvf' iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     AxonPacking Process  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(d_mean)

    % axons features
    axons.N{k}      = N;
    axons.d_mean{k} = d_mean(k);
    axons.d_var{k}  = d_var(k);
    axons.Delta{k}           = Delta(k);
    axons.threshold_high{k} = threshold_high;
    axons.threshold_low{k}  = threshold_low;
    
    % axon diameters sampling (under a gamma law or lognormal and initialization of positions 'x0' in a square area of length 'side')
    [d, x0, side] = axons_setup(axons,'gamma', k);
    axons.d{k} = d;
    axons.g_ratio{k} = compute_gratio(d);
    
    % packing process of the axons
    [final_positions, final_overlap, fvf_historic] = process_packing(x0, d, Delta(k), side, iter_max, iter_fvf);
    
    % store packing results
    % main results
    packing.initial_positions{k}    = reshape(x0,2,length(x0)/2);
    packing.final_positions{k}      = final_positions;
    % secondary results
    packing.final_overlap{k}        = final_overlap;
    packing.FVF_historic{k}         = fvf_historic;
    packing.iter_max{k}             = iter_max;
    
    % Statistics from the packing
    [~,~,~,~,AVF_mask,FVF_mask] = compute_statistics(axons.d{k}, axons.Delta{k}, packing.final_positions{k}, side, axons.g_ratio{k});
    
end

end
