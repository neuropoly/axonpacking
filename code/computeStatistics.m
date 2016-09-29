function [FVF, FR, MVF, AVF, g] = computeStatistics( D, gap, pts, side, resolution)
% stat evaluate from the results of the simulation :
%       - FVF : the fiber volume fraction e.g. the axon (=disk) density
%       - MVF : the myelin volume fraction
%       - AVF : the axon volume fraction
%       - FR  : the fraction of restricted water

N = length(D);
g =  compute_gratio(D);
t = 0:.1:2*pi+0.1;

% FVF mask
masksize = ceil(resolution*sqrt(N)/sqrt(500));
FVF_mask = false(masksize);
for id=1:N
    Xfibers = D(id)*cos(t) + pts(1,id);
    Yfibers = D(id)*sin(t) + pts(2,id);
    FVF_mask = FVF_mask | poly2mask(Xfibers/side*masksize, Yfibers/side*masksize, masksize, masksize);
end

% AVF mask
AVF_mask = false(masksize);
g_ratio=compute_gratio(D);
for id=1:N
    Xaxons = g_ratio(id)*D(id)*cos(t) + pts(1,id);
    Yaxons = g_ratio(id)*D(id)*sin(t) + pts(2,id);
    AVF_mask = AVF_mask | poly2mask(Xaxons/side*masksize, Yaxons/side*masksize, masksize, masksize);
end

% masks trunc
Ls = sqrt(sum(pi*(D+gap/2).^2))*(4/5)/side*masksize;
Xmin = round(mean(pts(1,:))/side*masksize - Ls/2);
Xmax = round(mean(pts(1,:))/side*masksize + Ls/2);
Ymin = round(mean(pts(2,:))/side*masksize - Ls/2);
Ymax = round(mean(pts(2,:))/side*masksize + Ls/2);

FVF_mask_trunc = FVF_mask(Xmin:Xmax,Ymin:Ymax);
AVF_mask_trunc = AVF_mask(Xmin:Xmax,Ymin:Ymax);

area1 = sum(FVF_mask_trunc(:));
area0 = Ls*Ls - area1;
area2 = sum(AVF_mask_trunc(:));

% FVF
FVF = area1 / (Ls*Ls); % = (AXON + myelin) / area

% FR
FR = area2 / (area2 + area0); % = AXON / (AXON + background)

% AVF
AVF = area2 / (Ls*Ls); % = AXON / area

% MVF
MVF = (area1 - area2) / (Ls*Ls); % = myelin / area

end

function g = compute_gratio(R)

% Ikeda M, Oka Y. Brain Behav, 2012. "The relationship between nerve conduction velocity and fiber morphology during peripheral nerve regeneration."
g = 0.220 .* log10(2*R) + 0.508;
% g = 0.76; % if you want a constant g-ratio

% figure
% plot(2.*R, g)

end

