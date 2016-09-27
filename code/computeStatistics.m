function [FVF, FR, MVF, AVF, g] = computeStatistics( D, gap, pts, side, resolution)
% stat evaluate from the results of the simulation :
%       - FVF : the fiber volume fraction e.g. the axon (=disk) density
%       - MVF : the myelin volume fraction
%       - AVF : the axon volume fraction
%       - FR  : the fraction of restricted water

N = length(D);
Ls = sqrt(sum(pi*(D+gap/2).^2))*(4/5);
g =  compute_gratio(D);

% masks
Xmin = max(1,round((mean(pts(1,:)) - Ls/2) * resolution/side));
Xmax = round((mean(pts(1,:)) + Ls/2) * resolution/side);
Ymin = max(1,round((mean(pts(2,:)) - Ls/2) * resolution/side));
Ymax = round((mean(pts(2,:)) + Ls/2) * resolution/side);

% mask1 : axon + myelin mask
mask1 = createCirclesMask(pts', D, side, resolution); mask1=double(mask1);
mask_trunc1 = mask1(Xmin:Xmax,Ymin:Ymax);
area1 = length(find(mask_trunc1==1));
area0 = length(find(mask_trunc1==0));

% mask2 : axon mask (without myelin)
mask2 = createCirclesMask(pts', g.*D, side, resolution); mask2=double(mask2);
mask_trunc2 = mask2(Xmin:Xmax,Ymin:Ymax);
area2 = length(find(mask_trunc2==1));

% background sizes
Lx = size(mask_trunc1,1);
Ly = size(mask_trunc1,2);

% FVF
FVF = area1 / (Lx*Ly); % = (AXON + myelin) / area

% FR
FR = area2 / (area2 + area0); % = AXON / (AXON + background)

% AVF
AVF = area2 / (Lx*Ly); % = AXON / area

% MVF
MVF = (area1 - area2) / (Lx*Ly); % = myelin / area

end

function g = compute_gratio(R)

% Ikeda M, Oka Y. Brain Behav, 2012. "The relationship between nerve conduction velocity and fiber morphology during peripheral nerve regeneration."
g = 0.220 .* log10(2*R) + 0.508;
% g = 0.76; % if you want a constant g-ratio

% figure
% plot(2.*R, g)

end

function mask = createCirclesMask(varargin)
% Create a binary mask from circle centers and radii

centers = varargin{1};
radi = varargin{2};
side = varargin{3};
M = varargin{4};

xc = centers(:,1);
yc = centers(:,2);

[xx,yy] = meshgrid(1:M,1:M);

mask = false(M,M);
for i = 1:numel(radi)
	mask = mask | hypot(xx - xc(i)*M/side, yy - yc(i)*M/side) <= radi(i)*M/side;
end
end

