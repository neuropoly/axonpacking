function [PHI, FR, MVF, AVF] = computeStatistics( R, x, side, g_ratio, resolution)
% stat evaluate from the results of the simulation :
%       - PHI : the extra-axonal density
%       - MVF : the Myelin Volume Fraction
%       - AVF : the Axon Volume Fraction
%       - FR  : the fraction restricted


N = length(R);
pts = reshape(x,2,length(x)/2);
Ls = sqrt(sum(pi*R.^2));

% masks
mask = createCirclesMask(pts',R,side,resolution); mask=double(mask);
Xmin = round((mean(pts(1,:))-2 * Ls/5) * resolution/side); 
Xmax = round((mean(pts(1,:))+2*Ls/5)*resolution/side);
Ymin = round((mean(pts(2,:))-Ls/3)*resolution/side); 
Ymax = round((mean(pts(2,:))+Ls/3)*resolution/side);
mask_trunc = mask(Xmin:Xmax,Ymin:Ymax);

% display masks
% figure 
% subplot(211); imagesc(mask)
% title(['Moyenne des diametres axonal : ',num2str(mean(R(:))),' (variance = 1)'],'FontSize',20,'FontWeight','bold');
% subplot(212); imagesc(mask_trunc)


% PHI
Lx = size(mask_trunc,1);
Ly = size(mask_trunc,2);
PHI = length(find(mask_trunc==1))/(Lx*Ly); % = (fiber + myelin) / background

% FR
FR = g_ratio^2*PHI/(1+(g_ratio^2-1)*PHI); % = fiber / (fiber + background)

% AVF
AVF = g_ratio^2*length(find(mask_trunc==1))/(Lx*Ly); % = fiber / area

% MVF
MVF = (1-g_ratio^2)*length(find(mask_trunc==1))/(Lx*Ly); % = myelin / area

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

