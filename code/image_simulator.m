
%% Specify simulation inputs

sim_params.N=100;            % number of axons i.e disks to pack  
sim_params.d_mean=4;         % theoretical mean of axon diameters in um
sim_params.d_var=1;         % theoretical variance of axon diameters in um
sim_params.Delta=1;         % gap between the edge of axons in um 
sim_params.iter_max=8000;    % number of iteration i.e migrations to perform. Example: iter_max = 30000 ok if N = 1000

% SECONDARY INPUTS
sim_params.threshold_high=10;     % no diameter above 'threshold_high'
sim_params.threshold_low=0.2;     % no diameter under 'threshold_low'
sim_params.iter_fvf=sim_params.iter_max/10; 

%% Simulate distribution & get axon (AVF) and fiber (FVF) masks

[AVF_mask,FVF_mask]=get_image_masks(sim_params);

%% Specify intensity values for axon, myelin, background

simulation_properties.axon_value=25;
simulation_properties.myelin_value=200;
simulation_properties.background_value=150;

%% Generate image

[img,groundtruth]=generate_img(FVF_mask,AVF_mask,simulation_properties,sim_params);
imshow(img);
imshow(groundtruth);

%% add noise

a1=imnoise(img,'salt & pepper',0.05);
imshow(a1);
imwrite(a1,'a1.tif');

a2=imnoise(img,'speckle',0.05);
imshow(a2);
imwrite(a2,'a2.tif');

h=fspecial('average');
a3=imfilter(img,h);
imshow(a3);
imwrite(a3,'a3.tif');

h=fspecial('gaussian',4,1);
a4=imfilter(img,h);
imshow(a4);
imwrite(a4,'a4.tif');

a5=imadjust(img);
imshow(a5);
imwrite(a5,'a5.tif');

%% add deformation


[x, y] = meshgrid(1:size(img,2), 1:size(img,1));
vx = 0.1*y+0.2*x;   
vy = 0.1*x-0.1*y;   
a6 = interp2(double(img), x-vx, y-vy);
a6=uint8(a6);
imshow(a6);

%% plot histograms

imhist(img);
imhist(a4);


%% synthetic images from axon packing




% 
% mask_stats=compute_stats_from_axonlist(axonlist,PixelSize,img);
% 
% [mean_gap_axon]=gap_axons(axonlist,PixelSize,3);
% good_pixelsize=0.05*1/mask_stats.axon_diam_std;
% 
% 
% 
% 


%% 






