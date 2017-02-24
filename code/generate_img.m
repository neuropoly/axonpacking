
function [img,groundtruth]=generate_img(FVF_mask,AVF_mask,simulation_properties,sim_params)
for i=1:size(FVF_mask,1)
    test=sum(FVF_mask(i,:));
    if test~=0
        max1=i;
        break;
    end
    
end

for i=size(FVF_mask,1):-1:1
    test=sum(FVF_mask(i,:));
    if test~=0
        max3=i;
        break;
    end
    
end

for j=1:size(FVF_mask,2)
    test=sum(FVF_mask(:,j));
    if test~=0
        max2=j;
        break;
    end
    
end

for j=size(FVF_mask,2):-1:1
    test=sum(FVF_mask(:,j));
    if test~=0
        max4=j;
        break;
    end
    
end

FVF_cropped=FVF_mask(max1-20:max3+20,max2-20:max4+20);
AVF_cropped=AVF_mask(max1-20:max3+20,max2-20:max4+20);

img=AVF_cropped+FVF_cropped;

tmp=img;

img(tmp==2)=simulation_properties.axon_value;
img(tmp==1)=simulation_properties.myelin_value;
img(tmp==0)=simulation_properties.background_value;

img=uint8(img);
groundtruth=AVF_cropped;

sim_params.pixelsize_in_um=0.2;

date_str=datestr(now, 'yyyymmdd-HH-MM-SS');
mkdir([date_str,'_SimulationResults']);
writetable(struct2table(sim_params), [date_str,'_SimulationResults/','Simulation_Parameters.txt']);
imwrite(img,[date_str,'_SimulationResults/','Simulation_img.tif']);
imwrite(groundtruth,[date_str,'_SimulationResults/','Simulation_groundtruth.tif']);

end
