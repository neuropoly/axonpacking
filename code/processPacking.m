function [finalPositions, overlap, FVF_historic] = processPacking(x0,D,gap,side,ITERmax, ITERfvf)
% author : Tom Mingasson

N = length(D);
FVF_historic = [];
x = x0;

progress = progressBar(ITERmax);

for iter=1:ITERmax
 
    if round(iter/ITERfvf)==iter/ITERfvf | iter==1

        % compute FVF
        [FVF,~,~,~,g_ratio] = computeStatistics(D, gap,  reshape(x,2,length(x)/2), side, 2048);
        FVF_historic = [FVF_historic, FVF];
      
        % display intermediate packing
        set(figure(201), 'Name', 'Disk migration'); clf
        subplot(1,2,1)
        pts = reshape(x,2,length(x)/2);
        t = 0:.1:2*pi+0.1;
        for i=1:N
            patch(D(i)*cos(t) + pts(1,i), D(i)*sin(t) + pts(2,i), 'k', 'EdgeColor', 'k');
            patch(g_ratio(i)*D(i)*cos(t) + pts(1,i), g_ratio(i)*D(i)*sin(t) + pts(2,i), 'w', 'EdgeColor', 'k');
            xlim([0 side])
            ylim([0 side])
            title({['Diameter Mean : ',num2str(mean(D(:))),' µm    ','Diameter Variance : ',num2str(var(D(:))),' µm    ','Gap : ',num2str(gap),' µm    '] ...
                ;['iteration : ',num2str(iter) ]} ,'FontSize',10,'FontWeight','bold');
        end
        axis square

        subplot(1,2,2)
        plot([1:length(FVF_historic)]*ITERfvf, FVF_historic, 'r*-')
        legend('FVF')
        drawnow
        
    end
    
    MyGrad = computeGrad(x, D+gap/2, side);
    x = x + MyGrad;
    
    progress.update();
    
end

progress.close();

finalPositions = reshape(x,2,length(x)/2);

% evaluate overlap in the packing
overlap = 0;
for i = 1:N-1
    for j = i+1:N
        [~,~,area] = areaIntersect(finalPositions(:,i),D(i),finalPositions(:,j),D(j));
        overlap = overlap + area;
    end
end
disp(['overlap in the packing: ', num2str(overlap)])

end

function MyGrad = computeGrad(x, D, side)
% author : Tom Mingasson

Kcenter0 = 0.01;     % center step coeff for disks withOUT overlapping
Kcenter1 = 0;        % center step coeff for disks with overlapping
Krep = 0.1;          % repulsion step coeff for disks with overlapping

pts = reshape(x,2,length(x)/2);
N=size(pts,2);

% intersection
P = squareform(pdist(pts','euclidean'))+eye(N);
Rsum = (repmat(D,1,N)' + repmat(D,1,N)).*(tril(ones(N,N),-1)+triu(ones(N,N),1));
L = (Rsum./P-1); % >0 if intersection
Lbin = L; Lbin(Lbin>0) = 1; Lbin(Lbin<=0) = 0;  F = floor(linspace(1,N+1,2*N+1)); F=F(1:end-1); Lbin2 = Lbin(:,F); 
inter1_index = repmat(sum(Lbin),2,1); inter1_index(inter1_index>0)=1;   % disks that overlap
inter0_index = 1 - inter1_index;                                        % disks that NOT overlap

% attraction
pts_centered = side/2-pts;
attraction_norm = sqrt(pts_centered(1,:).^2+pts_centered(2,:).^2);
attraction = pts_centered./repmat(attraction_norm,[2,1]);

% repulsion
U = repmat(x',N,1)-repmat(pts',1,N);
Usum  = sum(U.*Lbin2,1); 
Unorm = sqrt(Usum(1:2:end).^2 + Usum(2:2:end).^2); Unorm(Unorm==0) = 1;
Unormalization = repmat(Unorm,2,1); Unormalization = Unormalization(:)'; 
Usum_normed = Usum./Unormalization ;
repulsion = (Usum_normed.*inter1_index(:)')';

MyGrad = reshape(Kcenter0.*attraction.*inter0_index + Kcenter1.*attraction.*inter1_index,1,2*N)' + Krep.*repulsion;

end


function [p,q,areatotal,numberOfPoints] = areaIntersect(C1,r1,C2,r2)
% Circle centers C1 and C2 with radius r1 and r2.

% Compute x and y
% Solve in easy coordinates
a = norm(C1-C2);
% Check for divide by zero
if isequal(a,0)
    % Circles have same center. Return the area of the smaller circle.
    p = [NaN;NaN]; q = [NaN;NaN]; numberOfPoints = Inf;
    areatotal = pi*min(r1,r2)^2;
    return
end
x = 0.5*(a + (r1^2 - r2^2)/a);
if r1^2 < x^2 % Check for sqrt of negative
    p = [NaN;NaN]; q = [NaN;NaN];
    if r1 + r2 < a
        areatotal = 0; numberOfPoints = 0;
    else % One circle is inside the other
        areatotal = pi*min(r1,r2)^2;
        numberOfPoints = Inf;
    end
    return
end
y = sqrt(r1^2 - x^2);

% Original coordinates basis
i = (C2-C1)/norm(C2-C1);
j = null(i');

% Intersection points in original coordinates
p = C1 + i*x + j*y ;
q = C1 + i*x - j*y;

% Compute the angle theta between radius and x-axis
% in the easy coordinates
theta1 =  atan2(y,x);
theta2 =  atan2(y,a-x);

% Obtain A1 with x, y, r1
area1 = theta1*r1^2 - x*y;
% Obtain A2 with a-x, y, r2
area2 = theta2*r2^2 - (a-x)*y;

areatotal = area1 + area2;

if isequal(p,q)
    numberOfPoints = 1;
else
    numberOfPoints = 2;
end
end
