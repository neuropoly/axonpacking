function  [R,x0] = axonsSetup(axons,side,k)
% author : Tom Mingasson
% axonsSetup creates a randomly generated axon log-normale distribution and initilialize axons positions in a square area defined by 'side'.

% Radii distribution under a log-normale law
R = samplingHastingsMetropolis(axons,k);

% Random positions on a grid for the N axons
N = axons.nbAxons{k};
sqrt_N = round(sqrt(N))+1;
[Xgrid Ygrid] =  meshgrid(side/2-side/4:side/(2*sqrt_N):side/2+side/4, side/2-side/4:side/(2*sqrt_N):side/2+side/4);
Xgrid =Xgrid(:);
Ygrid =Ygrid(:);
Permutations = randperm(sqrt_N^2);
x0 = zeros(N,2);
for i=1:N
    x0(i,:) = [Xgrid(Permutations(i)) Ygrid(Permutations(i))];
end

x0 = reshape(x0',1,2*N)';

end

function R = samplingHastingsMetropolis(axons,k)

sigma_instru = 1;

N = axons.nbAxons{k};
x = zeros(N,1);
x(1)=1;
xmin = 0;
xmax = 15;
axe_hist =linspace(xmin,xmax,100);
delta_hist = axe_hist(2)-axe_hist(1);

for n=1:N-1
    % author : Tom Mingasson
    
    % drawing of x* from x(n) (= current point) under a gaussian intrumental law
    x_star = x(n) + sigma_instru*randn;
    
    % acceptation or reject
    proba_accept = (q_instru(x(n),x_star,sigma_instru)*pobj(x_star,axons,k))/(q_instru(x_star,x(n),sigma_instru)*pobj(x(n),axons,k));
    
    u=rand;
    if u<proba_accept
        x(n+1)=x_star;
    else
        x(n+1)=x(n);
    end
    
    % display
    %     figure(1)
    %     set(figure(1), 'WindowStyle', 'docked') ;
    %     [n_hist,bins]=hist(x(1:n),axe_hist);
    %     h = bar(bins,n_hist/(n*delta_hist));
    %     set(h,'barwidth',0.5)
    %     drawnow
    
end

% hold on
% xth = linspace(0,xmax,1000);
% yth = pobj(xth,param_obj);
% plot(xth,yth,'r*')
% legend('histogramme du tirage','loi log normale theorique')

R =x;

end

function [ pobj ] = pobj(x, axons,k)
% Probability pobj(x) where pobj is the function we want to sample
% Here pobj is a log normale function

meanAxons = axons.meanTheo{k};
varAxons = axons.varTheo{k};
threshold = axons.threshold{k};

% radii distribution
mu = log(meanAxons) - 1/2*log(1 + varAxons/(meanAxons)^2);
sigma = sqrt( log( 1 + varAxons/(meanAxons)^2) );

% troncature
if x>=threshold
    pobj=0;
else
    pobj = lognpdf(x,mu,sigma);
end

end

function [q] = q_instru( x_star, x_courant, sigma_instru )
% Evaluation of the instrumental function q for (x, x_current)

q = 1/(sigma_instru*sqrt(2*pi))*exp(-(x_star-x_courant)^2/(2*sigma_instru^2));
end

