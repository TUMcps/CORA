function res = test_zonotope_intersectStrip_volumeOptimum
% test_zonotope_intersectStrip_volumeOptimum - unit test function of 
% intersectStrip to check whether the optimum volume is found.
% According to [1], the volume of the resulting zonotope is a convex
% function.
%
% Syntax:
%    res = test_zonotope_intersectStrip_volumeOptimum
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035–1043,
%        2005.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       23-December-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Simple 2D example which can be easily visualized
% zonotope
Z = zonotope([-0.0889, -0.5526, -1.2237, -0.0184, -0.1200; ...
   -0.2158, 1.6579, 3.6711, -0.0447, 0.0200]);

% strip
C = [-2, 1];
y_strip = 4.83;
sigma = 0.2;

%% full factorial test of different lambda values
% step size and increments
delta = 0.1;
increment = 2/delta + 1;
% combinations
combs = combinator(increment,2,'p','r');
% lambda matrix
lambdaMat = -1 + (combs-1)*delta;

for i=1:length(lambdaMat(:,1))
    % obtain zonotopes 
    lambda = lambdaMat(i,:)';
    Zres = intersectStrip(Z,C,sigma,y_strip,lambda);
    % compute volume
    vol(i) = volume(Zres);
end

% compute result of optimization 
Zopt = intersectStrip(Z,C,sigma,y_strip,'alamo-volume');

% check if volume of op5timal solution is below brute force method
volOpt = volume(Zopt);
res = volOpt < 1.01*min(vol);

% % create plot
% figure;
% view([124 34]);
% grid('on');
% hold on;
% 
% % create surface
% [X,Y] = meshgrid(-1:delta:1,-1:delta:1);
% Z = reshape(vol,size(X));
% surf(X,Y,Z,'LineStyle','none');


% ------------------------------ END OF CODE ------------------------------
