function Z = priv_inscZonotope(E,m,mode)
% priv_inscZonotope - inner-approximates a non-degenerate ellipsoid by a
%    zonotope
%
% Syntax:
%    Z = priv_inscZonotope(E,m,mode)
%
% Inputs:
%    E - ellipsoid object
%    m - number of generators
%    mode - computation of zonotope norm (see zonotope/norm_)
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    Z = zonotope(E,'inner:norm',20);
% 
%    figure; hold on;
%    plot(E); plot(Z);
%
% References:
%    [1] P. Leopardi, "Distributing points on the sphere," School of Mathematics
%        and Statistics. PhD thesis. University of New South Wales, 2007
%    [2] The Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox,
%        https://de.mathworks.com/matlabcentral/fileexchange/13356-eq_sphere_partitions
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_encZonotope

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       18-September-2019
% Last update:   27-January-2021 (MA, degenerate case implemented)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out dimension
n = dim(E);
% extract center
c = center(E);
% compute transformation matrix s.t. T*E == unit hypersphere
Tinv = sqrtm(E.Q);
% compute "uniform" distribution of m points on unit hypersphere, see [1],[2]
if n == 1
    G = [-1,1];
else
    G = eq_point_set(n-1,m);
end

% compute zonotope norm
R = norm_(zonotope(zeros(n,1),G),2,mode);
% we want the zonotope to be enclosed in the ellipsoid, so we scale
% zonotope(.,G) such that is barely contained in unit hypersphere,
% and apply inverse transform
Z = zonotope(c,1/R*Tinv*G);

% ------------------------------ END OF CODE ------------------------------
