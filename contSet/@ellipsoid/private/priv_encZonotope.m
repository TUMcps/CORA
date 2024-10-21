function Z = priv_encZonotope(E,nrGen)
% priv_encZonotope - encloses a non-degenerate ellipsoid by a zonotope
%
% Syntax:
%    Z = priv_encZonotope(E,nrGen)
%
% Inputs:
%    E - ellipsoid object
%    nrGen - number of generators of resulting zonotope
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    Z = zonotope(E,'outer:norm',10);
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
% See also: priv_inscZonotope

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       18-September-2019
% Last update:   27-January-2021 (MA, degenerate case implemented)
%                08-June-2021 (VG, moved degeneracy to main file)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out dimension
n = dim(E);
% extract center
c = center(E);
% compute transformation matrix s.t. T*E == unit hyper-sphere
Tinv = sqrtm(E.Q);
% compute "uniform" distribution of m points on unit hyper-sphere, see [1],[2]
if n == 1
    G = [-1,1];
elseif n == nrGen
    Z = enc_parallelotope(E);
    return
else
    G = eq_point_set(n-1,nrGen);
end

L = minnorm(zonotope(zeros(n,1),G));
% we want the ellipsoid to be contained in the zonotope, so we scale
% zonotope(.,G) s.t. it touches E (for exact norm computation), then
% apply retransform
Z = zonotope(c,1/L*Tinv*G);
    
% ------------------------------ END OF CODE ------------------------------
