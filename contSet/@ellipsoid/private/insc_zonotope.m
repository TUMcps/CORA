function Z = insc_zonotope(E,m,comptype)
% insc_zonotope - inner-approximates an ellipsoid by a zonotope
%
% Syntax:
%    Z = insc_zonotope(E,m,comptype)
%
% Inputs:
%    E - ellipsoid object
%    m - number of generators
%    comptype - (optional) Specifies whether function uses an upper bound
%               on the maximum zonotope norm or the exact value
%               ('exact'/'approx')
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    Z = zonotope(E,20,'inner:norm');
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
% See also: enc_zonotope

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       18-September-2019
% Last update:   27-January-2021 (MA, degenerate case implemented)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isFullDim(E)
    throw(CORAerror('CORA:degenerateSet','Should be handled in main file'));
end

%check if respective norm should be bounded or computed exactly
if exist('comptype','var') 
    if strcmp(comptype,'exact')
        doExact = true;
    elseif ~strcmp(comptype,'approx')
        throw(CORAerror('CORA:wrongValue','third',"Must be either 'exact' or 'approx'."));
    end
else
    doExact = false;
end 
n = length(E.Q);
%extract center
c = center(E);
%compute transformation matrix s.t. T*E == unit hyper-sphere
Tinv = sqrtm(E.Q);
%compute "uniform" distribution of m points on unit hyper-sphere, see
%[1],[2]
if n==1
    G = [-1,1];
else
    G = eq_point_set(n-1,m);
end

% compute zonotope norm
if doExact
    R = norm_(zonotope([zeros(n,1),G]),2,'exact');
else
    R = norm_(zonotope([zeros(n,1),G]),2,'ub_convex');
end
%we want Z \in E, thus scale zonotope(.,G) such that is barely contained in
%unit hyper-sphere, and apply inverse transform
Z = zonotope([c,1/R*Tinv*G]);

% ------------------------------ END OF CODE ------------------------------
