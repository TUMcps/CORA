function Z = enc_zonotope(E,m,comptype)
% enc_zonotope - over-approximates an ellipsoid by a zonotope
%
% Syntax:
%    Z = enc_zonotope(E,m,comptype)
%
% Inputs:
%    E - ellipsoid object
%    m - ???
%    comptype - (optional) Specifies whether function uses a lower bound on
%               the minimum zonotope norm or the exact value
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    Z = zonotope(E,10,'outer:norm');
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
% See also: insc_zonotope

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       18-September-2019
% Last update:   27-January-2021 (MA, degenerate case implemented)
%                08-June-2021 (VG, moved degeneracy to main file)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isFullDim(E)
    throw(CORAerror('CORA:degenerateSet','Should be handled in main file'));
end

%check if respective norm should be bounded or computed exactly
if exist('comptype','var') && strcmp(comptype,'exact')
    doExact = true;
else
    doExact = false;
end 
%extract center
n = length(E.Q);
c = center(E);
%compute transformation matrix s.t. T*E == unit hyper-sphere
Tinv = sqrtm(E.Q);
%compute "uniform" distribution of m points on unit hyper-sphere, see
%[1],[2]
if n==1
    G = [-1,1];
elseif n == m
    Z = enc_parallelotope(E);
    return
else
    G = eq_point_set(n-1,m);
end

if doExact
    L = minnorm(zonotope([zeros(n,1),G]));
else
    throw(CORAerror('CORA:noops',E,m));
end
%we want E \in Z, thus scale zonotope([.,G]) s.t. it touches E (for exact
%norm computation) and apply retransform
Z = zonotope([c,1/L*Tinv*G]);
    
% ------------------------------ END OF CODE ------------------------------
