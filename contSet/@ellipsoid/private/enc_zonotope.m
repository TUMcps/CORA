function Z = enc_zonotope(E,m,comptype)
% enc_zonotope - overapproximates an ellipsoid by a zonotope
%
% Syntax:  
%    E = enc_zonotope(Z,m,comptype)
%
% Inputs:
%    E       - ellipsoid object
%    comptype- (Optional) Specifies whether function uses a lower bound on the 
%              minimum zonotope norm or the exact value
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid.generateRandom(0,2);
%    Z = enc_zonotope(E,10,'exact');
%    plot(E);
%    hold on
%    plot(Z);
%
% References:
%    [1] : P. Leopardi, "Distributing points on the sphere," School of Mathematics
%          and Statistics. PhD thesis. University of New South Wales, 2007
%    [2] : The Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox,
%          https://de.mathworks.com/matlabcentral/fileexchange/13356-eq_sphere_partitions
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: insc_zonotope

% Author:       Victor Gassmann, Matthias Althoff
% Written:      18-September-2019
% Last update:  27-January-2021 (MA, degenerate case implemented)
%               08-June-2021 (VG, moved degeneracy to main file)
% Last revision:---

%------------- BEGIN CODE --------------
assert(~E.isdegenerate,'Degeneracy should be handled in main file!');

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
else
    G = eq_point_set(n-1,m);
end

if doExact
    L = minnorm(zonotope([zeros(n,1),G]));
else
    error('lower bound on minimum norm not implemented yet!');
end
%we want E \in Z, thus scale zonotope([.,G]) s.t. it touches E (for exact
%norm computation) and apply retransform
Z = zonotope([c,1/L*Tinv*G]);
    
%------------- END OF CODE --------------