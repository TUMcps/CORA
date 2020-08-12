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

% Author:       Victor Gassmann
% Written:      18-September-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
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
T = inv(sqrtm(E.Q));
%compute "uniform" distribution of m points on unit hyper-sphere, see
%[1],[2]
G = eq_point_set(n-1,m);
%make sure rank(G)=n
counter = 1;
while rank(G)<n
    G = eq_point_set(n-1,m+counter);
    counter = counter + 1;
end
if counter > 1
    disp(['Increased m from ',num2str(m),' to ' , num2str(m+counter), ' since uniform equal partition algorithm returned rank-deficient matrix for your specification']);
end
if doExact
    L = minnorm(zonotope([zeros(n,1),G]));
else
    error('lower bound on minimum norm not implemented yet!');
end
%we want E \in Z, thus scale zonotope([.,G]) s.t. it touches E (for exact
%norm computation) and apply retransform
Z = zonotope([c,1/L*inv(T)*G]);
%------------- END OF CODE --------------