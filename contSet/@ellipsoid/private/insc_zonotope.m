function Z = insc_zonotope(E,m,comptype)
% insc_zonotope - underapproximates an ellipsoid by a zonotope
%
% Syntax:  
%    E = insc_zonotope(Z,N)
%
% Inputs:
%    E       - ellipsoid object
%    comptype- (Optional) Specifies whether function uses an upper bound on the 
%              maximum zonotope norm or the exact value
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid.generateRandom(0,2);
%    Z = insc_zonotope(E,20);
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
% See also: enc_zonotope

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
n = length(E.Q);
%extract center
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
    R = norm(zonotope([zeros(n,1),G]),2,'exact');
else
    R = norm(zonotope([zeros(n,1),G]),2,'ub_convex');
end
%we want Z \in E, thus scale zonotope(.,G) such that is barely contained in
%unit hyper-sphere, and apply inverse transform
Z = zonotope([c,1/R*inv(T)*G]);
%------------- END OF CODE --------------