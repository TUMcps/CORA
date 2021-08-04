function [E] = insc_ellipsoid(Z,comptype)
% insc_ellipsoid - Underapproximates a zonotope by an ellipsoid
%
% Syntax:  
%    E = insc_ellipsoid(Z)
%
% Inputs:
%    Z                  - zonotope object
%    (optional) comptype- specifies whether norm is computed exactly or not
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    Z = zonotope(rand(2,5));
%    E = insc_ellipsoid(Z);
%    plot(Z);
%    hold on
%    plot(E);
%
% References:
%    [1] : M. Cerny, "Goffin’s algorithm for zonotopes" Kybernetika, 
%          vol. 48, no. 5, pp. 890–906, 2012
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: enc_ellipsoid, MVIE, MVEE

% Author:       Victor Gassmann
% Written:      18-September-2019
% Last update:  08-June-2021 (degeneracy check in main file)
% Last revision:---

%------------- BEGIN CODE --------------
rankG = rank(generators(Z));
dimG = dim(Z);
assert(rankG==dimG,'Degeneracy should be handled in main file!');

%check if respective norm should be bounded or computed exactly
if exist('comptype','var') && strcmp(comptype,'exact')
    doExact = true;
else
    doExact = false;
end
%extract zonotope information, i.e. center and generator matrix
G = Z.Z(:,2:end);
c = Z.Z(:,1);
[n,m] = size(G);
%initial guess contained in Z analogous to [1], Sec. 4.1
E0 = G*G';

%cannot directly fit ellipsoid to zonotope
if n<m 
    %transformation matrix T s.t. T*ellipsoid(E0) == unit hyper sphere
    T = inv(sqrtm(E0));
    %transform Z
    Zt = T*Z;
    if doExact
        l = minnorm(Zt);
    else
        error('lower bound not implemented yet!');
    end
    %since transf. ellipsoid Et has radius 1 and Et \in Zt, scaling with l 
    %results in Et touching Zt (for exact norm computation at least at some 
    %point (thus optimal in that sense) (implicit application of inv(T))
    E = ellipsoid(l^2*E0,c);
%can directly compute the optimal E
else
    E = ellipsoid(E0,c);
end
%------------- END OF CODE --------------