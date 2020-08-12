function [E] = enc_ellipsoid(Z,comptype)
% enc_ellipsoid - Overapproximates a zonotope by an ellipsoid
%
% Syntax:  
%    E = ellipsoid(Z,comptype)
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
%    E = enc_ellipsoid(Z);
%    plot(Z);
%    hold on
%    plot(E);
%
% References:
%    [1] : M. Cerny, "Goffinís algorithm for zonotopes" Kybernetika, 
%          vol. 48, no. 5, pp. 890ñ906, 2012
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: insc_ellipsoid, MVEE, MVIE

% Author:       Victor Gaﬂmann
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
%extract zonotope information, i.e. center and generator matrix
G = Z.Z(:,2:end);
[n,m] = size(G);

c = Z.Z(:,1);

%there are more generators than dimensions: cannot "fit" the ellipse directly to
%generators
if n<m && rank(G)>=n
    %compute initial guess for ellipsoid containing Z ([1], Sec, 4.1) 
    Q0 = m*(G*G');
    TOL = 1e-8;
    Z0 = zonotope([zeros(n,1),G]);
    %compute transformation matrix T s.t. T*ellipsoid(Q0) == unit hyper
    %sphere
    T = inv(sqrtm(Q0));
    %apply same transform to zonotope
    Zt = T*Z0;
    if ~doExact
        lambda = norm(zonotope([zeros(n,1),generators(Zt)]),2,'ub_convex');
    else
        lambda = norm(zonotope([zeros(n,1),generators(Zt)]),2,'exact');
    end
    %since transformed zonotope is still contained in transf. ellipsoid and
    %since radius of said ell. =1, we can choose radius of ell =
    %norm(Zt) and apply inverse transform inv(T) (implicitly done here)
    E = ellipsoid((lambda+TOL)^2*Q0,c);
elseif n==m && rank(G)>=n
%direct mapping possible
    E = ellipsoid(n*(G*G'),c);
else
    error('Not implemented for degenerate zonotopes');
end
%------------- END OF CODE --------------