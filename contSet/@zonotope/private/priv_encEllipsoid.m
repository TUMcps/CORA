function E = priv_encEllipsoid(Z,mode)
% priv_encEllipsoid - encloses a non-degenerate zonotope by an ellipsoid
%
% Syntax:
%    E = priv_encEllipsoid(Z,mode)
%
% Inputs:
%    Z - zonotope object
%    mode - norm computation of zonotope: 'exact', 'ub_convex'
%
% Outputs:
%    E - ellipsoid object
%
% Example:
%    Z = zonotope([1;0],[2 -1 3; 0 1 2]);
%    E = ellipsoid(Z,'outer:norm_bnd');
%    
%    figure; hold on;
%    plot(Z);
%    plot(E);
%
% References:
%    [1] M. Cerny, "Goffin's algorithm for zonotopes" Kybernetika, vol. 48,
%        no. 5, pp. 890--906, 2012
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_inscEllipsoid, priv_MVEE, priv_MVIE

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       18-September-2019
% Last update:   22-December-2020 (MA, degenerate zonotopes added)
%                08-June-2021 (VG, moved degeneracy to main file)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% extract zonotope information, i.e. center and generator matrix
c = Z.c;
G = Z.G;
[n,m] = size(G);

% there are more generators than dimensions: cannot "fit" the ellipse
% directly to generators
if n<m 
    % compute initial guess for ellipsoid containing Z ([1], Sec. 4.1) 
    Q0 = m*(G*G');
    TOL = ellipsoid.empty(1).TOL;
    Z0 = zonotope([zeros(n,1),G]);
    %compute transformation matrix T s.t. T*ellipsoid(Q0) == unit hyper
    %sphere
    T = inv(sqrtm(Q0));
    %apply same transform to zonotope
    Zt = T*Z0;
    lambda = norm_(zonotope(zeros(n,1),generators(Zt)),2,mode);    
    %since transformed zonotope is still contained in transf. ellipsoid and
    %since radius of said ell. =1, we can choose radius of ell =
    %norm(Zt) and apply inverse transform inv(T) (implicitly done here)
    E = ellipsoid((lambda+TOL)^2*Q0,c);
else
    % direct mapping possible
    E = ellipsoid(n*(G*G'),c);
end

% ------------------------------ END OF CODE ------------------------------
