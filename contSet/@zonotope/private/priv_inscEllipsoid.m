function E = priv_inscEllipsoid(Z)
% priv_inscEllipsoid - underapproximates a non-degenerate zonotope by an
%    ellipsoid
%
% Syntax:
%    E = priv_inscEllipsoid(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    Z = zonotope([1;0],[2 -1 3; 0 1 2]);
%    E = ellipsoid(Z,'inner:norm');
%    
%    figure; hold on;
%    plot(Z); plot(E);
%
% References:
%    [1] M. Cerny, "Goffin's algorithm for zonotopes" Kybernetika, vol. 48,
%        no. 5, pp. 890-906, 2012
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_encEllipsoid, priv_MVIE, priv_MVEE

% Authors:       Victor Gassmann
% Written:       18-September-2019
% Last update:   08-June-2021 (degeneracy check in main file)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%extract zonotope information, i.e. center and generator matrix
c = Z.c;
G = Z.G;
[n,m] = size(G);
%initial guess contained in Z analogously to [1], Sec. 4.1
E0 = G*G';

%cannot directly fit ellipsoid to zonotope
if n<m 
    %transformation matrix T s.t. T*ellipsoid(E0) == unit hyper sphere
    T = inv(sqrtm(E0));
    %transform Z
    Zt = T*Z;
    l = minnorm(Zt);
    %since transf. ellipsoid Et has radius 1 and Et \in Zt, scaling with l 
    %results in Et touching Zt (for exact norm computation at least at some 
    %point (thus optimal in that sense) (implicit application of inv(T))
    E = ellipsoid(l^2*E0,c);
%can directly compute the optimal E
else
    E = ellipsoid(E0,c);
end

% ------------------------------ END OF CODE ------------------------------
