function E = priv_MVEE(Z)
% priv_MVEE - computes the Minimum-Volume-Enclosing ellipsoid (enclosing Z)
%
% Syntax:
%    E = priv_MVEE(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    Z = zonotope([1;0],[2 -1 3; 0 1 2]);
%    E = ellipsoid(Z,'outer:exact');
%    
%    figure; hold on;
%    plot(Z); plot(E);
%
% References:
%    [1] S. Boyd and L. Vandenberghe, Convex optimization. Cambridge
%        university press, 2004
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_encEllipsoid, priv_incEllipsoid

% Authors:       Victor Gassmann
% Written:       18-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[n,nrGen] = size(Z.G);

% check if efficient computation is possible
if n == nrGen && rank(Z.G) == n
    E = ellipsoid(n*(Z.G*Z.G'),Z.c);
    return;
end

% ATTENTION: Exact computation, requires vertices -> scales badly
% see [1], Sec. 8.4.1 for more details
% compute zonotope vertices
% Bug (?) in vertices: First and last element are the same
V = unique(vertices(Z)','stable','rows')';
assert(all(abs(mean(V,2)-Z.c) <= 1e-10),...
    'mean of vertices must result in center');
E = ellipsoid.enclosePoints(V,'min-vol');

% ------------------------------ END OF CODE ------------------------------
