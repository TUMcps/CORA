function E = MVEE(Z)
% MVEE - Computes the Minimum-Volume-Enclosing ellipsoid (enclosing Z)
%
% Syntax:
%    E = MVEE(Z)
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
%    plot(Z);
%    plot(E);
%
% References:
%    [1] S. Boyd and L. Vandenberghe, Convex optimization. Cambridge
%        university press, 2004
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: enc_ellipsoid, inc_ellipsoid

% Authors:       Victor Gassmann
% Written:       18-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = Z.c;
G = Z.G;
[n,m] = size(G);

%check if efficient computation is possible
if n==m && rank(G)==n
    E = ellipsoid(n*(G*G'),c);
    return;
end

%ATTENTION: Exact computation, requires vertices -> scales badly
%see [1], Sec. 8.4.1 for more details
%compute zonotope vertices
% Bug (?) in vertices: First and last element are the same
V = unique(vertices(Z)','stable','rows')';
assert(all(abs(mean(V,2)-c)<=1e-10),'mean of vertices must result in center');
E = ellipsoid.enclosePoints(V,'min-vol');

% ------------------------------ END OF CODE ------------------------------
