function E = priv_MVIE(Z)
% priv_MVIE - computes the Maximum-Volume-inscribed ellipsoid of a zonotope
%
% Syntax:
%    E = priv_MVIE(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    Z = zonotope([1;0],[2 -1 3; 0 1 2]);
%    E = ellipsoid(Z,'inner:exact');
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
% See also: priv_inscEllipsoid, priv_encEllipsoid

% Authors:       Victor Gassmann
% Written:       18-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to polytope, then to ellipsoid
E = ellipsoid(polytope(Z),'inner');

% ------------------------------ END OF CODE ------------------------------
