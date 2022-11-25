function E = MVIE(Z)
% MVEE - Computes the Maximum-Volume-inscribed ellipsoid (inscribed into Z)
%
% Syntax:  
%    E = MVIE(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    Z = zonotope(rand(2,5));
%    E = MVIE(Z);
%    plot(Z);
%    hold on
%    plot(E);
%
% References:
%    [1] : S. Boyd and L. Vandenberghe, Convex optimization. Cambridge
%          university press, 2004
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: insc_ellipsoid, enc_ellipsoid

% Author:       Victor Gassmann
% Written:      18-September-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

%get halfspace representation
P = polytope(Z);
P = get(P,'P');
A = P.A;
b = P.b;

E = ellipsoid(mptPolytope(A,b),'i');
%------------- END OF CODE --------------