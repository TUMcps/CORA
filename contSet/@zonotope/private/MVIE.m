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

% Author:       Victor Gaßmann
% Written:      18-September-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
%ATTENTION: Exact computation, requires halfspaces -> scales badly
%see [1], 8.4.2 for more details
n = dim(Z);
%get halfspace representation
P = polytope(Z);
P = get(P,'P');
A = P.A;
b = P.b;
%extract center
c = center(Z);

B = sdpvar(n);
F = cone([b-A*c,A*B]');
options = sdpsettings;
options.verbose = 0;
%uncomment if sdpt3 is installed
%options.solver = 'sdpt3';
%solve optimization problem
optimize(F, -logdet(B), options);
E = ellipsoid(value(B)^2, value(c));

%------------- END OF CODE --------------