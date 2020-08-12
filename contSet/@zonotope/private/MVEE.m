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
%    Z = zonotope(rand(2,5));
%    E = MVEE(Z);
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
% See also: enc_ellipsoid, inc_ellipsoid

% Author:        Victor Ga�mann
% Written:       18-September-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
c = center(Z);
G = generators(Z);
[n,m] = size(G);
%check if efficient computation is possible
if n==m && rank(G)==n
    E = ellipsoid(n*(G*G'),c);
    return;
end
%ATTENTION: Exact computation, requires vertices -> scales badly
%see [1], Sec. 8.4.1 for more details
%compute zonotope vertices
V = vertices(Z);
N = length(V);

Q = sdpvar(n);
X = V-repmat(c,1,N);
F = cone([ones(size(X,2),1),X'*Q]');
options = sdpsettings;
options.verbose = 0;
%uncomment if sdpt3 is installed
%options.solver = 'sdpt3';
%solve optimization problem
optimize(F, -logdet(Q), options);
E = ellipsoid(inv(value(Q)^2),c);

%------------- END OF CODE --------------