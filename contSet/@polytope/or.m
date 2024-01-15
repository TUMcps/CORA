function P_out = or(P1,P2)
% or - Computes an over-approximation for the union of polytopes
%
% Syntax:
%    P_out = or(P1,P2)
%
% Inputs:
%    P1 - polytope object
%    P2 - polytope object
%
% Outputs:
%    P_out - polytope object enclosing the union
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/or, zonotope/or, polytope/convHull

% Authors:       Niklas Kochdumper, Viktor Kotsev
% Written:       26-November-2019
% Last update:   31-August-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty cases
if representsa_(P1,'emptySet',eps)
    P_out = P2; return
elseif representsa_(P2,'emptySet',eps)
    P_out = P1; return
end

% call convex hull
P_out = convHull(P1,P2);
    
% ------------------------------ END OF CODE ------------------------------
