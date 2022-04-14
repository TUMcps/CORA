function P = or(P1,P2)
% or - Computes an over-approximation for the union of polytopes
%
% Syntax:  
%    P = or(P1, P2)
%
% Inputs:
%    P1 - first mptPolytop object
%    P2 - second mptPolytope object
%
% Outputs:
%    P - resulting mptPolytope object enclosing the union
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/or, zonotope/or, mptPolytope/convHull

% Author:       Niklas Kochdumper
% Written:      26-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    P = convHull(P1,P2);
    
%------------- END OF CODE --------------