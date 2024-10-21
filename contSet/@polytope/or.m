function P_out = or(P,S,varargin)
% or - computes an over-approximation for the union of polytopes using the
%    convex hull operation
%
% Syntax:
%    P_out = or(P,S)
%    P_out = or(P,S,mode)
%
% Inputs:
%    P - polytope object
%    S - contSet object, numeric
%    mode - 'exact', 'outer', 'inner'
%
% Outputs:
%    P_out - polytope object enclosing the union
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/or, zonotope/or, polytope/convHull_

% Authors:       Niklas Kochdumper, Viktor Kotsev
% Written:       26-November-2019
% Last update:   31-August-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call convex hull
P_out = convHull_(P,S,varargin{:});
    
% ------------------------------ END OF CODE ------------------------------
