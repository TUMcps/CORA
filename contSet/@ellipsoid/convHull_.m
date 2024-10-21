function S_out = convHull_(E,S,varargin)
% convHull_ - computes an overapproximation of the convex hull of an
%    ellipsoid and another set representation
%
% Syntax:
%    E = convHull_(E,S)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet class object (or class array)
%
% Outputs:
%    S_out - convex hull
%
% Example: 
%    E1 = ellipsoid([3 -1; -1 1],[1;0]);
%    E2 = ellipsoid([5 1; 1 2],[1;-1]);
%    E = convHull(E1,E2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull, ellipsoid/or

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   19-March-2021 (complete rewrite)
%                04-July-2022 (input checks)
% Last revision: 29-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% check inputs
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {S,'att',{'contSet','numeric'}}});

% ellipsoid is already convex
if nargin == 1
    S_out = E;
    return
end

% ensure that numeric is second input argument
[E,S] = reorderNumeric(E,S);

% check dimensions
equalDimCheck(E,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < E.precedence
    S_out = convHull(S,E,varargin{:});
    return
end

% convex hull with empty set
if representsa_(S,'emptySet',eps)
    return;
end

% simply call 'or'
if isa(S,'ellipsoid')
    S_out = or(E,S,'outer');
    return
end

% no other cases implemented
throw(CORAerror('CORA:noops',E,S));

% ------------------------------ END OF CODE ------------------------------
