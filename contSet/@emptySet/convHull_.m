function S_out = convHull_(O,S,varargin)
% convHull_ - computes the convex hull of an empty set and another set or a
%    point
%
% Syntax:
%    S_out = convHull_(O,S)
%
% Inputs:
%    O - emptySet object
%    S - contSet object or numeric
%
% Outputs:
%    S_out - convex hull
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull

% Authors:       Mark Wetzlinger
% Written:       29-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set is convex
if nargin == 1
    return
end

% ensure that numeric is second input argument
[O,S] = reorderNumeric(O,S);

% check dimensions
equalDimCheck(O,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < O.precedence
    S_out = convHull(S,O,varargin{:});
    return
end

% convex hull is always the other set
if isa(S,'contSet')
    S_out = S.copy();
else
    S_out = S;
end

% ------------------------------ END OF CODE ------------------------------
