function S_out = convHull(ls,S,varargin)
% convHull - computes an enclosure for the convex hull of a level set and
%    another set or a point
%
% Syntax:
%    S = convHull(ls,S)
%
% Inputs:
%    ls - levelSet object
%    S - contSet object or numeric
%
% Outputs:
%    S_out - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       29-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convex hull of set itself
if nargin == 1
    throw(CORAerror("CORA:noops",ls));
end

% check dimensions
equalDimCheck(ls,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < ls.precedence
    S_out = convHull(S,ls,varargin{:});
    return
end

% nothing supported
throw(CORAerror("CORA:noops",ls,S));

% ------------------------------ END OF CODE ------------------------------
