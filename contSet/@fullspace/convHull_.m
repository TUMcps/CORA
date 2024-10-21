function S_out = convHull_(fs,S,varargin)
% convHull_ - computes the convex hull of a fullspace and another set or a
%    point
%
% Syntax:
%    S_out = convHull_(fs,S)
%
% Inputs:
%    fs - fullspace object
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

% full space is convex
if nargin == 1
    return
end

% ensure that numeric is second input argument
[fs,S] = reorderNumeric(fs,S);

% check dimensions
equalDimCheck(fs,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < fs.precedence
    S_out = convHull(S,fs,varargin{:});
    return
end

% convex hull is full space
S_out = fs;

% ------------------------------ END OF CODE ------------------------------
