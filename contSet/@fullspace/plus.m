function S_out = plus(fs,S)
% plus - overloaded '+' operator for the Minkowski addition of a
%    full-dimensional space and another set or vector
%    case R^0: can only be added to another R^0, resulting in R^0;
%              or to the empty set, resulting in the empty set
%
% Syntax:
%    S_out = fs + S
%    S_out = plus(fs,S)
%
% Inputs:
%    fs - fullspace object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    res - fullspace object
%
% Example: 
%    fs = fullspace(2);
%    Z = zonotope([1;1],[2 1; -1 0]);
%    fs + Z
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call function with lower precedence
if isa(S,'contSet') && S.precedence < fs.precedence
    S_out = S + fs;
    return
end

% ensure that numeric is second input argument
[S_out,~] = reorderNumeric(fs,S);

if S_out.dimension == 0
    throw(CORAerror('CORA:notSupported','Minkowski sum of R^0 not supported'));
end

% check dimensions of ambient space
equalDimCheck(S_out,S);

% empty set case
if representsa_(S,'emptySet',eps)
    S_out = emptySet(dim(fs));
    return
end

% ------------------------------ END OF CODE ------------------------------
