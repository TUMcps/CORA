function res = plus(summand1,summand2)
% plus - overloaded '+' operator for the Minkowski addition of a
%    full-dimensional space and another set or vector
%    case R^0: can only be added to another R^0, resulting in R^0;
%              or to the empty set, resulting in the empty set
%
% Syntax:
%    res = plus(summand1,summand2)
%
% Inputs:
%    summand1 - fullspace object, contSet object, numerical vector
%    summand2 - fullspace object, contSet object, numerical vector
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

% determine fullspace object
[res,summand] = findClassArg(summand1,summand2,'fullspace');

if res.dimension == 0
    throw(CORAerror('CORA:notSupported','Minkowski sum of R^0 not supported'));
end

% check dimensions of ambient space
equalDimCheck(res,summand);

if representsa_(summand,'emptySet',eps)
    % return empty set
    res = emptySet(res.dimension);
else
    % return fs as before
end

% ------------------------------ END OF CODE ------------------------------
