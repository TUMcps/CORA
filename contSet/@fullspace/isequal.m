function res = isequal(fs,S,varargin)
% isequal - checks if two full-dimensional spaces are equal
%    case R^0: can only be equal to another R^0
%
% Syntax:
%    res = isequal(fs,S)
%    res = isequal(fs,S,tol)
%
% Inputs:
%    fs - fullspace object
%    S - contSet object or numerical vector
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    fs1 = fullspace(2);
%    fs2 = fullspace(3);
%    res1 = isequal(fs1,fs1);
%    res2 = isequal(fs1,fs2);
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

if isa(S,'fullspace')
    res = fs.dimension == S.dimension;

elseif isa(S,'interval')
    % only set that can cover R^n
    res = all(infimum(S) == -Inf) && all(supremum(S) == Inf);

elseif isa(S,'contSet') || isnumeric(S)
    % no other sets can cover R^n or represent R^0
    res = false;

else
    % unknown path...
    throw(CORAerror('CORA:noops',fs,S));
end

% ------------------------------ END OF CODE ------------------------------
