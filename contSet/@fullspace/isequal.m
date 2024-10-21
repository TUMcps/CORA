function res = isequal(fs,S,varargin)
% isequal - checks if a full-dimensional space is equal to another set or
%    point; case R^0: can only be equal to another R^0
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

narginchk(2,3);

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{fs,'att',{'fullspace','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% ensure that numeric is second input argument
[fs,S] = reorderNumeric(fs,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < fs.precedence
    res = isequal(S,fs,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(fs,S,true)
    res = false;
    return
end

if isa(S,'fullspace')
    res = fs.dimension == S.dimension;
    return
end

if isa(S,'interval')
    % only set that can cover R^n
    res = all(infimum(S) == -Inf) && all(supremum(S) == Inf);
    return
end

% no other sets can cover R^n or represent R^0
if isa(S,'contSet') || isnumeric(S)
    res = false;
    return
end

throw(CORAerror('CORA:noops',fs,S));

% ------------------------------ END OF CODE ------------------------------
