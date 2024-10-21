function res = isIntersecting_(fs,S,type,tol,varargin)
% isIntersecting_ - checks if the dimension of the affine hull of a
%    full-dimensional space is equal to the dimension of its ambient space
%    case R^0: only point in vector space is 0 (not representable in
%    MATLAB), so isIntersecting would always return true
%
% Syntax:
%    res = isIntersecting_(fs,S,type,tol)
%
% Inputs:
%    fs - fullspace object
%    S - contSet object or numerical vector
%    type - type of check
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    fs = fullspace(2);
%    p = [1;1]
%    res = isIntersecting(fs,p);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (MW, rename isIntersecting_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[fs,S] = reorderNumeric(fs,S);

if fs.dimension == 0
    throw(CORAerror('CORA:notSupported',...
        'Intersection check of R^0 not supported'));
end

% all singletons intersect the full space
if isnumeric(S)
    res = true;
    return
end

% intersection with empty set is empty
if isa(S,'emptySet')
    res = false;
    return
end

% other set must not be empty
if isa(S,'contSet')
    res = ~representsa_(S,'emptySet',tol);
    return
end

throw(CORAerror('CORA:noops',fs,S));

% ------------------------------ END OF CODE ------------------------------
