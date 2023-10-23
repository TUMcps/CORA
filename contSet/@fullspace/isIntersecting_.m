function res = isIntersecting_(fs,S,varargin)
% isIntersecting_ - checks if the dimension of the affine hull of a
%    full-dimensional space is equal to the dimension of its ambient space
%    case R^0: only point in vector space is 0 (not representable in
%    MATLAB), so isIntersecting would always return true
%
% Syntax:
%    res = isIntersecting_(fs,S)
%    res = isIntersecting_(fs,S,type)
%
% Inputs:
%    fs - fullspace object
%    S - contSet object or numerical vector
%    type - (optional) type of check
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

[fs,S] = findClassArg(fs,S,'fullspace');

if fs.dimension == 0
    throw(CORAerror('CORA:notSupported','Intersection check of R^0 not supported'));
end

if isa(S,'emptySet')
    % intersection with empty set is empty
    res = false;

elseif isa(S,'contSet')
    % other set must not be empty
    res = ~representsa_(S,'emptySet',eps);

elseif isnumeric(S)
    % all singletons intersect the full space
    res = true;

end

% ------------------------------ END OF CODE ------------------------------
