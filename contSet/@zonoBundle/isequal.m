function res = isequal(zB1,zB2,varargin)
% isequal - checks if two zonotope bundles are equal
%
% Syntax:
%    res = isequal(zB1,zB2)
%    res = isequal(zB1,zB2,tol)
%
% Inputs:
%    zB1 - zonoBundle object
%    zB2 - zonoBundle object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    Z1 = zonotope([1;1], [1 1; -1 1]);
%    Z2 = zonotope([-1;1], [1 0; 0 1]);
%    zB = zonoBundle({Z1,Z2});
%    isequal(zB,zB)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{zB1,'att','zonoBundle'};
                {zB2,'att','zonoBundle'};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% assume false
res = false;

% check number of parallelSets
if zB1.parallelSets ~= zB2.parallelSets
    return
end

% compare individual zonotopes
for z=1:zB1.parallelSets
    if ~isequal(zB1.Z{z},zB2.Z{z})
        return
    end
end

% all comparisons ok
res = true;

% ------------------------------ END OF CODE ------------------------------
