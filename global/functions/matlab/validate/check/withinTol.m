function res = withinTol(a,b,tol)
% withinTol - checks whether two numeric values (scalars, vectors, arrays)
%    are within a given tolerance
%
% Syntax:
%    res = withinTol(a,b)
%    res = withinTol(a,b,TOL)
%
% Inputs:
%    a,b - double (scalar, vector, matrix)
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false for each comparison
%
% Example: 
%    res = withinTol(1,1+1e-12)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       19-July-2021
% Last update:   03-December-2023 (MW, handling of Inf values)
% Last revision: 20-July-2023 (TL, speed up input parsing)

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 3
    tol = 1e-8;
elseif nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% allow scalar values to be expanded
if ~all(size(a)==size(b)) && ~isscalar(a) && ~isscalar(b) && ...
        ~(size(a,1) == size(b,1) && size(a,2) == 1) && ...
        ~(size(a,2) == size(b,2) && size(a,1) == 1) && ...
        ~(size(a,1) == size(b,1) && size(b,2) == 1) && ...
        ~(size(a,2) == size(b,2) && size(b,1) == 1)
    throw(CORAerror('CORA:dimensionMismatch',a,b));
end

% direct check for speed reasons
if ~isa(a,'double')
    throw(CORAerror('CORA:wrongValue','first','double'));
elseif ~isa(b,'double')
    throw(CORAerror('CORA:wrongValue','second','double'));
elseif ~isscalar(tol) || tol < 0
    throw(CORAerror('CORA:wrongValue','third','nonnegative scalar'));
end

% absolute tolerance
res_abs = abs(a-b) <= tol;

% relative tolerance
res_rel = abs(a-b) ./ min(abs(a),abs(b)) <= tol;

% handling of Inf values
res_inf = isinf(a) & isinf(b) & (sign(a) == sign(b));

% joint result
res = res_abs | res_rel | res_inf;

% ------------------------------ END OF CODE ------------------------------
