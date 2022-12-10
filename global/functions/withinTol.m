function res = withinTol(a,b,varargin)
% withinTol - checks whether two numeric values (scalars, vectors, arrays)
%    are within a given tolerance
%
% Syntax:  
%    res = withinTol(a,b)
%    res = withinTol(a,b,TOL)
%
% Inputs:
%    a,b - double (scalar, vector, matrix)
%    tol - tolerance
%
% Outputs:
%    res - true/false for each comparison
%
% Example: 
%    res = withinTol(1,1+1e-12)
%
% Other m-files required: -
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      19-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = setDefaultValues({1e-8},varargin{:});

% allow scalar values to be expanded
if ~all(size(a)==size(b)) && ~isscalar(a) && ~isscalar(b) && ...
        ~(size(a,1) == size(b,1) && size(a,2) == 1) && ...
        ~(size(a,2) == size(b,2) && size(a,1) == 1) && ...
        ~(size(a,1) == size(b,1) && size(b,2) == 1) && ...
        ~(size(a,2) == size(b,2) && size(b,1) == 1)
    throw(CORAerror('CORA:dimensionMismatch',a,b));
end

if ~isa(a,'double')
    throw(CORAerror('CORA:wrongValue','first','double'));
elseif ~isa(b,'double')
    throw(CORAerror('CORA:wrongValue','second','double'));
end

% absolute tolerance
res_abs = abs(a-b) <= tol;

% relative tolerance
res_rel = abs(a-b) ./ min(abs(a),abs(b)) <= tol;

% joint result
res = res_abs | res_rel;

%------------- END OF CODE --------------