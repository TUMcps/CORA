function res = withinTol(a,b,tol)
% withinTol - checks whether two numeric values (scalars, vectors, arrays)
%    are within a given tolerance
%
% Syntax:
%    res = withinTol(a,b)
%    res = withinTol(a,b,tol)
%
% Inputs:
%    a,b - double (scalar, vector, matrix, n-d arrays)
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
%                22-April-2024 (LK, isnumeric check)
%                18-October-2024 (TL, allow n-d arrays)
% Last revision: 20-July-2023 (TL, speed up input parsing)

% ------------------------------ BEGIN CODE -------------------------------

% for speed:
if nargin < 3
    tol = 1e-8;
end

% check dimension
aux_checkDims(a,b)

% direct check for speed reasons
if ~isnumeric(a)
    throw(CORAerror('CORA:wrongValue','first','double'));
elseif ~isnumeric(b)
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

end


% Auxiliary functions -----------------------------------------------------

function aux_checkDims(a,b)
    % check scalar
    if isscalar(a) || isscalar(b)
        % quick exit
        return
    end

    % read size
    size_a = size(a);
    size_b = size(b);
    n = max(numel(size_a),numel(size_b));

    % extend to match size
    size_a = [size_a,ones(1,n-numel(size_a))];
    size_b = [size_b,ones(1,n-numel(size_b))];

    % mismatching dimensions must be scalar in either array
    idxMiss = size_a ~= size_b;
    if ~all(size_a(idxMiss) == 1 | size_b(idxMiss) == 1)
        throw(CORAerror('CORA:dimensionMismatch',a,b));
    end
end

% ------------------------------ END OF CODE ------------------------------
