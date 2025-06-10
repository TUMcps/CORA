function v = vecleftmtimes(A,v)
% vecleftmtimes - left mtimes for a vectorized matrix
%
% Syntax:
%    v = vecleftmtimes(A,v)
%
% Inputs:
%    A - numeric matrix (n x k)
%    v - vectorized matrix (k*m x h)
%
% Outputs:
%    v - vectorized result
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       16-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init dimensions
[n,k] = size(interval(A));
[km,h] = size(v);
m = km/k;

% reshape vectorized matrix back to matrix
M = reshape(v,k,m,h);

% do matrix multiplication (including broadcasting)
if isa(M,'interval')
    M = interval(A) * M;
else
    M = pagemtimes(A,M);
end

% reshape back to vector
try
    v = reshape(M,n*m,h);
catch ME
    keyboard
end

end

% ------------------------------ END OF CODE ------------------------------
