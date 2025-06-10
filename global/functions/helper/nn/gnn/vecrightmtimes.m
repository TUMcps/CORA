function v = vecrightmtimes(v,A)
% vecrightmtimes - right mtimes for a vectorized matrix
%
% Syntax:
%    v = vecrightmtimes(v,A)
%
% Inputs:
%    v - nectorized matrix (n*k x h)
%    A - numeric matrix (k x m)
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
[nk,h] = size(v);
[k,m] = size(A);
n = nk/k;

% reshape vectorized matrix back to matrix
M = reshape(v,n,k,h);

% do matrix multiplication (including broadcasting)
if isa(M,'interval')
    M = M * A;
else
    M = pagemtimes(M,A);
end

% reshape back to vector
try
    v = reshape(M,n*m,h);
catch ME
    keyboard
end

end

% ------------------------------ END OF CODE ------------------------------
