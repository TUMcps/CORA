function res = testLong_interval_norm
% testLong_interval_norm - unit test function of norm
%
% Syntax:
%    res = testLong_interval_norm
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       12-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% random cases
res = true;
nrOfTests = 100;

for i=1:nrOfTests

    % random dimension
    n = randi(50);

    % init random interval, where upper bound vertex defines max-norm
    lb = -rand(n,1);
    ub = 1+rand(n,1);
    I = interval(lb,ub);

    % compute norms: 1, 2, Inf
    normI_1 = norm(I,1);
    normI_2 = norm(I,2);
    normI_Inf = norm(I,Inf);

    % check with correct result
    if abs(normI_1 - sum(ub)) > tol || ...
            abs(normI_2 - norm(ub)) > tol || ...
            abs(normI_Inf - max(ub)) > tol
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
