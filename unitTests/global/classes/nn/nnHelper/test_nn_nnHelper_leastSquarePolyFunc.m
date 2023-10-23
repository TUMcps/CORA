function res = test_nn_nnHelper_leastSquarePolyFunc()
% test_nn_nnHelper_leastSquarePolyFunc - tests nnHelper.leastSquarePolyFunc
%
% Syntax:
%    res = test_nn_nnHelper_leastSquarePolyFunc()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnHelper/leastSquarePolyFunc

% Authors:       Tobias Ladner
% Written:       17-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
rng(1)

for order=1:5
    % exact
    x = rand(1, order+1);
    y = rand(1, order+1);

    coeffs = nnHelper.leastSquarePolyFunc(x, y, order);
    y_p = polyval(coeffs, x);

    res = res && all(withinTol(y, y_p));

    % with noise
    n = 50;

    x = rand(1, order+1);
    y = rand(1, order+1);
    noise = interval(-1, 1) * 0.0001;
    noise = reshape(noise.randPoint(n * (order+1)), n, []);
    y = y + noise;
    y = reshape(y, 1, []);

    x = repmat(x, n, 1);
    x = reshape(x, 1, []);

    coeffs = nnHelper.leastSquarePolyFunc(x, y, order);
    y_p = polyval(coeffs, x);

    res = res && all(withinTol(y, y_p, 1e-3));
end


end

% ------------------------------ END OF CODE ------------------------------
