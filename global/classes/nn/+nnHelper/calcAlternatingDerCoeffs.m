function coeffs = calcAlternatingDerCoeffs(l, u, order, layer)
% calcAlternatingDerCoeffs - calculate coefficients for polynomial using
% 'throw-catch'
%
% Syntax:
%    res = nnHelper.calcAlternatingDerCoeffs(l, u, order, f, dfs)
%
% Inputs:
%    l - lower bound of input domain
%    u - upper bound of input domain
%    order - order of the resulting polynomial
%    layer - nnSShapeLayer
%
% Outputs:
%    res - output of the neural network
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   17-February-2023
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


if layer.df(l) > layer.df(u)
    coeffs = nnHelper.calcAlternatingDerCoeffs(u, l, order, layer);
    return
end

X = [];
Y = [];

exponents = fliplr(0:order);
coeffs_full = ones(size(exponents));
coeffs = ones(size(exponents));

xi = coeffs_full .* (l.^exponents);
yi = [layer.f(l)];

X = [X; xi];
Y = [Y; yi];

exponents = max(0, exponents-1);
coeffs = polyder(coeffs);
coeffs_full = [coeffs, zeros(1, 1+order-length(coeffs))];

xi = coeffs_full .* (l.^exponents);
yi = [layer.df(l)];

X = [X; xi];
Y = [Y; yi];

for i = 2:ceil((order + 1)/2)
    xi = coeffs_full .* (u.^exponents);
    df_i = layer.getDf(i);
    yi = [df_i(u)];

    X = [X; xi];
    Y = [Y; yi];

    exponents = max(0, exponents-1);
    coeffs = polyder(coeffs);
    coeffs_full = [coeffs, zeros(1, 1+order-length(coeffs))];

    xi = coeffs_full .* (l.^exponents);
    yi = [df_i(l)];

    X = [X; xi];
    Y = [Y; yi];
end

% disp([X, Y])

% determine if X is square (and thus also invertible by construction)
if diff(size(X)) == 0
    % X\Y is more efficient than inv(X) * Y
    % https://de.mathworks.com/help/matlab/ref/inv.html#bu6sfy8-1
    coeffs = X \ Y;
else
    coeffs = pinv(X) * Y;
end
coeffs = coeffs';
end

% ------------------------------ END OF CODE ------------------------------
