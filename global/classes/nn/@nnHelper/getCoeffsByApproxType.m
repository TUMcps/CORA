function coeffs = getCoeffsByApproxType(approx_type, order, l, u, f, dfs)
% getCoeffsByApproxType - finds coefficients using the given approx_type
%
% Syntax:
%    coeffs = nnHelper.getCoeffsByApproxType(approx_type, order, l, u, f, dfs)
%
% Inputs:
%    approx_type - 'regression', 'ridgeregression', 'throw-catch'
%    order - polynomial order
%    l - lower bound
%    u - upper bound
%    f - function handle to approximate
%    dfs - cell array containing derivative function handles
%
% Outputs:
%    coeffs - polynomial coefficients
%             (decending order, e.g. contant term is at coeffs(end))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: respective coeffs calculation

% Author:       Tobias Ladner
% Written:      22-July-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin < 6
    throw(CORAerror('CORA:notEnoughInputArgs',6));
end

% use 10 points within [l, u] for regression
x = linspace(l, u, 10);
y = f(x);

if strcmp(approx_type, "regression")
    % compute polynomial that best fits the activation function
    coeffs = nnHelper.leastSquarePolyFunc(x, y, order);
    coeffs = fliplr(coeffs'); % descending order

elseif strcmp(approx_type, "ridgeregression")
    % compute polynomial that best fits the activation function
    coeffs = nnHelper.leastSquareRidgePolyFunc(x, y, order);
    coeffs = fliplr(coeffs'); % descending order

elseif strcmp(approx_type, "throw-catch")
    if df(l) < df(u)
        coeffs = nnHelper.calcAlternatingDerCoeffs(l, u, order, f, dfs);
    else
        coeffs = nnHelper.calcAlternatingDerCoeffs(u, l, order, f, dfs);
    end
else
    throw(CORAerror('CORA:wrongValue','first',...
        'has to be ''regression'', ''ridgeregression'', or ''throw-catch''.'));
end
coeffs(isnan(coeffs)) = 0;

%------------- END OF CODE --------------