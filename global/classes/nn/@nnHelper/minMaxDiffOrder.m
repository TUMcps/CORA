function int = minMaxDiffOrder(order, coeffs, l, u, f, der1)
% minMaxDiffOrder - compute the maximum and the minimum difference between the activation
% function and a polynomial fit
%
% Syntax:
%    int = nnHelper.minMaxDiffOrder(order, coeffs, l, u, f, der1)
%
% Inputs:
%    order - order of polynomial
%    coeffs - coefficients of polynomial
%    l - lower bound of input domain
%    u - upper bound of input domain
%    f - function handle of activation function
%    der1 - bounds for derivative of activation functions
%
% Outputs:
%    int - interval bounding the lower and upper error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Tobias Ladner
% Written:       28-March-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

tol = 0.0001;

% calculate bounds for derivative of polynomial

der2 = nnHelper.getDerInterval(coeffs, l, u);
der2 = -der2; % '-' as we calculate f(x) - p(x)

der = supremum(abs(der1-der2));

%------------- END OF CODE --------------

% determine function bounds by sampling
dx = tol / der;
x = linspace(l, u, ceil((u - l)/dx));
x = [l, x, u]; % add l, u in case x is empty (der = 0)
y = f(x) - polyval(coeffs, x);
int = interval(min(y)-tol, max(y)+tol);
end
