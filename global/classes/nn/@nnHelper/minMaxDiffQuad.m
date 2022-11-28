function int = minMaxDiffQuad(a, b, c, l, u, f, df)
% minMaxDiffQuad - compute the maximum and the minimum difference between 
%    the activation function and a quadratic fit
%
% Syntax:
%    int = nnHelper.minMaxDiffQuad(a, b, c, l, u, f, df)
%
% Inputs:
%    a-c - polynomial coefficients
%    l - lower bound of input domain
%    u - upper bound of input domain
%    u - polynomial coefficients
%    f - function handle of activation function
%    df - function handle of derivative
%
% Outputs:
%    int - interval bounding the lower and upper error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 28-March-2022 (TL)

%------------- BEGIN CODE --------------

tol = 0.001;

% compute rough bounds for the derivative
x = linspace(-10, 10, 1000);
y = df(x);
% round to nearest 2nd decimal place.
% TODO can be improved by using the exact bounds [l, u]
der1 = interval(floor(min(y)*100)/100, ceil(100*max(y))/100);

der2 = -(b + 2 * a * interval(l, u)); % quadratic fit

der = supremum(abs(der1-der2));

% determine function bounds by sampling
dx = tol / der;
x = linspace(l, u, ceil((u - l)/dx));
y = f(x) - (a * x.^2 + b * x + c);
int = interval(min(y)-tol, max(y)+tol);
end

%------------- END OF CODE --------------
