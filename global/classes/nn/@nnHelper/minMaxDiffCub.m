function int = minMaxDiffCub(a, b, c, d, l, u, f, df)
% minMaxDiffCub - compute the maximum and the minimum difference between 
%    the activation function and a cubic fit
%
% Syntax:
%    int = nnHelper.minMaxDiffCub(a, b, c, d, l, u, f, df)
%
% Inputs:
%    a-d - polynomial coefficients
%    l - lower bound of input domain
%    u - upper bound of input domain
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
der1 = interval(floor(min(y)*100)/100, ceil(100*max(y))/100);

der2 = nnHelper.minMaxQuadFun(3*a, 2*b, c, l, u); % cubic fit

der = supremum(abs(der1-der2));

% determine function bounds by sampling
dx = tol / der;
x = linspace(l, u, ceil((u - l)/dx));
y = f(x) - a * x.^3 - b * x.^2 - c * x - d;
int = interval(min(y)-tol, max(y)+tol);
end

%------------- END OF CODE --------------