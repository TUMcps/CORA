function int = minMaxQuadFun(a, b, c, l, u)
% minMaxQuadFun - compute the maximimum and minimum value of a quadratic 
%    function on an interval
%
% Syntax:
%    int = nnHelper.minMaxQuadFun(a, b, c, l, u)
%
% Inputs:
%    a-c - polynomial coefficients
%    l - lower bound of input domain
%    u - upper bound of input domain
%
% Outputs:
%    int - interval bounding the min and max value
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

% quadratic function
f = @(x) a * x^2 + b * x + c;

% extremal point of the function
x_ = -b / (2 * a);

% maximum and minimum
if x_ > l && x_ < u
    val = [f(l), f(x_), f(u)];
else
    val = [f(l), f(u)];
end

int = interval(min(val), max(val));
end

%------------- END OF CODE --------------