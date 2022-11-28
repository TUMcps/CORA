function int = minMaxCubFun(a, b, c, d, l, u)
% minMaxCubFun - compute the maximimum and minimum value of a 
%    cubic function on an interval
%
% Syntax:
%    int = nnHelper.minMaxCubFun(a, b, c, d, l, u)
%
% Inputs:
%    a-d - polynomial coefficients
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

% cubic function
f = @(x) a * x^3 + b * x^2 + c * x + d;

% potential extremal points of the function
x_ = roots([3 * a, 2 * b, c]);

% maximum and minimum
val = [f(l), f(u)];
for i = 1:length(x_)
    if isreal(x_(i)) && x_(i) > l && x_(i) < u
        val = [val, f(x_(i))];
    end
end
int = interval(min(val), max(val));
end

%------------- END OF CODE --------------