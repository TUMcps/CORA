function der = getDerInterval(coeffs, l, u)
% getDerInterval - compute the maximum and the minimum derivative of 
%    the given polynomial
%
% Syntax:
%    int = nnHelper.getDerInterval(coeffs, l, u)
%
% Inputs:
%    coeffs - coefficients of polynomial
%    l - lower bound of input domain
%    u - upper bound of input domain
%
% Outputs:
%    der - interval bounding the derivative within [l,u]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Tobias Ladner
% Written:       13-May-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% find extreme points of derivate of polynomial
p = coeffs;
dp = polyder(p);
dp2 = polyder(dp);
dp2_roots = roots(dp2);
dp2_roots = dp2_roots(imag(dp2_roots) == 0); % filter imaginary roots

% evaluate extreme points of derivative
points = [l, dp2_roots', u];
points = points(l <= points & points <= u);
dp_y = polyval(dp, points);
der = interval(min(dp_y), max(dp_y));

end

%------------- END OF CODE --------------
