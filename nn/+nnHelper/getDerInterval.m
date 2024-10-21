function [derl,deru] = getDerInterval(coeffs, l, u)
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
%    derl,deru - interval bounding the derivative within [l,u]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       13-May-2022
% Last update:   30-May-2023 (TL, fpolyder, output bounds)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% find extreme points of derivate of polynomial
p = coeffs;
dp = fpolyder(p);
dp2 = fpolyder(dp);
dp2_roots = roots(dp2);
dp2_roots = dp2_roots(imag(dp2_roots) == 0); % filter imaginary roots

% evaluate extreme points of derivative
points = [l, dp2_roots', u];
points = points(l <= points & points <= u);
dp_y = polyval(dp, points);

derl = min(dp_y);
deru = max(dp_y);

end

% ------------------------------ END OF CODE ------------------------------
