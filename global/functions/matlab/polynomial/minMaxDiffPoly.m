function [diffl,diffu] = minMaxDiffPoly(coeffs1, coeffs2, l, u)
% minMaxDiffPoly - compute the maximum and the minimum difference between 
%     two polynomials on the given domain: min/max_x p_1(x) - p_2(x)
%
% Syntax:
%    L = minMaxDiffPoly(coeffs1, coeffs2, l, u)
%
% Inputs:
%    coeffs1 - coefficients of first polynomial
%    coeffs1 - coefficients of second polynomial
%    l - lower bound of domain
%    u - upper bound of domain
%
% Outputs:
%    L - interval bounding the lower and upper error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       26-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute difference polynomial: p_1(x) - p_2(x)
p = zeros(1,max(length(coeffs1),length(coeffs2)));
p(end-length(coeffs1)+1:end) = coeffs1;
p(end-length(coeffs2)+1:end) = p(end-length(coeffs2)+1:end)-coeffs2;

% determine extreme points
dp = fpolyder(p);
dp_roots = roots(dp);
dp_roots = dp_roots(imag(dp_roots) == 0); % filter imaginary roots
dp_roots = dp_roots(l < dp_roots & dp_roots < u);
extrema = [l, dp_roots', u]; % extrema or boundary
diff = polyval(p, extrema);
            
% compute final approx error
diffl = min(diff);
diffu = max(diff);

end

% ------------------------------ END OF CODE ------------------------------
