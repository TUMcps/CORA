function L = minMaxDiffPoly(p1, p2, l, u)
% minMaxDiffPoly - compute the min/max difference between two polynomials within [l, u]
%
% Syntax:
%    int = nnHelper.minMaxDiffPoly(p1, p2, l, u)
%
% Inputs:
%    p1 - coeffs of p1
%    p2 - coeffs of p2
%    l - lower bound of input domain
%    u - upper bound of input domain
%
% Outputs:
%    L - interval bounding the min/max difference
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Tobias Ladner
% Written:       13-June-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% init p = p1-p2
size = max(length(p1), length(p2));
p = zeros(1, size);

% substract polynomials p1 - p2
p(1, 1+size-length(p1):end) = p1;
p(1, 1+size-length(p2):end) = p(1, 1+size-length(p2):end) - p2;

% get extreme points
dp = polyder(p);
x = roots(dp);

x = x(imag(x) == 0); % filter imaginary roots
x = [l, x', u]; % add bounds
x = x(l <= x & x <= u); % filter bounds

% get min/max
y = polyval(p, x);
L = interval(min(y), max(y));

end

%------------- END OF CODE --------------