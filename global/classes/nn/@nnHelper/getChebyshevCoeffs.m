function coeff = getChebyshevCoeffs(l, u, order)
% getChebyshevCoeffs - determine the optimal polynomial function fit that 
%    minimizes the squared distance to the data points
%
% Syntax:
%    coeff = nnHelper.getChebyshevCoeffs(l, u, order)
%
% Inputs:
%    l - lower bound
%    u - upper bound
%    order - polynomial order
%
% Outputs:
%    coeff - coefficients of resulting polynomial
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Tobias Ladner
% Written:       22-June-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

b = max(abs([l, u]));

mid = 0;
rad = 8;

syms x
x = (x - mid) / rad;

p = floor((order - 1)/2);
q = ceil((order - 1)/2);
y = cheby(p, q, x);
y = simplify(y);

coeff = coeffs(y, 'All');
coeff = double(coeff);

end

function y = cheby(p, q, x)
y = zeros(size(x), 'like', x);
for i = 1:length(x)
    x_i = x(i);
    y_i = 0;
    for mu = 0:p
        s = nchoosek(mu+q, mu) * ((1 - x_i) / 2)^mu;
        y_i = y_i + s;
    end

    y_i = y_i * ((1 + x_i) / 2)^(q + 1);
    y(i) = y_i;
end
end

%------------- END OF CODE --------------