function coeff = calcTaylorTanh(order)
% calcTaylorTanh - find polynomial coefficients for tanh using taylor series
%
% Syntax:
%    coeff = nnHelper.calcTaylorTanh(order)
%
% Inputs:
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
% Written:       27-June-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

coeff = zeros(1, order+1);
for n = 1:order
    if 2 * n - 1 > order
        break
    end
    coeff((2 * n - 1)+1) = bernoulli(2*n) * 4^n * (4^n - 1) / factorial(2*n);
end

end

%------------- END OF CODE --------------