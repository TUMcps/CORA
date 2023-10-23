function coeffs = leastSquarePolyFunc(x, y, n)
% leastSquarePolyFunc - determine the optimal polynomial function fit that 
%    minimizes the squared distance to the data points
%
% Syntax:
%    coeffs = nnHelper.leastSquarePolyFunc(x, y, n)
%
% Inputs:
%    x - x values
%    y - y values
%    n - polynomial order
%
% Outputs:
%    coeffs - coefficients of resulting polynomial
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 28-March-2022 (TL)

% ------------------------------ BEGIN CODE -------------------------------

A = x'.^(0:n);
coeffs = pinv(A) * y';
coeffs = fliplr(coeffs');

end

% ------------------------------ END OF CODE ------------------------------
