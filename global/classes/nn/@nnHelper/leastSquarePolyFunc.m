function coeff = leastSquarePolyFunc(x, y, n)
% leastSquarePolyFunc - determine the optimal polynomial function fit that 
%    minimizes the squared distance to the data points
%
% Syntax:
%    coeff = nnHelper.leastSquarePolyFunc(x, y, n)
%
% Inputs:
%    x - x values
%    y - y values
%    n - polynomial order
%
% Outputs:
%    coeff - coefficients of resulting polynomial
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

A = x'.^(0:n);
coeff = pinv(A) * y';
end

%------------- END OF CODE --------------