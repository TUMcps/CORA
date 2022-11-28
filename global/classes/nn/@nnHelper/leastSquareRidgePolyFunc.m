function coeff = leastSquareRidgePolyFunc(x, y, n, lambda)
% leastSquareRidgePolyFunc - determine the optimal polynomial function fit 
%    that minimizes the squared distance to the data points while using 
%    Tikhonov regularization
%
% Syntax:
%    coeff = nnHelper.leastSquareRidgePolyFunc(x, y, n, lambda)
%
% Inputs:
%    x - x values
%    y - y values
%    n - polynomial order
%    lambda - coefficient of Tikhonov regularization, default to 0.001
%
% Recommended values: -Tanh layer: lambda = 0.002
%                     -Sigmoid layer: lambda = 0.000075
%                     -LeakyReLU layer with any alpha: lambda = 0
%
% Outputs:
%    coeff - coefficients of resulting polynomial
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Sebastian Sigl
% Written:       7-July-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin < 4
    lambda = 0.001;
end

A = x'.^(0:n);
coeff = (A' * A + lambda * eye(n+1)) \ A' * y';

end

%------------- END OF CODE --------------