function [coeffs, d] = computeApproxPoly(obj, l, u, varargin)
% computeApproxPoly - computes an approximating polynomial and respective
%   error bound
%
% Syntax:
%    [coeffs, d] = computeApproxPoly(obj, l, u, order, poly_method)
%
% Inputs:
%    obj - nnActivationLayer
%    l,u - bounds to compute the polynomial within
%    order - polynomial order
%    poly_method -  how approximation polynomial is found in nonlinear
%        layers, e.g. 'regression', 'singh', ... 
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Kochdumper, Niklas, et al. "Open-and closed-loop neural network 
%       verification using polynomial zonotopes." 
%       arXiv preprint arXiv:2207.02715 (2022).
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner
% Written:       16-February-2023
% Last update:   ---
% Last revision: 26-May-2023

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 3
    throw(CORAerror('CORA:notEnoughInputArgs', 3))
elseif nargin > 5
    throw(CORAerror('CORA:tooManyInputArgs', 5))
end
[order, poly_method] = setDefaultValues({1, 'regression'}, varargin);

validPolyMethods = {'regression', 'ridgeregression', 'bernstein', ...
    'throw-catch', 'taylor', 'singh'};
inputArgsCheck({ ...
    {l, 'att', 'numeric', 'scalar'}, ...
    {u, 'att', 'numeric', 'scalar'}, ...
    {order, 'att', 'numeric', 'scalar'}, ...
    {poly_method, 'str', validPolyMethods} ...
})

% trivial case
if l == u
    coeffs = [0, obj.f(l)];
    d = 0;
    return
elseif l > u
    throw(CORAerror("CORA:wrongValue", "second/third", 'l <= u'))
end

% init
coeffs = [];
d = [];

% compute approximation polynomial ----------------------------------------

% use at 10 points per coeff within [l, u] for regression
% https://en.wikipedia.org/wiki/One_in_ten_rule
numPoints = 10*(order+1);

if strcmp(poly_method, 'regression')
    x = linspace(l, u, numPoints);
    y = obj.f(x);

    % compute polynomial that best fits the activation function
    coeffs = nnHelper.leastSquarePolyFunc(x, y, order);

elseif strcmp(poly_method, 'ridgeregression')
    x = linspace(l, u, numPoints);
    y = obj.f(x);

    coeffs = nnHelper.leastSquareRidgePolyFunc(x, y, order);

elseif strcmp(poly_method, 'taylor')
    % taylor series expansion at middle point
    c = (l+u)/2;

    % init
    P = [1]; % pascal's triangle
    coeffs = zeros(1,order+1);
    
    % taylor series expansion
    for i=0:order
        df_i = obj.getDf(i);
        coeffs(end-i:end) = coeffs(end-i:end) + ...
            (P .* (-c).^(0:i)) .* df_i(c) ./ factorial(i);

        % prepare for next iteration
        P = [1 P(1:end-1)+P(2:end) 1];
    end

    % TODO lagrange remainder
    % d = ?

elseif strcmp(poly_method, 'bernstein')

    coeffs = nnHelper.findBernsteinPoly(obj.f,l,u,order);

else
    % check custom computation in subclass
    [coeffs, d] = computeApproxPolyCustom(obj, l, u, order, poly_method);
end
   
% parse coeffs and d
if isempty(coeffs)
    % unable to determine coeffs
    throw(CORAerror('CORA:nnLayerNotSupported', obj, ...
        sprintf("'%s' for polynomial of order=%d", poly_method, order)))
    
elseif isempty(d)
    % compute approx error if not found already
    [coeffs, d] = computeApproxError(obj, l, u, coeffs);
end

end

% ------------------------------ END OF CODE ------------------------------
