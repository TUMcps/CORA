function [diffl,diffu] = minMaxDiffOrder(layer, coeffs, l, u)
% minMaxDiffOrder - compute the maximum and the minimum difference between the activation
% function and a polynomial fit
%
% Syntax:
%    L = nnHelper.minMaxDiffOrder(layer, coeffs, l, u)
%
% Inputs:
%    layer - nnActivationLayer
%    coeffs - coefficients of polynomial
%    l - lower bound of input domain
%    u - upper bound of input domain
%
% Outputs:
%    [diffl,diffu] - interval bounding the lower and upper error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnActLayerFromHandle

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   31-August-2022 (adjust tol)
%                30-May-2023 (output bounds)
%                02-May-2025 (added maxPoints)
%                18-December-2025 (refactor,monotonicity)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% settings
tol = 1e-4;
minPoints = 1e4;
maxPoints = 5e9; % requires 40GB

if l == u
    % compute exact result directly
    diff = f(l);
    yp = polyval(coeffs, l);
    diff = diff-yp;
    diffl = diff;
    diffu = diff;
    return
end

% check if monotonicity can be exploited
if ~isempty(layer.monotonicity) && layer.monotonicity >= 2 && numel(coeffs) == 2
    % f is convex/concave, and linear polynomial
    % -> check end points, and extrema
    xs = [l,u];
    ddiff = @(x) layer.df(x) - coeffs(1);
    if ddiff(l) * ddiff(u) > 0
        try
            xs = [xs fzero(ddiff, [l u])];
        catch ME
            keyboard
        end
    end
    ys = layer.f(xs) - polyval(coeffs,xs);
    
    % find bounds
    diffl = min(ys)-eps;
    diffu = max(ys)+eps;
    return
end

% calculate bounds for derivative of f and polynomial
[der1l,der1u] = layer.getDerBounds(l, u);
[der2l,der2u] = nnHelper.getDerInterval(coeffs, l, u);

% der = der1 - -der2; % '-' as we calculate f(x) - p(x)
der = max(abs([ ...
    der1l - -der2l; ...
    der1l - -der2u; ...
    der1u - -der2l; ...
    der1u - -der2u; ...
]));

% determine number of points to sample
dx = tol / der;
reqPoints = ceil((u - l)/dx);
numPoints = min(max(reqPoints, minPoints), maxPoints);

% re-calculate tolerance with number of used points
dx = (u-l)/numPoints;
tol = der * dx;

% sample points
xs = linspace(l, u, numPoints);
xs = [l, xs, u]; % add l, u in case x is empty (der = 0)
diff = layer.f(xs) - polyval(coeffs, xs);

% find bounds
diffl = min(diff)-tol;
diffu = max(diff)+tol;
end

% ------------------------------ END OF CODE ------------------------------
