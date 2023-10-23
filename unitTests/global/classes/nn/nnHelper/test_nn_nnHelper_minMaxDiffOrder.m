function res = test_nn_nnHelper_minMaxDiffOrder()
% test_nn_nnHelper_minMaxDiffOrder - tests nnHelper.minMaxDiffOrder
%
% Syntax:
%    res = test_nn_nnHelper_minMaxDiffOrder()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnHelper/minMaxDiffOrder

% Authors:       Tobias Ladner
% Written:       17-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

layer = nnSigmoidLayer();
% compute error
l = 1; u = 2;
[coeffs, ~] = layer.computeApproxPoly(1, 2, 1, 'regression');
[der1l,der1u] = layer.getDerBounds(l, u);
[diff2l,diff2u] = nnHelper.minMaxDiffOrder(coeffs, l, u, layer.f, der1l, der1u);
L = interval(diff2l,diff2u);
% check containment
x = linspace(l, u);
y_p = polyval(coeffs, x);
y_f = layer.f(x);
res = res && all(contains(L, y_f-y_p));

l = 1; u = 2;
[coeffs, ~] = layer.computeApproxPoly(1, 2, 1, 'singh');
[der1l,der1u] = layer.getDerBounds(l, u);
[diff2l,diff2u] = nnHelper.minMaxDiffOrder(coeffs, l, u, layer.f, der1l, der1u);
L = interval(diff2l,diff2u);
% check containment
x = linspace(l, u);
y_p = polyval(coeffs, x);
y_f = layer.f(x);
res = res && all(contains(L, y_f-y_p));

layer = nnTanhLayer();
% compute error
l = 1; u = 2;
[coeffs, ~] = layer.computeApproxPoly(1, 2, 1, 'regression');
[der1l,der1u] = layer.getDerBounds(l, u);
[diff2l,diff2u] = nnHelper.minMaxDiffOrder(coeffs, l, u, layer.f, der1l, der1u);
L = interval(diff2l,diff2u);
% check containment
x = linspace(l, u);
y_p = polyval(coeffs, x);
y_f = layer.f(x);
res = res && all(contains(L, y_f-y_p));

l = 1; u = 2;
[coeffs, ~] = layer.computeApproxPoly(1, 2, 1, 'singh');
[der1l,der1u] = layer.getDerBounds(l, u);
[diff2l,diff2u] = nnHelper.minMaxDiffOrder(coeffs, l, u, layer.f, der1l, der1u);
L = interval(diff2l,diff2u);
% check containment
x = linspace(l, u);
y_p = polyval(coeffs, x);
y_f = layer.f(x);
res = res && all(contains(L, y_f-y_p));


end

% ------------------------------ END OF CODE ------------------------------
